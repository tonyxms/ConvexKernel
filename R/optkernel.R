# Ordered search dictionary and the original initialization pool.
.optkernel_theta <- c(
  as.vector(outer(c(1, 3, 5, 7, 9), 10^(-2:3))), 10000, 30000)
.optkernel_initial_theta <- .optkernel_theta[seq_len(25L)]

.optkernel_prepare <- function(x, x.norm) {
  cache_mb <- 128
  n <- nrow(x)
  p <- ncol(x)
  if (x.norm) {
    x <- scale(x, center = TRUE, scale = TRUE)
    x.mean <- attr(x, "scaled:center")
    x.scale <- attr(x, "scaled:scale")
  } else {
    x.mean <- rep(0, p)
    x.scale <- rep(1, p)
  }
  # Preserve the original pair order and arithmetic, including rowSums order.
  D <- matrix(0, nrow = n * (n - 1) / 2, ncol = p)
  ll <- 0
  for (k in 1:(n-1)) {
    ll <- tail(ll, 1) + (1:(n-k))
    D[ll, ] <- (x[rep(k, n-k), ] - x[(k+1):n, ])^2
  }
  cache <- new.env(parent = emptyenv())
  cache$distance <- new.env(hash = TRUE, parent = emptyenv())
  cache$kernel <- new.env(hash = TRUE, parent = emptyenv())
  cache$distance_bytes <- 0
  cache$kernel_bytes <- 0
  cache$distance_limit <- cache_mb * 1024^2 / 4
  cache$kernel_limit <- cache_mb * 1024^2 * 3 / 4
  lower <- lower.tri(matrix(0, n, n))
  # Same column-major pair order as D and the cached lower triangle.
  pair_index <- which(lower, arr.ind = TRUE)
  list(n = n, p = p, D = D,
       x.mean = x.mean, x.scale = x.scale,
       lower = lower, pair_row = pair_index[, 1L],
       pair_col = pair_index[, 2L], cache = cache)
}

.optkernel_key <- function(act_dim, theta = NULL) {
  key <- paste(which(act_dim), collapse = ",")
  if (!is.null(theta)) key <- paste(key, sprintf("%.17g", theta), sep = ":")
  key
}

.optkernel_packed <- function(prepared, act_dim, theta) {
  cache <- prepared$cache
  key <- .optkernel_key(act_dim, theta)
  if (exists(key, cache$kernel, inherits = FALSE)) return(cache$kernel[[key]])
  dkey <- .optkernel_key(act_dim)
  if (exists(dkey, cache$distance, inherits = FALSE)) {
    distance <- cache$distance[[dkey]]
  } else {
    if (sum(act_dim) == 1) distance <- drop(prepared$D[, act_dim])
    else distance <- rowSums(prepared$D[, act_dim])
    bytes <- as.numeric(object.size(distance))
    if (cache$distance_bytes + bytes <= cache$distance_limit) {
      cache$distance[[dkey]] <- distance
      cache$distance_bytes <- cache$distance_bytes + bytes
    }
  }
  packed <- exp(-drop(theta * distance))
  bytes <- as.numeric(object.size(packed))
  # Bounded first-admission cache: no thrashing if the candidate bank is large.
  if (cache$kernel_bytes + bytes <= cache$kernel_limit) {
    cache$kernel[[key]] <- packed
    cache$kernel_bytes <- cache$kernel_bytes + bytes
  }
  packed
}

.optkernel_warm <- function(prepared, inputs, theta) {
  # Warm in the parent, so forked workers can read/reuse cached immutable data.
  # Once full, leave the remaining candidates lazy instead of computing twice.
  bytes <- as.numeric(object.size(numeric(nrow(prepared$D))))
  cache <- prepared$cache
  for (i in seq_along(theta)) {
    if (cache$kernel_bytes + bytes > cache$kernel_limit) break
    .optkernel_packed(prepared, inputs[i, ], theta[i])
  }
  invisible(NULL)
}

.optkernel_matrix <- function(prepared, act_dim, theta) {
  K <- matrix(0, nrow = prepared$n, ncol = prepared$n)
  K[prepared$lower] <- .optkernel_packed(prepared, act_dim, theta)
  K <- t(K) + K
  diag(K) <- rep(1, prepared$n)
  K
}

.optkernel_score_state <- function(prepared, K, alpha) {
  alpha <- drop(alpha)
  # Rebuilt once per selection round, after the current kernel/alpha is known.
  # Keep the actual diagonal: weighted kernels need not have a bitwise unit diag.
  list(pair_weight = alpha[prepared$pair_row] * alpha[prepared$pair_col],
       current = K[prepared$lower],
       diagonal = sum(alpha^2 * (1 - diag(K))))
}

.optkernel_phi_packed <- function(prepared, act_dim, theta, score_state, mu) {
  packed <- .optkernel_packed(prepared, act_dim, theta)
  # Symmetry: a' (K_new-K_now) a = diagonal + 2 * lower-triangle sum.
  # Fuse only the inner product. Keep diagonal correction and scaling in R.
  inner <- .optkernel_dot_difference_cpp(
    score_state$pair_weight, packed, score_state$current)
  drop(-mu * (score_state$diagonal + 2 * inner))
}

.optkernel_state <- function(K, y, mu) {
  A <- K
  diag(A) <- diag(A) + mu
  R <- chol(A)
  alpha <- backsolve(R, forwardsolve(t(R), y))
  # Match the original n x 1 coefficient shape.
  alpha <- matrix(alpha, ncol = 1L)
  list(R = R, alpha = alpha, loss = drop(mu * crossprod(y, alpha)))
}

.optkernel_precision <- function(state) {
  # Exact dense inverse, once per fitted mu; no repeated factorization.
  chol2inv(state$R)
}

.optkernel_smoother <- function(Q, mu) {
  S <- -mu * Q
  diag(S) <- 1 + diag(S)
  S
}

# Fit a convex combination of squared-exponential kernels to matrix x and response y.
# mu is eta; drop.tol is DEL; stop.tol/stop2.tol control weight/search convergence.
# Returns an opt.kernel list with aligned supports, weights, coefficients, K and S.
# .prepared and .return_state reuse distances and Cholesky factors during LOOCV.
optkernel3_fun <- function(x, y, mu, theta=.optkernel_theta, max_active=4,
    heredity="strong", drop.tol=0.005, drop.tol2=0.05, stop.tol=0.005,
    stop2.tol=0.005, max_iter_1=1000, max_iter_2=1000, x.norm=FALSE,
    .prepared=NULL, .return_state=FALSE, final.prune=FALSE, final.drop.tol=0.02) {
  y <- drop(y)
  p <- ncol(x)
  if (is.null(.prepared)) .prepared <- .optkernel_prepare(x, x.norm)
  x.mean <- .prepared$x.mean
  x.scale <- .prepared$x.scale
  K.matrix <- function(act_dim, theta) .optkernel_matrix(.prepared, act_dim, theta)

  # The same directional derivative, evaluated directly on the cached triangle.
  phi_func <- function(act_dim, theta, score_state) {
    .optkernel_phi_packed(.prepared, act_dim, theta, score_state, mu)
  }
  d_func <- function(K, alpha) drop(crossprod(alpha, K %*% alpha))

  # Multiplicative simplex update, initialized uniformly for each support set.
  cal_weight=function(K_sel){
    if(length(dim(K_sel))==2){
      lambda=1
      return(list(weight=lambda, K=K_sel, state=.optkernel_state(K_sel,y,mu)))
    }
    nn=dim(K_sel)[3]
    lambda=rep(1/nn,nn)
    K_now=0
    l=1 # The iteration for algorithm 2
    for (i in 1:nn){
      K_now=K_now+lambda[i]*K_sel[ , ,i]
    }
    state <- .optkernel_state(K_now,y,mu)
    while (l<=max_iter_2){
      lambup=rep(0,nn)
      alpha <- state$alpha
      if (numCores == 1L) {
        for (i in 1:nn) {
          lambup[i] <- lambda[i]*d_func(K_sel[ , ,i],alpha)
        }
      } else {
        lambup=foreach(i=1:nn, .combine=c,.export="d_func") %dopar% {
          lambda[i]*d_func(K_sel[ , ,i],alpha)
        }
      }
      den=sum(lambup)
      lambda=lambup/den
      K_now=0
      for (i in 1:nn){
        K_now=K_now+lambda[i]*K_sel[ , ,i]
      }
      next_state <- .optkernel_state(K_now,y,mu)
      converged <- abs((state$loss-next_state$loss)/state$loss)<=stop.tol
      state <- next_state
      if (converged) break
      l=l+1
    }
    return(list(weight=lambda, K=K_now, state=state))
  }

  numCores <- getOption("optkernel.workers", 1L)
  if (numCores > 1L) {
    registerDoParallel(numCores)
    on.exit(stopImplicitCluster(), add=TRUE)
  }
  # Start with one univariate kernel, using the frozen 25-value sampling pool.
  kernel_input=matrix(FALSE,nrow=p,ncol=p)
  diag(kernel_input)=TRUE
  C_size=p
  B_size=length(theta)
  kernel_input_all=kernel_input[rep(1:C_size,each=B_size), ]
  theta_all=rep(theta,times=C_size)
  left=1:(C_size*B_size)
  int_k=sample(p,size=1)
  int_kernel=rep(FALSE,p)
  int_kernel[int_k]=TRUE
  pick=logical(C_size*B_size)
  for (i in left){
    pick[i]=all(kernel_input_all[i, ]==int_kernel)
  }
  original_candidate <- sample(x=((int_k-1L)*25L+
    seq_along(.optkernel_initial_theta)), size=1)
  initial_theta <- .optkernel_initial_theta[((original_candidate-1L) %% 25L)+1L]
  pick_int <- left[pick][theta_all[left[pick]] == initial_theta]
  left=setdiff(left,pick_int)
  sel=pick_int
  K_now=K.matrix(act_dim=kernel_input_all[sel, ],theta=theta_all[sel])
  state <- .optkernel_state(K_now,y,mu)
  K_sel=K_now
  size_sel=1
  weight=1
  t=0 # Search iterations; terminal reoptimization is excluded.

  K_order=1
  sel_kernel_input=NULL
  sel_theta=NULL

  while (K_order<=max_active){
    .optkernel_warm(.prepared, kernel_input_all, theta_all)
    t0=0 # The iteration for algorithm 1
    while (t0<=max_iter_1 && length(left)>0){
      t=t+1
      t0=t0+1
      possible_phi=rep(0,length(left))
      alpha <- state$alpha
      score_state <- .optkernel_score_state(.prepared, K_now, alpha)
      if (numCores == 1L) {
        for (i in 1:length(left)) {
          possible_phi[i] <- phi_func(kernel_input_all[left[i], ],theta_all[left[i]],score_state)
        }
      } else {
        possible_phi=foreach(i=1:length(left), .combine=c) %dopar% {
          phi_func(kernel_input_all[left[i], ],theta_all[left[i]],score_state)
        }
      }
      better_ix=which.min(possible_phi)
      if(possible_phi[better_ix]>=0) break
      add=left[better_ix]
      sel=c(sel,add)
      K_add=K.matrix(act_dim=kernel_input_all[add, ],theta=theta_all[add])
      K_sel=abind(K_sel,K_add,along=3)
      weight_old=weight
      fitted <- cal_weight(K_sel)
      weight <- fitted$weight
      # Reject an added kernel below DEL and restore the previous fit.
      if (tail(weight,1)<drop.tol){
        if (length(sel)==1) {
          sel=integer(0)
        } else {
          sel=sel[1:length(sel)-1]
        }
        weight=weight_old
        K_sel=K_sel[ , ,1:(dim(K_sel)[3]-1)]
        break
      }
      else{
        left=left[-better_ix]
        size_sel=size_sel+1
      }
      # cal_weight already formed this exact weighted sum in the original order.
      old_loss <- state$loss
      K_now <- fitted$K
      state <- fitted$state
      if (abs((old_loss-state$loss)/old_loss)<=stop2.tol) break
    }
    # Delete small weights and renormalize the retained support.
    drop=(weight<drop.tol)
    if(any(drop)){
      keep=(!drop)
      drop.w=weight[drop]
      weight.scale=1-sum(drop.w)
      weight=weight[keep]
      weight=weight/weight.scale
      K_sel=K_sel[ , ,keep]
      if(K_order>1){
        sel_kernel_input=sel_kernel_input[keep[1:size_sel_previous_stages], ]
        sel_theta=sel_theta[keep[1:size_sel_previous_stages]]
        sel=sel[keep[(size_sel_previous_stages+1):size_sel]]
        if (size_sel_previous_stages==1){
          sel_kernel_input=t(as.matrix(sel_kernel_input))
        }
        if (size_sel_previous_stages==0){
          sel_kernel_input=NULL
          sel_theta=NULL
        }
        size_sel_previous_stages=length(sel_theta)
      }
      else{
        sel=sel[keep]
      }
      size_sel=size_sel-length(drop.w)
      K_now=0
      if (length(weight)==1 && length(dim(K_sel))==2) {
        K_now=K_now+weight*K_sel
      } else{
        for (i in 1:length(weight)){
          K_now=K_now+weight[i]*K_sel[ , ,i]
        }
      }
      state <- .optkernel_state(K_now,y,mu)
    }
    # Update active inputs from the retained supports.
    if (K_order==1){
      act_dim=integer(0)
      for (j in 1:length(sel)){
        next_temp=(1:p)[kernel_input_all[sel[j], ]]
        act_dim=c(act_dim,next_temp)
      }
      act_dim=unique(act_dim)
      size_act_dim=length(act_dim)
    }
    if (length(sel)!=0) {
      sel_kernel_input=rbind(sel_kernel_input,kernel_input_all[sel, ])
      sel_theta=c(sel_theta,theta_all[sel])
      size_sel_previous_stages=length(sel_theta)
    }
    if (K_order==2){
      act_dim=integer(0)
      for (j in 1:length(sel_theta)){
        next_temp=(1:p)[sel_kernel_input[j, ]]
        act_dim=c(act_dim,next_temp)
      }
      act_dim=unique(act_dim)
      size_act_dim=length(act_dim)
    }
    # Apply the chosen heredity rule when lower-order active inputs shrink.
    if (K_order>2){
      act_dim_old=act_dim
      act_dim=integer(0)
      for (j in 1:nrow(sel_kernel_input)) {
        if (sum(sel_kernel_input[j, ])>K_order-1) break
        next_temp=(1:p)[sel_kernel_input[j, ]]
        act_dim=c(act_dim,next_temp)
      }
      if (length(act_dim)!=0) {

        act_dim=unique(act_dim)
        size_act_dim=length(act_dim)
        if(length(act_dim)<length(act_dim_old)){
          keep=rep(TRUE,nrow(sel_kernel_input))
          ss=setdiff(act_dim_old,act_dim)
          if (heredity=="strong") {
            for (j in 1:nrow(sel_kernel_input)) {
              if (any(sel_kernel_input[j,ss]) && weight[j]<drop.tol2){
                keep[j]=FALSE
              }
            }
          }
          if (heredity=="weak") {
            for (j in 1:nrow(sel_kernel_input)) {
              if ((!any(sel_kernel_input[j,act_dim])) && weight[j]<drop.tol2){
                keep[j]=FALSE
              }
            }
          }
          drop=(!keep)
          if(any(drop)){
            drop.w=weight[drop]
            weight.scale=1-sum(drop.w)
            weight=weight[keep]
            weight=weight/weight.scale
            sel_kernel_input=sel_kernel_input[keep, ]
            sel_theta=sel_theta[keep]
            size_sel_previous_stages=length(sel_theta)
            K_sel=K_sel[ , ,keep]
            size_sel=size_sel-length(drop.w)
          }
          K_now=0
          if (length(weight)==1 && length(dim(K_sel))==2) {
            K_now=K_now+weight*K_sel
          } else{
            for (i in 1:length(weight)){
              K_now=K_now+weight[i]*K_sel[ , ,i]
            }
          }
          state <- .optkernel_state(K_now,y,mu)
        }
      }
      act_dim=integer(0)
      for (j in 1:length(sel_theta)){
        next_temp=(1:p)[sel_kernel_input[j, ]]
        act_dim=c(act_dim,next_temp)
      }
      act_dim=unique(act_dim)
      size_act_dim=length(act_dim)
    }

    # Preserve combination order for the chosen heredity candidate space.
    K_order=K_order+1
    if(K_order>max_active || K_order>p) break
    if (heredity=="strong") {
      if(K_order>size_act_dim) break
      kernel_input=matrix(FALSE, nrow=choose(size_act_dim,K_order),ncol=p)
      K_act=t(combn(act_dim,K_order))
      for (i in 1:nrow(K_act)){
        kernel_input[i,K_act[i, ]]=TRUE
      }
    }
    if (heredity=="weak") {
      kernel_input=NULL
      K_act=t(combn(1:p,K_order))
      for (i in 1:nrow(K_act)){
        if(length(intersect(K_act[i, ],act_dim))>0) {
          next_temp=rep(FALSE,p)
          next_temp[K_act[i, ]]=TRUE
          kernel_input=rbind(kernel_input,next_temp)
        }
      }
      row.names(kernel_input)=NULL
    }

    C_size=nrow(kernel_input)
    kernel_input_all=kernel_input[rep(1:C_size,each=B_size), ]
    theta_all=rep(theta,times=C_size)
    left=1:(C_size*B_size)
    sel=integer(0)
  }
  # Optional one-pass terminal pruning precedes conditional LOOCV scoring.
  if (final.prune) {
    terminal <- .optkernel_terminal_reopt(weight, K_sel, sel_kernel_input,
      sel_theta, K_now, state, final.drop.tol, cal_weight)
    weight <- terminal$weight
    K_sel <- terminal$K_sel
    sel_kernel_input <- terminal$kernel_input
    sel_theta <- terminal$theta
    K_now <- terminal$K
    state <- terminal$state
    final_pruning <- terminal$final_pruning
  }
  coeff <- state$alpha
  S <- if (.return_state) NULL else .optkernel_smoother(.optkernel_precision(state),mu)
  optK <- list(S=S, K=K_now, coeff=coeff, weight=weight, K_sel=K_sel,
    kernel_input=sel_kernel_input, theta=sel_theta, heredity=heredity, Iter=t,
    dropweight=drop.tol, dropactdim=drop.tol2, mu=mu,
    normx=list(x.norm=x.norm, x.mean=x.mean, x.scale=x.scale))
  if (final.prune) optK$final_pruning <- final_pruning
  class(optK)="opt.kernel"
  if (.return_state) return(list(model=optK, state=state))
  return(optK)
}
