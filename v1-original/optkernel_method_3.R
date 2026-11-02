optkernel3_fun=function(x,y,mu,theta=seq(from=0,to=10,length=50),max_active=2,heredity=c("strong","weak"),drop.tol=0.01,drop.tol2=0.05,stop.tol=0.1,stop2.tol=0.1,max_iter_1=100,max_iter_2=100,x.norm=TRUE)
### x is the training data and y is the test data
### mu is the nugget effect (which is \eta in the paper)
### theta is the pre-defined set for the values of the kernel parameter
### max_active is the maximum dimensions of input variables (which is MaxDim in the paper)
### heredity depends on how we construct the basic kernel space
### drop.tol is the parameter DEL for Algorithm 1 in the paper
### drop.tol2 is used to check the the validation of sel_kernel_input and set new act_dim
### stop.tol is the parameter Tol for Algorithm 2 in the paper
### stop2.tol is the parameter Tol for Algorithm 1 in the paper
### max_iter_1 is the parameter MaxIter for Algorithm 1 in the paper (need check)
### max_iter_2 is the parameter MaxIter_0 for Algorithm 2 in the paper
### x.norm is the parameter whether to standardize the input data or not
  
{
  call <- match.call()
  
  
  if (!is.matrix(x)) x <- as.matrix(x)
  y <- drop(y)
  n <- nrow(x)	
  p <- ncol(x)
  one <- drop(rep(1,n)) # The vector of 1's
  I <- diag(n) # The identity matrix
  
  ### Standardize x if asked. 
  if (x.norm) {
    x <- scale(x,center=T,scale=T)
    x.mean <- attributes(x)$'scaled:center'
    x.scale <- attributes(x)$'scaled:scale'
  }
  else {
    x.mean<-rep(0,p)
    x.scale<-rep(1,p)	     	
  }
  
  ### Compute the Distance matrix 
  mzmax <- n*(n-1)/2
  D <- matrix(0,nrow=mzmax,ncol=p)
  ll <- 0
  for (k in 1:(n-1)){
    ll <- tail(ll,1)+(1:(n-k))
    D[ll,] <- (x[rep(k,n-k),]-x[(k+1):n,])^2
  }
  
  ###Compute the kernel matrix
  K.matrix <- function(act_dim,theta){
    # ix <- (1:p)[act_dim]
    if (sum(act_dim)==1) 
      Dis2 <- drop(theta*D[,act_dim])
    else 
      Dis2 <- drop(theta*rowSums(D[,act_dim]))
    Dis <- exp(-Dis2)
    K <- matrix(0,nrow=n,ncol=n)
    K[lower.tri(K)] <- Dis
    K <- t(K)+K
    diag(K) <- one
    return(K)
  }
  
  ###Compute the loss function
  loss=function(K){
    c=solve(K+mu*I)%*%y
    ll=t(y-K%*%c)%*%(y-K%*%c)+mu*t(c)%*%K%*%c
    return(drop(ll))
  }
  
  ###Compute the phi function of K_new and K_old which is the directional derivative in the paper
  phi_func=function(K_new,K_old,q){
    ll=mu*t(y)%*%(q%*%(K_new-K_old)%*%q)%*%y
    return(drop(-ll))
  }
  
  ###Compute the d_delta function which is the formula 10 in the paper
  d_func=function(K,q){
    ll=t(y)%*%(q%*%(K)%*%q)%*%y
    return(drop(ll))
  }
  
  ###update weight (find the stationary point)
  cal_weight=function(K_sel,weight_old){
    if(length(dim(K_sel))==2){
      lambda=1
      return(lambda)
    }
    nn=dim(K_sel)[3]
    # Set the initial value for weight
    lambda=rep(1/nn,nn)
    #lambda=c(weight_old*((nn-1)/nn),1/nn)
    K_now=0
    l=1 # The iteration for algorithm 2
    for (i in 1:nn){
      K_now=K_now+lambda[i]*K_sel[ , ,i]
    }
    while (l<=max_iter_2){
      lambup=rep(0,nn)
      q=solve(K_now+mu*I)
      lambup=foreach(i=1:nn, .combine=c,.export="d_func") %dopar% {
        lambda[i]*d_func(K_sel[ , ,i],q)
      }
      den=sum(lambup)
      lambda=lambup/den
      K_old=K_now
      K_now=0
      for (i in 1:nn){
        K_now=K_now+lambda[i]*K_sel[ , ,i]
      }
      if (abs((loss(K_old)-loss(K_now))/loss(K_old))<=stop.tol) break
      l=l+1
    }
    return(lambda)
  }
  
  
  numCores=detectCores()
  registerDoParallel(numCores)
  ### Starting the 1-dim kernel
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
  pick_int=sample(x=left[pick],size=1)
  left=setdiff(left,pick_int)
  sel=pick_int
  K_now=K.matrix(act_dim=kernel_input_all[sel, ],theta=theta_all[sel])
  K_sel=K_now # The kernel matrix we have selected
  size_sel=1
  weight=1
  t=0 # The total iteration for the method
  
  K_order=1
  sel_kernel_input=NULL
  sel_theta=NULL
  
  while (K_order<=max_active){
    t0=0 # The iteration for algorithm 1
    while (t0<=max_iter_1 && length(left)>0){
      t=t+1
      t0=t0+1
      size_left=length(left)
      possible_phi=rep(0,length(left))
      q=solve(K_now+mu*I)
      possible_phi=foreach(i=1:length(left), .combine=c) %dopar% {
        phi_func(K.matrix(act_dim=kernel_input_all[left[i], ],theta=theta_all[left[i]]),K_now,q)
      }
      better_ix=which.min(possible_phi)
      if(possible_phi[better_ix]>=0) break
      add=left[better_ix]
      sel=c(sel,add)
      K_add=K.matrix(act_dim=kernel_input_all[add, ],theta=theta_all[add]) # The kernel we want to add
      K_sel=abind(K_sel,K_add,along=3)
      weight_old=weight
      weight=cal_weight(K_sel,weight_old)
      ### If the weight for the added kernel is too small, we drop it and end the iteration.
      if (tail(weight,1)<drop.tol){
        if (length(sel)==1) {
          sel=integer(0)
        } else {
          sel=sel[1:length(sel)-1]
        }
        weight=weight_old
        K_sel=K_sel[ , ,1:(dim(K_sel)[3]-1)]
        #cat("drop tail","\n")
        break
      }
      else{
        left=left[-better_ix]
        size_sel=size_sel+1
      }
      K_old=K_now
      K_now=0
      if (length(weight)==1 && length(dim(K_sel))==2) {
        K_now=K_now+weight*K_sel
      } else{
        for (i in 1:length(weight)){
          K_now=K_now+weight[i]*K_sel[ , ,i]
        } 
      }
      if (abs((loss(K_old)-loss(K_now))/loss(K_old))<=stop2.tol) break
    }
    # remove the support kernels whose weight is smaller than drop.tol (DEL)
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
      K_old=K_now
      K_now=0
      if (length(weight)==1 && length(dim(K_sel))==2) {
        K_now=K_now+weight*K_sel
      } else{
        for (i in 1:length(weight)){
          K_now=K_now+weight[i]*K_sel[ , ,i]
        } 
      }
    }
    ### Get the active dimensions
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
    ### check the validation of sel_kernel_input and set new act_dim
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
          if(heredity=="strong"){
            for (j in 1:nrow(sel_kernel_input)) {
              if (any(sel_kernel_input[j,ss]) && weight[j]<drop.tol2){
                keep[j]=FALSE
              }
            }
          }
          if(heredity=="weak"){
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
          K_old=K_now
          K_now=0
          if (length(weight)==1 && length(dim(K_sel))==2) {
            K_now=K_now+weight*K_sel
          } else{
            for (i in 1:length(weight)){
              K_now=K_now+weight[i]*K_sel[ , ,i]
            } 
          }
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
    
    ### Construct the higher order candidate kernel space depends on the heredity we choose
    K_order=K_order+1
    if(K_order>max_active) break
    if(heredity=="strong" && K_order>size_act_dim) break
    if(heredity=="strong" && K_order<=size_act_dim){
      kernel_input=matrix(FALSE, nrow=choose(size_act_dim,K_order),ncol=p)
      K_act=t(combn(act_dim,K_order))
      for (i in 1:nrow(K_act)){
        kernel_input[i,K_act[i, ]]=TRUE
      }
    }
    if (heredity=="weak"){
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
  stopImplicitCluster()
  coeff=solve(K_now+mu*I)%*%y
  S=K_now%*%solve(K_now+mu*I)
  optK=list(S=S,K=K_now,coeff=coeff,weight=weight,K_sel=K_sel,kernel_input=sel_kernel_input,theta=sel_theta,heredity=heredity,Iter=t,dropweight=drop.tol,dropactdim=drop.tol2,mu=mu,normx=list(x.norm=x.norm,x.mean=x.mean,x.scale=x.scale))
  class(optK)="opt.kernel"
  return(optK)
}