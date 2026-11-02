# Fit each eta independently; score the fixed learned kernel with exact LOO.
# x: n-by-p matrix; y: length n; mu_seeds: one predetermined seed per eta.
# Returns the winning opt.kernel list, without a second fit.
optkernel_fun_loocv <- function(x, y, mu = c(.005, .01, .02, .05, .1, .5),
    theta = .optkernel_theta, max_active = 4, heredity = "strong",
    drop.tol = .005, drop.tol2 = .05, stop.tol = .005, stop2.tol = .005,
    max_iter_1 = 1000, max_iter_2 = 1000, x.norm = FALSE,
    diagnostics = TRUE, final.prune = FALSE, final.drop.tol = .02,
    mu_seeds) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  prepared <- .optkernel_prepare(x, x.norm)
  scores <- numeric(length(mu))
  best_score <- Inf
  best_model <- NULL
  best_Q <- NULL
  best_index <- NA_integer_
  for (i in seq_along(mu)) {
    set.seed(mu_seeds[i])
    fit <- optkernel3_fun(x, y, mu = mu[i], theta = theta,
      max_active = max_active, heredity = heredity, drop.tol = drop.tol,
      drop.tol2 = drop.tol2, stop.tol = stop.tol, stop2.tol = stop2.tol,
      max_iter_1 = max_iter_1, max_iter_2 = max_iter_2, x.norm = x.norm,
      .prepared = prepared, .return_state = TRUE,
      final.prune = final.prune, final.drop.tol = final.drop.tol)
    Q <- .optkernel_precision(fit$state)
    residual <- drop(fit$model$coeff) / diag(Q)
    scores[i] <- mean(residual^2)
    # A strict comparison retains the first eta in an exact tie.
    if (scores[i] < best_score) {
      best_score <- scores[i]
      best_model <- fit$model
      best_Q <- Q
      best_index <- i
    }
    rm(fit, Q)
  }
  best_model$S <- .optkernel_smoother(best_Q, mu[best_index])
  if (diagnostics) {
    attr(best_model, "loocv") <- data.frame(mu = mu, score = scores,
      selected = seq_along(mu) == best_index)
    attr(best_model, "cache") <- c(
      distance_bytes = prepared$cache$distance_bytes,
      kernel_bytes = prepared$cache$kernel_bytes)
  }
  best_model
}
