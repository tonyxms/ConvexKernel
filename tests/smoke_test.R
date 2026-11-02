# Run in a clean session from the repository root.
source("R/load_all.R")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(701)
x <- matrix(runif(48), nrow = 16, ncol = 3)
y <- sin(2 * x[, 1]) + x[, 2]^2
newx <- matrix(runif(21), nrow = 7, ncol = 3)
model <- optkernel_fun_loocv(x, y, mu_seeds = 3101:3106)
prediction <- optkernel_pred(x, model, newx)
L <- length(model$weight)
loo <- attr(model, "loocv")

stopifnot(inherits(model, "opt.kernel"), is.list(model),
  identical(dim(model$coeff), c(16L, 1L)),
  identical(dim(model$K), c(16L, 16L)),
  identical(dim(model$S), c(16L, 16L)),
  identical(dim(model$kernel_input), c(as.integer(L), 3L)),
  is.logical(model$kernel_input), length(model$theta) == L,
  all(is.finite(model$weight)), all(model$weight >= 0),
  abs(sum(model$weight) - 1) < 1e-10,
  is.numeric(prediction), is.null(dim(prediction)), length(prediction) == 7L,
  all(is.finite(prediction)), nrow(loo) == 6L, all(is.finite(loo$score)),
  sum(loo$selected) == 1L, model$mu == loo$mu[which.min(loo$score)])

# Verify the retained kernel, solve, and LOO score independently of their helpers.
K <- matrix(0, 16, 16)
for (j in seq_len(L)) {
  Kj <- if (L == 1L && length(dim(model$K_sel)) == 2L)
    model$K_sel else model$K_sel[, , j]
  K <- K + model$weight[j] * Kj
}
Q <- solve(K + diag(model$mu, 16))
score <- mean((drop(Q %*% y) / diag(Q))^2)
stopifnot(max(abs(K - model$K)) < 1e-10,
  max(abs(model$coeff - Q %*% y)) < 1e-8,
  max(abs(model$S %*% y - K %*% model$coeff)) < 1e-8,
  abs(score - loo$score[loo$selected]) < 1e-10,
  identical(optkernel_pred(x, model, newx[1, , drop = FALSE]), prediction[1]))

# Weak heredity can introduce a variable through an interaction.
x_weak <- as.matrix(expand.grid(x1 = c(0, 1), x2 = c(0, 1), x3 = (0:3)/3))
y_weak <- (2*x_weak[, 1] - 1) * (1 + .75*(2*x_weak[, 2] - 1))
set.seed(1)
strong <- optkernel3_fun(x_weak, y_weak, mu = .05)
set.seed(1)
weak <- optkernel3_fun(x_weak, y_weak, mu = .05, heredity = "weak")
# No interactions exist beyond the number of input variables.
set.seed(1)
weak_bounded <- optkernel3_fun(x_weak, y_weak, mu = .05,
  heredity = "weak", max_active = ncol(x_weak))
weak_prediction <- optkernel_pred(x_weak, weak, x_weak)
stopifnot(!any(strong$kernel_input[, 2]),
  any(weak$kernel_input[, 1] & weak$kernel_input[, 2]),
  identical(weak$heredity, "weak"), identical(weak, weak_bounded),
  length(weak_prediction) == nrow(x_weak), all(is.finite(weak_prediction)),
  all(weak$weight >= 0), abs(sum(weak$weight) - 1) < 1e-10,
  max(abs(weak_prediction - drop(weak$S %*% y_weak))) < 1e-8)

cat(sprintf("Smoke test passed: n=16, p=3, predictions=7, eta=%g, supports=%d\n",
            model$mu, L))
cat("Weak heredity checks passed: interaction selection and dimension limit.\n")
