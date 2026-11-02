# Benchmark responses take unit-cube inputs; remaining columns are irrelevant.
borehole <- function(x) {
  lower <- c(.05, 100, 63070, 990, 63.1, 700, 1120, 9855)
  upper <- c(.15, 50000, 115600, 1110, 116, 820, 1680, 12045)
  q <- sweep(sweep(x[, 1:8, drop = FALSE], 2, upper - lower, "*"),
             2, lower, "+")
  2 * pi * q[, 3] * (q[, 4] - q[, 6]) /
    (log(q[, 2] / q[, 1]) *
       (1 + 2 * q[, 7] * q[, 3] /
          (log(q[, 2] / q[, 1]) * q[, 1]^2 * q[, 8]) + q[, 3] / q[, 5]))
}

# Positive-sign, two-active-input Michalewicz function, with x_physical = pi*x.
michalewicz <- function(x) {
  sin(pi * x[, 1]) * sin(pi * x[, 1]^2)^20 +
    sin(pi * x[, 2]) * sin(2 * pi * x[, 2]^2)^20
}

# Returns training x/y and independent test newx/noiseless truth.
prepare_borehole <- function(n = 50, d = 10, m = 100) {
  x <- lhs::maximinLHS(n, d, dup = 1, optimize.on = "grid")
  lower <- apply(x, 2, min)
  width <- apply(x, 2, max) - lower
  x <- sweep(sweep(x, 2, lower, "-"), 2, width, "/")
  y <- borehole(x) + rnorm(n, sd = .02)
  # The paper rescales training columns only; test inputs remain on [0,1]^d.
  newx <- lhs::randomLHS(m, d)
  list(x = x, y = y, newx = newx, truth = borehole(newx))
}

prepare_michalewicz <- function(n = 60, d = 6, m = 100) {
  x <- lhs::maximinLHS(n, d, dup = 1, optimize.on = "grid")
  y <- michalewicz(x) + runif(n, -.02, .02)
  newx <- lhs::randomLHS(m, d)
  list(x = x, y = y, newx = newx, truth = michalewicz(newx))
}

# Original-response-scale errors, using the paper's test-SST denominator.
example_metrics <- function(prediction, truth) {
  c(RMSE = sqrt(mean((prediction - truth)^2)),
    SRMSE = sqrt(sum((prediction - truth)^2) /
                   sum((truth - mean(truth))^2)))
}
