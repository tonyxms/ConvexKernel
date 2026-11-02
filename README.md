# ConvexKernel

R implementation of convex kernel learning for prediction and input selection.
The method combines Gaussian kernels using strong or weak heredity and selects the regularization parameter μ (`mu` in the code; η in the paper) by leave-one-out (LOO) scoring.

The original GitHub release is preserved unchanged in [`v1-original/`](v1-original/).
The repository root contains the revised implementation.

## Setup and run

Tested with R 4.5.3 on macOS. 
A C++ toolchain for R is required: Command Line Tools on macOS, Rtools on Windows, or C++ build tools on Linux.
Install the R dependencies once (`lhs` is used only by the examples):

```r
install.packages(c("Rcpp", "abind", "foreach", "doParallel", "lhs"))
```

Run from the repository root:

```sh
Rscript --vanilla examples/borehole.R
Rscript --vanilla examples/michalewicz.R
Rscript --vanilla tests/smoke_test.R
```

The examples generate their own data and print the six LOO scores, selected η, RMSE, SRMSE, and selected input columns.

## Fit and predict

```r
source("R/load_all.R")
source("R/benchmarks.R")
set.seed(101)
data <- prepare_borehole()
model <- optkernel_fun_loocv(data$x, data$y, mu_seeds = 1101:1106)
prediction <- optkernel_pred(data$x, model, data$newx)
attr(model, "loocv")
```

Input data must contain finite numeric values, and `mu` must be positive.

- `x`: n × p training matrix (`n >= 2`, `p >= 2`).
- `y`: length-n response vector or n × 1 matrix.
- `newx`: new inputs with the same variables, column order, units, and preprocessing as `x`.

`x.norm = FALSE` is the default and paper setting. 
Setting it to `TRUE` standardizes inputs using training means and sample standard deviations.
The response is not rescaled.

`mu_seeds` supplies one random seed per candidate η. 
The LOO routine fits each candidate and returns the model with the lowest LOO score. 
Scores are computed with each learned kernel fixed.
Use `optkernel3_fun(x, y, mu, ...)` to fit a single η. 
Optionally call `set.seed()` beforehand for reproducibility.

## Returned model

Both fitters return a **list of class `opt.kernel`**. 
Main fields are listed below, where L is the number of retained kernels.

| Field | Meaning |
|---|---|
| `weight`, `theta` | Kernel weights and Gaussian kernel parameters |
| `kernel_input` | L × p logical matrix identifying each kernel's variables |
| `K_sel` | Individual kernels: n × n × L array, or n × n matrix when L = 1 |
| `K` | Weighted training kernel matrix, n × n |
| `S` | `S = K (K + mu I)^(-1)`, n × n; `S %*% y` gives fitted values at the training inputs |
| `coeff` | Prediction coefficients `(K + mu I)^(-1) y`, n × 1 |
| `mu` | Fitted regularization parameter |

The model also records search settings, input scaling, and optional pruning details. 
The LOO fitter provides its score table through `attr(model, "loocv")` when `diagnostics = TRUE` (default). 
Prediction requires the original training matrix in its original row order and returns one value per row of `newx`.

## Examples and paper settings

| Example | Training samples | Input variables | Test samples |
|---|---:|---:|---:|
| Borehole | 50 | 10 | 100 |
| Michalewicz | 60 | 6 | 100 |

Each example uses one replication, inputs in [0,1], and all six η candidates `(0.005, 0.01, 0.02, 0.05, 0.1, 0.5)` with the full kernel search.
Borehole has eight active inputs. 
Michalewicz uses the positive-sign function with two active inputs. 
Data generation follows the paper's conventions.
SRMSE is defined on the test set as

$$
\mathrm{SRMSE}
= \sqrt{\frac{\sum_{i=1}^{m}(\hat y_i-y_i)^2}
                  {\sum_{i=1}^{m}(y_i-\bar y)^2}}.
$$

Here, $m$ is the number of test samples, $y_i$ is the true response, and $\hat y_i$ is its prediction.

The paper uses 100 replications, larger datasets, and additional scenarios.

The defaults retain the Section 6 settings. 
The paper uses `heredity = "strong"` (default). 
`heredity = "weak"` is also available.
The three method variants use the same fitting function:

| Variant | Arguments to `optkernel_fun_loocv()` |
|---|---|
| optK 0.005 | Defaults |
| optK-pruned | `final.prune = TRUE` (terminal threshold 0.02) |
| optK 0.05 | `drop.tol = 0.05` |

## Files

```text
v1-original/             Original GitHub release
R/                       Fitting, prediction, and example data functions
src/fused_dot.cpp        C++ helper for kernel selection
examples/borehole.R      Borehole example
examples/michalewicz.R   Michalewicz example
tests/smoke_test.R       Minimal fitting and prediction checks
```
