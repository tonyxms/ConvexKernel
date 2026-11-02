# Run from the repository root: Rscript --vanilla examples/michalewicz.R
source("R/load_all.R")
source("R/benchmarks.R")

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(201)
data <- prepare_michalewicz(n = 60, d = 6, m = 100)

# Six eta candidates, the full kernel dictionary, and the paper's search settings.
model <- optkernel_fun_loocv(data$x, data$y, mu_seeds = 2101:2106)
prediction <- optkernel_pred(data$x, model, data$newx)
metrics <- example_metrics(prediction, data$truth)
selected_inputs <- which(colSums(model$kernel_input[model$weight > 0, , drop = FALSE]) > 0)

print(attr(model, "loocv"))
print(metrics)
cat("Selected input columns:", selected_inputs, "\n")
cat("Analytically active input columns:", 1:2, "\n")
