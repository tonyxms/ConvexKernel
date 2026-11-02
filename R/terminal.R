# One terminal threshold pass, then uniform-start weight optimization if needed.
# Indices in final_pruning refer to the support before terminal pruning.
.optkernel_terminal_reopt <- function(weight, K_sel, inputs, theta, K, state,
                                     threshold, cal_weight) {
  original_weight <- weight
  keep <- weight >= threshold
  # Preserve the frozen algorithm's nonempty-support boundary rule.
  if (!any(keep)) keep[which.max(weight)] <- TRUE
  retained <- which(keep)
  applied <- any(!keep)
  subset_kernel <- function(kernels, indices) {
    if (length(dim(kernels)) == 2L) return(kernels)
    out <- kernels[, , indices, drop = FALSE]
    if (length(indices) == 1L) dim(out) <- dim(out)[1:2]
    out
  }
  if (applied) {
    K_sel <- subset_kernel(K_sel, retained)
    inputs <- inputs[retained, , drop = FALSE]
    theta <- theta[retained]
    fitted <- cal_weight(K_sel)
    weight <- fitted$weight
    K <- fitted$K
    state <- fitted$state
    # Do not threshold again after reoptimization; remove exact zeros only.
    positive <- which(weight > 0)
    if (length(positive) < length(weight)) {
      K_sel <- subset_kernel(K_sel, positive)
      inputs <- inputs[positive, , drop = FALSE]
      theta <- theta[positive]
      weight <- weight[positive]
      retained <- retained[positive]
    }
  }
  removed <- setdiff(seq_along(original_weight), retained)
  list(weight = weight, K_sel = K_sel, kernel_input = inputs, theta = theta,
       K = K, state = state,
       final_pruning = list(enabled = TRUE, applied = applied,
         threshold = threshold, pre_support_count = length(original_weight),
         post_support_count = length(weight), retained_original_indices = retained,
         removed_original_indices = removed,
         removed_weight_mass = sum(original_weight[removed]),
         terminal_final_objective = state$loss))
}
