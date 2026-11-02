#include <Rcpp.h>

// Subtract before multiplying, as in the reference R expression.
// Read-only inputs; no temporary vector, parallelism or RNG calls.
// [[Rcpp::export(rng = false)]]
double optkernel_dot_difference_cpp(const Rcpp::NumericVector& weight,
                                   const Rcpp::NumericVector& packed,
                                   const Rcpp::NumericVector& current) {
  const R_xlen_t n = weight.size();
  const double* w = weight.begin();
  const double* p = packed.begin();
  const double* c = current.begin();
  double total = 0.0;
  for (R_xlen_t j = 0; j < n; ++j) {
    const double delta = p[j] - c[j];
    total += w[j] * delta;
  }
  return total;
}
