# Source from the repository root; preserve the frozen floating-point settings.
.optkernel_old_cxxflags <- Sys.getenv("PKG_CXXFLAGS", unset = NA_character_)
Sys.setenv(PKG_CXXFLAGS = paste(
  if (is.na(.optkernel_old_cxxflags)) "" else .optkernel_old_cxxflags,
  "-ffp-contract=off -fno-fast-math"))
Rcpp::sourceCpp("src/fused_dot.cpp", env = environment())
if (is.na(.optkernel_old_cxxflags)) Sys.unsetenv("PKG_CXXFLAGS") else
  Sys.setenv(PKG_CXXFLAGS = .optkernel_old_cxxflags)
rm(.optkernel_old_cxxflags)
.optkernel_dot_difference_cpp <- optkernel_dot_difference_cpp
abind <- abind::abind
foreach <- foreach::foreach
`%dopar%` <- foreach::`%dopar%`
registerDoParallel <- doParallel::registerDoParallel
stopImplicitCluster <- doParallel::stopImplicitCluster
source("R/optkernel.R", local = TRUE)
source("R/terminal.R", local = TRUE)
source("R/loocv.R", local = TRUE)
source("R/predict.R", local = TRUE)
