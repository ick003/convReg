#' convReg : Package for convolutive mixture models.
#'
#' convReg allows you to infer convolutive mixture models.
#' @useDynLib convReg
#' @importFrom Rcpp sourceCpp
#' @importFrom grDevices hcl rgb
#' @importFrom graphics abline arrows axis filled.contour grid hist layout legend par plot plot.new points polygon text
#' @importFrom utils combn
#' @importFrom COMPoissonReg glm.cmp rcmp dcmp
#' @import RcppArmadillo compiler data.table doMC doParallel foreach inline lattice nnet numDeriv parallel pscl stats MASS
"_PACKAGE"
