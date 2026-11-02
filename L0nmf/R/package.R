#' Scalable Sparse Nonnegative Matrix Factorization with Automatic Rank Selection and Flexible Noise Models
#'
#' @keywords internal
#' @import L0glm
#' @useDynLib L0nmf, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @import Matrix
#' @importFrom stats AIC BIC Gamma approx coef cor dgamma dnorm fitted gaussian integrate lm logLik mad median nobs optim pchisq plogis pnorm poisson predict qnorm quantile quasipoisson rbinom rgamma rlnorm rnorm rpois runif sd setNames var
#' @importFrom graphics abline axis box layout legend lines matplot mtext par plot points polygon segments text title
#' @import grDevices methods
#' @importFrom utils data flush.console
"_PACKAGE"

#' @importFrom L0glm GIC
#' @export
L0glm::GIC
#' @importFrom L0glm EBIC
#' @export
L0glm::EBIC
#' @importFrom L0glm MBIC
#' @export
L0glm::MBIC

utils::globalVariables("simulated_spectra")
