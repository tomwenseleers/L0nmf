#' Log-Concave Nonnegative Matrix Factorization with Automatic Rank Selection
#'
#' @details
#' For a scans-by-channels signal, each component is the product of an
#' ordered nonnegative profile and a nonnegative vector of channel loadings.
#' The positive part of the profile has a concave logarithm, which accommodates
#' asymmetric peaks and tails. Chromatographic elution profiles and mass
#' spectra are a motivating example.
#'
#' Use \code{\link{estimate_psf}} to build a sparse design from a standard
#' or overlapping signals. \code{\link{logconNMF}} estimates the number of
#' components, their profiles and their loadings. It fits one matrix or shares
#' spectra across multiple samples with their own profiles, areas and
#' backgrounds.
#'
#' @keywords internal
#' @import L0glm
#' @useDynLib logconNMF, .registration=TRUE
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
