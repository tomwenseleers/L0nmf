#' Scalable Sparse Nonnegative Matrix Factorization with Automatic Rank Selection and Flexible Noise Models
#'
#' @details
#' For a scans-by-channels signal, each component is the product of an
#' ordered nonnegative profile and a nonnegative vector of channel loadings.
#' The positive part of the profile has a concave logarithm, which accommodates
#' asymmetric peaks and tails. Chromatographic elution profiles and mass
#' spectra are a motivating example.
#'
#' Use \code{\link{estimate_psf}} to build a sparse design from a standard
#' or overlapping signals. \code{\link{splogconNMF}} starts from that design
#' and estimates the number of components, their profiles and their loadings.
#' \code{\link{coupledNMF}} shares loadings across samples and allows their
#' profiles, areas and backgrounds to differ.
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
