#' Estimate a shared noise model for ordered multichannel signals
#'
#' @param y Observations-by-channels matrix.
#' @param family Identity-link family accepted by [splogconNMF()].
#' @param dispersion,theta Fixed noise parameters, or NULL to estimate them.
#' @param weights Scalar, row vector or matrix of observation precisions.
#' @param correlation Correlation convention used by [splogconNMF()].
#' @param censor Gamma detection limit.
#' @param ncores Number of physical worker cores; NULL detects them.
#' @return A list containing the fitted noise specification, pilot means and
#'   equivalent Poisson working precision for PSF diagnostics.
#' @details Pooled marginal variances are estimated from temporal differences.
#'   A shared information factor accounts for estimated lag-one correlation.
#'   The same calibration can be reused during PSF estimation and factorization.
#' @export
estimate_signal_noise <- function(y,family=quasipoisson("identity"),dispersion=NULL,
                                  theta=NULL,weights=NULL,correlation="none",
                                  censor=.5,ncores=NULL) {
  Y<-as.matrix(y)
  if (!is.numeric(Y)||nrow(Y)<9L||!ncol(Y)||any(!is.finite(Y))) stop("Use a finite numeric signal matrix with at least nine rows.")
  ncores<-nmf_threads(ncores)
  spec<-slcn_family(family,dispersion,theta)
  W<-slcn_weights(weights,nrow(Y),ncol(Y))
  if (any(!is.finite(W)|W<0)) stop("Use finite nonnegative weights.")
  noise<-slcn_noise(Y,spec,W,ncores,slcn_correlation(correlation),censor)
  spec<-noise$spec;factor<-spec$information_factor
  if(spec$name %in% c("poisson","quasipoisson","nb1","Gamma_scale")) {
    working_weight<-1/(spec$dispersion*factor)
  } else {
    level<-if(is.null(noise$pilot)) slcn_pilot(Y) else noise$pilot
    total<-mean(rowSums(level))
    variance<-switch(spec$name,gaussian=ncol(Y)*spec$dispersion,
      nb=mean(rowSums(level+level^2/spec$theta)),Gamma=spec$dispersion*mean(rowSums(level^2)))
    working_weight<-total/max(variance,1e-15)/factor
  }
  noise$working_weight<-working_weight
  noise
}
