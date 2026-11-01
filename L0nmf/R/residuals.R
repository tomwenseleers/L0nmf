#' @rdname splogconNMF
#' @param type Residual scale: Pearson or response.
#' @param rows,channels Optional row and channel indices for residual extraction.
#' @param mu Optional fitted means for the selected block.
#' @export
residuals.splogconNMF <- function(object,y=object$y,type=c("pearson","response"),
                                  rows=seq_len(object$n),channels=seq_len(ncol(object$S)),mu=NULL,...) {
  type<-match.arg(type);cols<-channels
  if(is.null(y)) stop("Supply the observed signal matrix as y.")
  Y<-as.matrix(y)
  if(identical(dim(Y),c(object$n,ncol(object$S)))) Y<-Y[rows,cols,drop=FALSE]
  if(!identical(dim(Y),c(length(rows),length(cols)))) stop("y has incompatible dimensions.")
  O<-object$offset;if(length(O)>1L) O<-O[rows,cols,drop=FALSE]
  if(is.null(mu)) mu<-as.matrix(object$H[rows,,drop=FALSE]%*%object$S[,cols,drop=FALSE])+
    object$background[rows,,drop=FALSE]%*%object$alpha[,cols,drop=FALSE]+O
  if(!identical(dim(mu),dim(Y))) stop("mu has incompatible dimensions.")
  if(type=="response") return(Y-mu)
  spec<-object$spec;phi<-spec$dispersion;theta<-spec$theta
  variance<-switch(spec$name,gaussian=matrix(phi,nrow(mu),ncol(mu)),poisson=mu,
    quasipoisson=phi*mu,nb=mu+mu^2/theta,nb1=phi*mu,
    Gamma=phi*mu^2,Gamma_scale=phi*mu)
  expected<-mu
  detection<-object$control$censor
  if(is.null(detection)) detection<-0
  if(spec$name%in%c("Gamma","Gamma_scale")&&detection>0) {
    positive<-which(mu>0)
    a<-if(spec$name=="Gamma") rep(1/phi,length(positive)) else mu[positive]/phi
    b<-if(spec$name=="Gamma") phi*mu[positive] else rep(phi,length(positive))
    expected[]<-variance[]<-0
    expected[positive]<-mu[positive]*stats::pgamma(detection,a+1,scale=b,lower.tail=FALSE)
    second<-a*(a+1)*b^2*stats::pgamma(detection,a+2,scale=b,lower.tail=FALSE)
    variance[positive]<-pmax(0,second-expected[positive]^2)
  }
  W<-object$weights
  if(is.null(W)) W<-1
  if(length(W)>1L) W<-W[rows,cols,drop=FALSE]
  delta<-Y-expected
  residual<-sqrt(W)*delta/sqrt(variance)
  residual[variance==0&delta==0]<-0
  residual[!is.finite(W)|W==0]<-NA_real_
  residual
}

