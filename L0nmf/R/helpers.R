nmf_threads <- function(nthreads) {
  if (is.null(nthreads)) nthreads <- parallel::detectCores(logical=FALSE)
  if (length(nthreads)==1L && is.na(nthreads)) nthreads <- 1L
  if (!is.numeric(nthreads)||length(nthreads)!=1L||!is.finite(nthreads)||
      nthreads<1||nthreads!=floor(nthreads)||nthreads>.Machine$integer.max)
    stop("nthreads must be a positive integer.",call.=FALSE)
  as.integer(nthreads)
}
nmf_matrix <- function(value,n,m,default,name) {
  if (is.null(value)) return(matrix(default,n,m))
  if (length(value)==1L) return(matrix(as.numeric(value),n,m))
  if (length(value)==n && (is.null(dim(value))||NCOL(value)==1L)) return(matrix(as.numeric(value),n,m))
  if (length(value)!=n*m) stop(name," has incompatible dimensions.",call.=FALSE)
  matrix(as.numeric(value),n,m)
}
nmf_ic_name <- function(IC) {
  if (!is.character(IC)||length(IC)!=1L||is.na(IC)||!tolower(IC)%in%c("aic","bic","gic","ebic","mbic"))
    stop("Choose aic, bic, gic, ebic or mbic.",call.=FALSE)
  tolower(IC)
}
nmf_model_criterion <- function(object,criterion,mu_floor=0) {
  ll <- stats::logLik(object,mu_floor=mu_floor)
  -2*as.numeric(ll)+preset.lambda(criterion,p=max(1,object$p),n=attr(ll,"nobs"))*attr(ll,"df")
}
