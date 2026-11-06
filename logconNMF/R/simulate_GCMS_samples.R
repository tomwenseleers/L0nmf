cpl_emg_density <- function(t, center, sigma, tail) {
  if(tail<=1e-10) return(dnorm(t,center,sigma))
  z <- (t-center)/sigma-sigma/tail
  exp(-log(tail)+(center-t)/tail+sigma^2/(2*tail^2)+pnorm(z,log.p=TRUE))
}

#' Simulate related chromatograms with shared compound spectra
#'
#' Generate exponentially modified Gaussian and optional triangular acid peaks,
#' smooth backgrounds, irregular scan grids, sample drift and concentration-
#' dependent peak widths and tails. Missing compounds can be simulated explicitly.
#' @param samples Number of samples.
#' @param n Number of scans, recycled across samples.
#' @param k Number of distinct compounds.
#' @param S Optional spectral library with one compound per row.
#' @param preset Peak-shape scenario: identical shapes, time-dependent shapes,
#'   time and concentration-dependent shapes, or the latter with acid peaks.
#' @param presence Probability that each compound occurs in each sample.
#' @param absent Optional logical compounds-by-samples absence matrix.
#' @param irregular Relative scan-width variation, from zero to 0.5.
#' @param drift Drift per sample, in seconds.
#' @param area Typical integrated compound intensity.
#' @param seed Random-number seed.
#' @param width_time,width_area,tail_time,tail_area,rt_area Effects of time and
#'   concentration on peak width, tail and retention time.
#' @param acid_fraction Fraction of broad triangular acid profiles.
#' @param family Signal distribution name supported by [coupledNMF()].
#' @param dispersion,theta Noise parameters.
#' @return Simulation truth and a list of samples with observations, original
#'   bin edges and times (minutes), integrated-count means, true spectra,
#'   profiles, areas, presence indicators and acquisition dates.
#'   A simulated observation is integrated over its bin; `raw_exposure` is the
#'   original bin duration in seconds. For [coupledNMF()], supply it as exposure.
#' @examples
#' sim <- simulate_GCMS_samples(samples=3L,n=c(180L,190L,200L),k=8L,
#'   presence=0.5,seed=3)
#' sim$present
#' @seealso [coupledNMF()], [simulate_GCMS()]
#' @md
#' @export
simulate_GCMS_samples <- function(samples=3L,n=600L,k=26L,S=NULL,
    preset=c("time_area_acid","same","time","time_area"),presence=1,
    absent=NULL,irregular=.12,drift=2,area=2e5,seed=811,
    width_time=.5,width_area=.45,tail_time=1,tail_area=1.5,rt_area=.25,acid_fraction=.08,
    family=c("poisson","quasipoisson","nb","gaussian","Gamma","nb1","Gamma_scale"),dispersion=1,theta=20) {
  preset<-match.arg(preset);family<-match.arg(family)
  if((family=="nb1"&&dispersion<1)||samples<1||k<2||any(rep_len(n,samples)<50)||presence<0||presence>1||
    irregular<0||irregular>.5||dispersion<=0||theta<=0) stop("Invalid simulation parameters.")
  set.seed(seed);ns<-as.integer(rep_len(n,samples));duration<-max(ns)
  if(is.null(S)) {e<-new.env();data("simulated_spectra",package="L0glm",envir=e);S<-as.matrix(e$simulated_spectra)}
  ids<-sample(nrow(S),k);P<-S[ids,,drop=FALSE];P<-P/rowSums(P);m<-ncol(P)
  # Include a close pair with distinct spectra; the matching algorithm never sees truth IDs.
  centers<-seq(30,duration-35,length.out=k)
  if(k>=10L) centers[10]<-centers[9]+5
  basearea<-area*exp(rnorm(k,0,.5));A<-outer(basearea,rep(1,samples))*exp(matrix(rnorm(k*samples,0,.65),k,samples))
  present<-matrix(runif(k*samples)<presence,k,samples)
  if(!is.null(absent)) {stopifnot(all(dim(absent)==c(k,samples)));present<-!absent}
  A[!present]<-0
  acid<-if(preset=="time_area_acid"&&acid_fraction>0)utils::tail(seq_len(k),min(k,max(1L,round(k*acid_fraction)))) else integer()
  backgrounds<-matrix(runif(4L*m,.5,2),4L,m)
  result<-vector("list",samples)
  for(s in seq_len(samples)) {
    # The native/reference time bins differ between injections.
    e<-runif(ns[s],1-irregular,1+irregular);e<-e/sum(e)*duration
    edges<-c(0,cumsum(e));u<-(utils::head(edges,-1L)+utils::tail(edges,-1L))/2
    shift<-function(v) drift*(s-1L)*(1+.3*sin(2*pi*v/duration))
    rawedges<-edges+shift(edges);time<-(utils::head(rawedges,-1)+utils::tail(rawedges,-1))/2
    H<-matrix(0,ns[s],k);sigmas<-tails<-numeric(k)
    for(j in seq_len(k)) {
      large<-if(A[j,s]>0)max(0,log(A[j,s]/area)) else 0
      sigma<-1.8*(if(preset=="same")1 else 1+width_time*centers[j]/duration)
      tau<-.8*(if(preset=="same")1 else 1+tail_time*centers[j]/duration)
      if(preset%in%c("time_area","time_area_acid")){sigma<-sigma*(1+width_area*large);tau<-tau*(1+tail_area*large)}
      center<-centers[j]+if(preset%in%c("time_area","time_area_acid")) rt_area*large else 0
      # Composite midpoint integration gives known truth on each original bin.
      for(q in seq_len(8L)) {
        t<-utils::head(edges,-1L)+(q-.5)/8*e
        h<-if(j%in%acid) ifelse(t>=center-22 & t<=center+6,
          pmax(0,(t-center+22)/22)*exp(-pmax(t-center,0)/.8),0) else cpl_emg_density(t,center,sigma,tau)
        H[,j]<-H[,j]+h/8
      }
      H[,j]<-H[,j]/sum(H[,j]*e);sigmas[j]<-sigma;tails[j]<-tau
    }
    B<-splines::bs(u,df=4L,intercept=TRUE);background<-B%*%backgrounds
    rates<-H%*%(P*A[,s])+background;mu<-rates*e
    Y<-switch(family,poisson=matrix(rpois(length(mu),mu),ns[s],m),
      quasipoisson=matrix(rpois(length(mu),mu/dispersion)*dispersion,ns[s],m),
      nb=matrix(stats::rnbinom(length(mu),mu=mu,size=theta),ns[s],m),
      gaussian=mu+matrix(rnorm(length(mu),sd=sqrt(dispersion)),ns[s],m),
      Gamma=matrix(rgamma(length(mu),shape=1/dispersion,scale=dispersion*mu),ns[s],m),
      nb1=if(dispersion==1)matrix(rpois(length(mu),mu),ns[s],m)else
        matrix(stats::rnbinom(length(mu),mu=mu,size=mu/(dispersion-1)),ns[s],m),
      Gamma_scale=matrix(rgamma(length(mu),shape=mu/dispersion,scale=dispersion),ns[s],m))
    if(family%in%c("Gamma","Gamma_scale"))Y[Y<.5]<-0
    storage.mode(Y)<-"double";colnames(Y)<-colnames(P);rownames(Y)<-format(time/60,digits=12)
    result[[s]]<-list(Y=Y,time=time/60,edges=rawedges/60,reference_time=u/60,
      reference_edges=edges/60,exposure=e,raw_exposure=diff(rawedges),
      mean=mu,rates=rates,H_true=H,S_true=P*A[,s],P_true=P,present=present[,s],
      area_true=A[,s],centers=centers/60,sigma=sigmas,tail=tails,acid=acid,
      background=background*e,family=family,dispersion=dispersion,theta=theta,
      sample_id=sprintf("sample_%03d",s),acquired_at=as.POSIXct("2026-11-01 08:00:00",tz="UTC")+3600*(s-1L))

  }
  structure(list(samples=result,P=P,compound_ids=ids,present=present,areas=A,
    preset=preset,seed=seed,settings=list(presence=presence,irregular=irregular,drift=drift,
      family=family,dispersion=dispersion,theta=theta)),class="GCMS_simulated_samples")
}
