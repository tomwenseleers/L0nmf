#' Coupled sparse log-concave factorization of multiple samples
#'
#' Align related multichannel samples and estimate common spectra with
#' sample-specific nonnegative areas, log-concave profiles and backgrounds.
#' Observations remain on their original time grids. Samples are read from a
#' disk-backed store one at a time during fitting.
#'
#' @param samples List of sample descriptions. Each contains `y` (a numeric
#'   scans-by-channels matrix, or an RDS file containing one), increasing `time`,
#'   and optionally `X`, `background`, `fit`, `weights`, `offset`, `exposure`,
#'   `sample_id`, `group`, `small_width` and `edges`. Large `weights` and `offset`
#'   matrices can also be supplied as RDS paths. `fit` is a [splogconNMF()]
#'   initializer or its RDS path. Without a fit, `X` supplies candidate profiles;
#'   otherwise [estimate_psf()] estimates them. Channel names and ordering must
#'   agree across samples. Times and `small_width` use the same units.
#'   Shared scan-wise weights and offsets remain vectors throughout coupled
#'   fitting; each mass channel reuses the same values.
#' @param family Identity-link signal family supported by [splogconNMF()]:
#'   Gaussian, Poisson, quasipoisson, NB2, NB1, censored Gamma or common-scale Gamma.
#' @param lambda Information-criterion preset or nonnegative penalty per
#'   selection degree of freedom. The selection count is
#'   `(channels - 1) * shared_compounds + detected_occurrences`, plus estimated
#'   noise parameters. Background and profile complexity are unpenalized.
#' @param path New directory for the sample store and fitted factors. Keep this
#'   directory with the saved R object; its paths are used by `fitted()` and plots.
#' @param precluster Maximum outcome groups for independent initialization.
#'   Coupled refinement and selection always use all original channels.
#' @param dispersion,theta Shared-across-channels noise parameters, as in
#'   [splogconNMF()]. They may be supplied per sample in its description.
#'   Resolved dispersions are frozen; a negative NB2 theta enables bounded
#'   likelihood-improving shape updates in the coupled iterations.
#' @param method Spectral distance for anchor matching and grouping.
#' @param align Estimate monotone sample-level time warps.
#' @param ncores Number of native OpenMP threads; defaults to physical cores.
#' @param control Named controls: `cycles=10L`, `initial_cycles=5L`, `rounds=3L`,
#'   `q=64L`, `rt_tolerance=NULL`, `alignment_window=NULL`, `ridge=TRUE`,
#'   `censor=0.5`, `correlation="none"`, and `initial=list()` for controls passed
#'   to independent [splogconNMF()] fits. Default matching windows are derived
#'   from the supplied or estimated narrow-peak widths.
#' @param verbose Print initialization and coupled objective progress.
#' @param object,x A fitted coupled model.
#' @param sample One sample index or name; `NULL` in `coef()` returns consensus
#'   spectra, and in `plot()` selects all samples.
#' @param normalize Return unit-sum spectra or amplitude-scaled spectra.
#' @param type Plot sample areas, common spectra or aligned compound profiles.
#' @param compounds Compound indices to display.
#' @param ... Additional method arguments.
#' @param k Information-criterion multiplier for `AIC()`.
#'
#' @details
#' Let `u = w_s(t)` be a monotone alignment map and `J` its binwise Jacobian.
#' With supplied exposure `e`, the original-grid mean is
#' \deqn{\mu_{sij}=e_{si}\{J_{si}\sum_g a_{sg}h_{sg}(u_{si})p_{gj}
#' +(B_s A_s)_{ij}\}+o_{sij}.}
#' Spectra satisfy \eqn{\sum_j p_{gj}=1}; profiles satisfy
#' \eqn{\sum_i\Delta u_{si}h_{sg}(u_{si})=1} and have concave logarithms in
#' actual aligned time. Thus areas integrate the fitted rate over original time.
#' `exposure=1` models the supplied intensities. For integrated counts, known
#' integration durations can instead be supplied as exposure. Scan spacing,
#' per-mass dwell time and vendor abundance units are distinct quantities.
#' No exposure is inferred from the number of nonzero mass traces.
#'
#' Reciprocal spectral matches supply leading-edge anchors: the 1%-height start
#' plus half the local narrow-peak width. A robust low-dimensional drift spline
#' includes log area ratio as a nuisance covariate; only drift enters the warp.
#' Positive-slope projection gives a monotone map. Grouping pools compatible
#' spectra algebraically, with at most one occurrence per sample and an
#' area-ordered profile-overlap check. It performs no model refits. Subsequent
#' local profile, pooled spectral and background updates are backtracked against
#' the full original-channel likelihood. A grouped proposal is retained only
#' when the complete penalized objective improves.
#'
#' The grouping criterion counts shared spectral parameters and detected
#' sample-specific amplitudes; hard shape constraints control profile geometry.
#' It is a rank-selection criterion with unpenalized profile and baseline
#' complexity. Quasipoisson uses the scaled likelihood kernel described in
#' [splogconNMF()]. Noise scales are held fixed when comparing groupings.
#' Undetected sample-compound combinations have `NA` area and are not filled in.
#'
#' @return A `coupledNMF` object with consensus spectra `P`, integrated `areas`,
#'   occurrence table `compounds`, monotone `alignment`, disk-backed factor
#'   `paths` and sample `store`, objective history, likelihood and timings.
#'   `fitted(object, sample=...)` reconstructs one sample on its original grid.
#' @seealso [splogconNMF()], [estimate_psf()], [simulate_GCMS_samples()]
#' @examples
#' set.seed(19)
#' n <- 100L; m <- 6L
#' X <- convmatrix(dnorm(-6:6, sd=1.5), n)
#' H <- X[, c(25, 65)]; H <- H %*% Matrix::Diagonal(x=1/colSums(H))
#' S <- matrix(rexp(2*m, 1/2000), 2, m)
#' samples <- lapply(1:3, function(s) list(
#'   y=matrix(rpois(n*m, as.matrix(H %*% (S*s))+1), n, m),
#'   time=seq_len(n), X=X, background=matrix(1,n,1), small_width=9))
#' fit <- coupledNMF(samples, ncores=1L,
#'   control=list(cycles=2L, initial_cycles=1L, rounds=1L,
#'     initial=list(profile_cycles=0L, original_cleanup=FALSE)))
#' fit$areas
#' GIC(fit)
#' @md
#' @export
coupledNMF <- function(samples, family=poisson("identity"), lambda="gic",
                       path=tempfile("coupledNMF_"), precluster=32L,
                       dispersion=NULL, theta=NULL, method=c("cosine","hellinger"),
                       align=TRUE, ncores=NULL, control=list(), verbose=FALSE) {
  begun<-proc.time()[["elapsed"]];method<-match.arg(method)
  if(!is.list(samples)||length(samples)<2L)stop("Supply at least two sample descriptions.")
  if(!is.logical(align)||length(align)!=1L||is.na(align))stop("align must be TRUE or FALSE.")
  if(is.null(ncores))ncores<-parallel::detectCores(logical=FALSE)
  if(length(ncores)!=1L||!is.finite(ncores)||ncores<1||ncores!=floor(ncores))stop("Invalid ncores.")
  defaults<-list(cycles=10L,initial_cycles=5L,rounds=3L,q=64L,rt_tolerance=NULL,
    alignment_window=NULL,ridge=TRUE,censor=.5,correlation="none",initial=list())
  if(!is.list(control)||length(control)&&(is.null(names(control))||anyDuplicated(names(control))||
      length(setdiff(names(control),names(defaults)))))stop("Use named coupled controls.")
  ctl<-utils::modifyList(defaults,control)
  if(is.list(family)&&!is.null(family$scalar)&&!"censor"%in%names(control))ctl$censor<-family$scalar$censor
  for(nm in c("cycles","initial_cycles","rounds","q")) {
    v<-ctl[[nm]];if(length(v)!=1L||!is.finite(v)||v<(if(nm=="q")3 else 0)||v!=floor(v))stop("Invalid ",nm)
  }
  if(!is.logical(ctl$ridge)||length(ctl$ridge)!=1||is.na(ctl$ridge))stop("Invalid ridge control.")
  if(!is.finite(ctl$censor)||length(ctl$censor)!=1L||ctl$censor<0)stop("Invalid censor.")
  if(dir.exists(path))stop("Use a new, empty store path.")
  dir.create(path,recursive=TRUE);path<-normalizePath(path,winslash="/",mustWork=TRUE)
  initial<-system.time(store<-cpl_prepare(samples,path,family,dispersion,theta,precluster,ncores,ctl,verbose))[["elapsed"]]
  state<-cpl_initialize(store,method,align,ncores,lambda,ctl)
  state<-cpl_optimize_coupled(state,cycles=ctl$initial_cycles,q=ctl$q)
  reference<-cpl_coupled_score(state);candidate<-state
  for(round in seq_len(ctl$rounds)) {
    before<-length(unique(candidate$map))
    candidate<-cpl_batch_spectral_groups(candidate,method,ctl$ridge)
    if(round>1L&&length(unique(candidate$map))==before)break
    candidate<-cpl_optimize_coupled(candidate,cycles=ctl$cycles,q=ctl$q)
    if(verbose)cat("Coupling round",round,":",length(unique(candidate$map)),"compounds; objective",cpl_coupled_score(candidate),"\n")
  }
  candidate$batch_delta<-cpl_coupled_score(candidate)-reference
  accepted<-is.finite(candidate$batch_delta)&&candidate$batch_delta<=0
  if(accepted)state<-candidate
  state$batch_delta<-candidate$batch_delta;state$batch_accepted<-accepted
  state<-cpl_compact_coupled(state)
  state$loglik<- -cpl_coupled_score(state,penalty=FALSE)/2
  state$ic_df<-(ncol(state$P)-1L)*nrow(state$P)+state$L+state$nuisance_df
  state$objective<- -2*state$loglik+state$lambda*state$ic_df
  state$seconds<-c(initial=initial,coupled=proc.time()[["elapsed"]]-begun-initial,total=proc.time()[["elapsed"]]-begun)
  # Only immutable checkpoints of this private store are eligible for cleanup.
  for(d in store$manifest$directory) {
    obsolete<-setdiff(list.files(d,pattern="^fit_.*\\.rds$",full.names=TRUE),state$paths)
    if(length(obsolete))unlink(obsolete)
    meta<-readRDS(file.path(d,"metadata.rds"));meta$initial<-NULL
    saveRDS(meta,file.path(d,"metadata.rds"),compress=FALSE)
  }
  saveRDS(state,file.path(path,"model.rds"))
  state
}

cpl_edges <- function(time) {
  d<-diff(time);c(time[1L]-d[1L]/2,(utils::head(time,-1L)+utils::tail(time,-1L))/2,utils::tail(time,1L)+utils::tail(d,1L)/2)
}
cpl_prepare <- function(samples,path,family,dispersion,theta,precluster,ncores,ctl,verbose) {
  entries<-vector("list",length(samples));axis<-NULL
  for(s in seq_along(samples)) {
    x<-samples[[s]];if(!is.list(x)||is.null(x$y)||is.null(x$time))stop("Each sample needs y and time.")
    Y<-if(is.character(x$y))readRDS(x$y)else x$y;Y<-as.matrix(Y);storage.mode(Y)<-"double"
    n<-nrow(Y);m<-ncol(Y);time<-as.numeric(x$time)
    if(n<9L||!m||length(time)!=n||any(!is.finite(time))||any(diff(time)<=0)||any(!is.finite(Y)))stop("Invalid sample data/time.")
    mass<-colnames(Y);if(is.null(mass))mass<-as.character(seq_len(m))
    if(is.null(axis))axis<-mass
    if(!identical(mass,axis))stop("Sample channel names and ordering must agree.")
    edges<-if(is.null(x$edges))cpl_edges(time)else x$edges
    if(length(edges)!=n+1L||any(!is.finite(edges))||any(diff(edges)<=0)||any(time<=utils::head(edges,-1L)|time>=utils::tail(edges,-1L)))stop("Invalid bin edges.")
    E<-if(is.null(x$exposure))rep(1,n)else rep_len(x$exposure,n)
    if(!length(x$exposure)%in%c(0L,1L,n)||any(!is.finite(E)|E<=0))stop("Exposure must be positive and scalar or scan-specific.")
    value<-function(x)if(is.character(x))readRDS(x)else x
    W<-cpl_observation(value(x$weights),n,m,1);O<-cpl_observation(value(x$offset),n,m,0)
    if(any(!is.finite(W)|W<0)||!any(W>0)||any(!is.finite(O)|O<0))stop("Invalid weights/offset.")
    B<-if(is.null(x$background))matrix(1,n,1L)else as.matrix(x$background)
    if(nrow(B)!=n||any(!is.finite(B)|B<0))stop("Invalid background basis.")
    f<-if(is.character(x$fit))readRDS(x$fit)else x$fit
    if(!is.null(f$model))f<-f$model
    if(!is.null(f)) {
      if(is.null(x$weights))W<-cpl_observation(f$weights,n,m,1)
      if(is.null(x$offset))O<-cpl_observation(f$offset,n,m,0)
      if(any(!is.finite(W)|W<0)||!any(W>0)||any(!is.finite(O)|O<0))stop("Invalid initial-fit weights/offset.")
    }
    if(is.null(f)) {
      X<-x$X
      if(is.null(X))X<-estimate_psf(pmax(Y,0),time=time,family=family,nwindows=1L,precluster=1L,ncores=ncores)$X
      f<-splogconNMF(X,Y,family=family,background=B*E,lambda="gic",precluster=precluster,
        dispersion=if(is.null(x$dispersion))dispersion else x$dispersion,
        theta=if(is.null(x$theta))theta else x$theta,weights=W,offset=O,
        correlation=ctl$correlation,ncores=ncores,control=ctl$initial,verbose=verbose)
    }
    if(nrow(f$H)!=n||ncol(f$S)!=m||ncol(f$H)!=nrow(f$S)||any(!is.finite(f$S)|f$S<0)||
        any(!is.finite(f$H@x)|f$H@x<0))stop("Invalid initial factors.")
    if(is.null(f$spec))stop("Initial fits must contain resolved noise settings in spec.")
    if(f$spec$kind!=0L&&any(Y<0))stop("Nonnegative observations are required by this family.")
    if(f$spec$kind%in%c(1L,3L,5L)&&any(abs(Y-round(Y))>1e-7))stop("Count models require integer observations.")
    if(f$spec$kind%in%c(4L,6L)&&any(Y>0&Y<ctl$censor))stop("Represent values below the Gamma detection limit by zero.")
    id<-if(is.null(x$sample_id))paste0("sample",s)else as.character(x$sample_id)
    if(length(id)!=1||is.na(id)||!nzchar(id))stop("Invalid sample_id.")
    d<-file.path(path,sprintf("sample_%06d",s));dir.create(d)
    cpl_write_matrix_file(Y,file.path(d,"observed.bin"))
    small<-x$small_width
    if(is.null(small))small<-2*sqrt(2*log(100))*1.6*median(diff(time))
    if(length(small)!=1L||!is.finite(small)||small<=0)stop("small_width must be positive.")
    # Compact initial factors; observed means and reduced-space fits are released.
    f<-f[c("H","S","background","alpha","spec")]
    if(is.null(f$background))f$background<-B*E
    initial<-cpl_checkpoint_sample(f,d,"initial")
    wfile<-if(length(W)>1L)file.path(d,"weights.rds")else NULL
    ofile<-if(length(O)>1L)file.path(d,"offset.rds")else NULL
    if(!is.null(wfile))saveRDS(W,wfile,compress=FALSE)
    if(!is.null(ofile))saveRDS(O,ofile,compress=FALSE)
    meta<-list(n=n,m=m,time=time,edges=edges,raw_exposure=E,mass=mass,sample_id=id,
      group=if(is.null(x$group))"samples"else as.character(x$group),small_width=small,
      initial=initial,weights=if(is.null(wfile))W else wfile,
      offset=if(is.null(ofile))O else ofile,censor=ctl$censor)
    saveRDS(meta,file.path(d,"metadata.rds"),compress=FALSE)
    entries[[s]]<-data.frame(sample_id=id,directory=d,n=n,m=m,group=meta$group)
    if(verbose)cat("Initialized",id,":",ncol(f$H),"occurrences\n")
  }
  manifest<-do.call(rbind,entries)
  if(anyDuplicated(manifest$sample_id))stop("Use unique sample_id values.")
  saveRDS(manifest,file.path(path,"manifest.rds"),compress=FALSE)
  list(path=path,manifest=manifest)
}

cpl_initialize <- function(store,method,align,ncores,lambda,ctl) {
  meta<-lapply(store$manifest$directory,function(d)readRDS(file.path(d,"metadata.rds")))
  # Retain spectra and anchor summaries globally; load profile matrices singly.
  fits<-lapply(meta,function(m) {
    f<-readRDS(m$initial)
    list(S=f$S,spec=f$spec,geometry=cpl_profile_anchors(f$H,f$S,m$time,m$small_width))
  })
  if(is.null(ctl$rt_tolerance))ctl$rt_tolerance<-max(vapply(meta,`[[`,0,"small_width"))
  if(is.null(ctl$alignment_window))ctl$alignment_window<-5*ctl$rt_tolerance
  for(key in c("rt_tolerance","alignment_window"))if(length(ctl[[key]])!=1L||!is.finite(ctl[[key]])||ctl[[key]]<=0)stop("Invalid ",key)
  alignment<-if(align)cpl_estimate_alignment(fits,meta,method,lapply(meta,`[[`,"small_width"),ctl$alignment_window,
    geometry=lapply(fits,`[[`,"geometry"))else
    list(maps=lapply(meta,function(x)list(evaluate=identity,method="identity")),reference=1L,method="none")
  N<-sum(store$manifest$n);m<-store$manifest$m[1L]
  factor<-max(vapply(fits,function(f)if(is.null(f$spec$information_factor))1 else f$spec$information_factor,0))
  scalar<-if(is.character(lambda))preset.lambda(lambda,p=N,n=N)else lambda
  if(length(scalar)!=1L||!is.finite(scalar)||scalar<0)stop("Invalid lambda.")
  total<-sum(vapply(fits,function(f)nrow(f$S),0L))
  state<-list(store=store,paths=character(length(meta)),P=matrix(0,total,m),map=seq_len(total),
    scale=rep(1,total),serial=0L,lambda=scalar*factor,L=0L,nuisance_df=0L,alignment=alignment,
    history=data.frame(),control=ctl,mass=meta[[1]]$mass,n=N,p=N,quasi=FALSE,information_factor=factor)
  for(s in seq_along(meta)) {
    z<-meta[[s]];a<-readRDS(z$initial);fits[s]<-list(NULL);warp<-alignment$maps[[s]]$evaluate
    quadrature<-diff(warp(z$edges));J<-quadrature/diff(z$edges);E<-z$raw_exposure*J
    time<-warp(z$time)
    if(any(!is.finite(E)|E<=0)||any(diff(time)<=0))stop("Alignment must have positive bin widths.")
    H<-a$H;K<-ncol(H);area<-rowSums(a$S);P<-a$S/area
    for(k in seq_len(K)) {
      at<-which(H[,k]>max(H[,k])*1e-10)
      h<-as.numeric(H[at,k])/E[at]
      if(length(at)>2L)h<-cpl_project_irregular_lc(h,time[at])
      scale<-sum(quadrature[at]*h);H[,k]<-0;H[at,k]<-h/scale;area[k]<-area[k]*scale
    }
    gid<-state$L+seq_len(K);spec<-a$spec
    state$nuisance_df<-state$nuisance_df+as.integer(spec$dispersion_estimated)+as.integer(spec$theta_estimated)
    state$quasi<-state$quasi||spec$kind==2L
    f<-list(H=drop0(as(H,"dgCMatrix")),area=area,gid=gid,aligned_time=time,
      exposure=E,quadrature=quadrature,raw_exposure=z$raw_exposure,
      background=a$background/z$raw_exposure,alpha=a$alpha,offset=z$offset,
      kind=spec$kind,phi=spec$dispersion,theta=spec$theta,spec=spec,censor=z$censor,
      weights=if(is.character(z$weights))readRDS(z$weights)else z$weights,ncores=as.integer(ncores),sample_id=z$sample_id,
      small_width=z$small_width,time=z$time,edges=z$edges)
    if(is.character(f$offset))f$offset<-readRDS(f$offset)
    state$P[gid,]<-P
    state$L<-state$L+K;state<-cpl_coupled_checkpoint(state,s,f,"aligned")
  }
  colnames(state$P)<-state$mass
  state
}

cpl_sample_index <- function(object,sample) {
  s<-if(is.character(sample))match(sample,object$store$manifest$sample_id)else sample
  if(length(s)!=1L||is.na(s)||s<1||s>length(object$paths)||s!=floor(s))stop("Select one sample by index or name.")
  s
}
#' @rdname coupledNMF
#' @export
fitted.coupledNMF <- function(object,sample=1L,...) {
  s<-cpl_sample_index(object,sample);f<-readRDS(object$paths[s])
  out<-cpl_coupled_mean(f,object);colnames(out)<-object$mass;out
}
#' @rdname coupledNMF
#' @export
coef.coupledNMF <- function(object,sample=NULL,normalize=c("none","l1"),...) {
  normalize<-match.arg(normalize)
  if(is.null(sample))return(if(normalize=="l1")object$P else object$consensus_counts)
  f<-readRDS(object$paths[cpl_sample_index(object,sample)]);ids<-object$map[f$gid]
  out<-object$P[ids,,drop=FALSE]
  if(normalize=="none")out<-out*(f$area*object$scale[ids])
  out
}
#' @rdname coupledNMF
#' @export
logLik.coupledNMF <- function(object,...)structure(object$loglik,class="logLik",df=object$ic_df,nobs=object$n,quasi=object$quasi)
#' @rdname coupledNMF
#' @export
nobs.coupledNMF <- function(object,...)object$n
#' @rdname coupledNMF
#' @export
AIC.coupledNMF <- function(object,...,k=2)-2*object$loglik+k*object$information_factor*object$ic_df
#' @rdname coupledNMF
#' @export
BIC.coupledNMF <- function(object,...)AIC(object,k=log(object$n))
#' @rdname coupledNMF
#' @export
GIC.coupledNMF <- function(object,...)AIC(object,k=preset.lambda("gic",p=object$p,n=object$n))
#' @rdname coupledNMF
#' @export
EBIC.coupledNMF <- function(object,...)AIC(object,k=preset.lambda("ebic",p=object$p,n=object$n))
#' @rdname coupledNMF
#' @export
MBIC.coupledNMF <- function(object,...)AIC(object,k=preset.lambda("mbic",p=object$p,n=object$n))
#' @rdname coupledNMF
#' @method print coupledNMF
#' @export
print.coupledNMF <- function(x,...) {
  cat("Coupled log-concave NMF:",length(x$paths),"samples;",nrow(x$P),"shared compounds;",x$L,"occurrences\n")
  cat("Objective:",format(x$objective,digits=9),"; time:",round(x$seconds[["total"]],2),"s\n")
  invisible(x)
}
#' @rdname coupledNMF
#' @export
plot.coupledNMF <- function(x,type=c("areas","spectra","profiles"),sample=NULL,compounds=seq_len(min(6L,nrow(x$P))),...) {
  type<-match.arg(type)
  if(!length(compounds)||any(!compounds%in%seq_len(nrow(x$P))))stop("Select valid compounds.")
  if(type=="areas")return(invisible(matplot(compounds,t(x$areas[,compounds,drop=FALSE]),type="b",xlab="Compound",ylab="Integrated area",...)))
  axis<-suppressWarnings(as.numeric(x$mass));if(anyNA(axis))axis<-seq_along(x$mass)
  if(type=="spectra")return(invisible(matplot(axis,t(x$P[compounds,,drop=FALSE]),type="h",xlab="Channel",ylab="Spectral fraction",...)))
  samples<-if(is.null(sample))seq_along(x$paths)else cpl_sample_index(x,sample)
  old<-par(mfrow=c(length(compounds),1L));on.exit(par(old))
  colour<-grDevices::hcl.colors(length(x$paths),"Dark 3")
  for(g in compounds) {
    domain<-range(unlist(lapply(samples,function(s) {
      f<-readRDS(x$paths[s]);k<-which(x$map[f$gid]==g)
      if(!length(k))return(numeric())
      h<-as.numeric(f$H[,k]);f$aligned_time[h>max(h)*1e-4]
    })))
    margin<-max(diff(domain)*.1,.Machine$double.eps)
    plot(domain+c(-margin,margin),c(0,1.05),type="n",xlab="Aligned time",ylab="Relative profile",main=paste("Compound",g),...)
    for(s in samples) {
      f<-readRDS(x$paths[s]);k<-which(x$map[f$gid]==g)
      if(length(k)){h<-as.numeric(f$H[,k]);lines(f$aligned_time,h/max(h),col=colour[s],lwd=1.5)}
    }
    legend("topleft",legend=x$store$manifest$sample_id[samples],col=colour[samples],lty=1,lwd=1.5,bty="n",cex=.8)
  }
  invisible(x)
}

cpl_noise_pass <- function(state,samples) {
  for(s in samples) {
    f<-readRDS(state$paths[s])
    if(f$kind!=3L||!isTRUE(f$spec$theta_update))next
    Y<-cpl_sample_data(state$store,s);mu<-cpl_coupled_mean(f,state)
    for(iteration in 1:3) {
      score<-slcn_nb_score(Y,mu,f$theta,as.numeric(f$weights),f$ncores)
      move<-if(is.finite(score[2L])&&score[2L]<0)-score[1L]/score[2L]else sign(score[1L])
      move<-max(-1,min(1,move));theta0<-f$theta;reference<-cpl_coupled_ll(f,Y,mu)
      for(half in 0:12) {
        f$theta<-max(1e-6,min(1e8,theta0*exp(move/2^half)))
        if(cpl_coupled_ll(f,Y,mu)>=reference)break
        f$theta<-theta0
      }
      if(abs(log(f$theta/theta0))<1e-4)break
    }
    f$spec$theta<-f$theta
    state<-cpl_coupled_checkpoint(state,s,f,"noise")
  }
  state
}

#' Extract one sample's factorization from a coupled fit
#'
#' Reconstruct a standard [splogconNMF()] object on a sample's original grid.
#' Its profiles include alignment and exposure, and have unit column sums;
#' the matching spectra carry their original-grid amplitudes. `compound_id`
#' maps its local components to the coupled consensus compounds.
#' @param object A [coupledNMF()] fit.
#' @param sample One sample index or name.
#' @return A `splogconNMF` object without an attached observed or fitted matrix.
#' @export
sampleNMF <- function(object,sample=1L) {
  s<-cpl_sample_index(object,sample);f<-readRDS(object$paths[s]);ids<-object$map[f$gid]
  H<-f$H*f$exposure;sums<-as.numeric(colSums(H))
  S<-object$P[ids,,drop=FALSE]*(f$area*object$scale[ids]*sums)
  if(length(sums))H<-H%*%Diagonal(x=1/sums)
  model<-list(H=as(H,"dgCMatrix"),S=S,background=f$background*f$raw_exposure,alpha=f$alpha,
    offset=f$offset,weights=f$weights,spec=f$spec,n=length(f$time),p=length(f$time),
    ncores=f$ncores,rank=length(ids),compound_id=ids,integrated_area=f$area*object$scale[ids],
    time=f$time,control=list(censor=f$censor),seconds=object$seconds)
  model$loglik<-cpl_coupled_ll(f,cpl_sample_data(object$store,s),cpl_coupled_mean(f,object))
  model$nuisance_df<-as.integer(f$spec$dispersion_estimated)+as.integer(f$spec$theta_estimated)
  model$ic_df<-ncol(S)*nrow(S)+model$nuisance_df
  model$lambda<-object$lambda*ncol(S);model$objective<- -2*model$loglik+model$lambda*model$rank
  model<-slcn_temporal_order(model)
  # slcn_temporal_order orders H/S; carry the global IDs through the same ordering.
  apex<-vapply(seq_len(ncol(f$H)),function(k)which.max(f$H[,k]*f$exposure),0L)
  order<-order(apex);model$compound_id<-ids[order];model$integrated_area<-model$integrated_area[order]
  class(model)<-"splogconNMF"
  model
}


# Scalar and row-wise inputs are reused across channels without materialization.
cpl_observation <- function(x,n,m,default=1) {
  if(is.null(x))return(default)
  if(length(x)==1L)return(as.numeric(x))
  if(length(x)==n&&(is.null(dim(x))||ncol(x)==1L))return(as.numeric(x))
  if(!identical(dim(x),c(n,m)))stop("Observation input needs length n or dimensions n x m.")
  as.matrix(x)
}
cpl_observation_block <- function(x,rows) {
  if(length(x)==1L)x else if(is.null(dim(x)))x[rows]else x[rows,,drop=FALSE]
}
cpl_weight_matrix <- function(x,n) {
  if(length(x)==1L)matrix(x,n,1L)else if(is.null(dim(x)))matrix(x,n,1L)else x
}
