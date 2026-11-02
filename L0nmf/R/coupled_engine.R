cpl_coupled_mean <- function(f,state) {
  ids<-state$map[f$gid]
  S<-state$P[ids,,drop=FALSE]*(f$area*state$scale[ids])
  as.matrix(f$H%*%S)*f$exposure + f$background%*%f$alpha*f$raw_exposure + f$offset
}
cpl_coupled_ll <- function(f,Y,mu) {
  slcn_loglik(Y,mu,f$kind,f$phi,f$theta,as.numeric(f$weights),as.integer(f$ncores),f$censor)
}
cpl_coupled_score <- function(state,samples=seq_along(state$paths),penalty=TRUE) {
  ll<-0
  for(s in samples) {f<-readRDS(state$paths[s]);Y<-cpl_sample_data(state$store,s);ll<-ll+cpl_coupled_ll(f,Y,cpl_coupled_mean(f,state))}
  -2*ll + if(penalty)state$lambda*((ncol(state$P)-1L)*length(unique(state$map))+state$L+state$nuisance_df) else 0
}
cpl_coupled_checkpoint <- function(state,s,f,tag) {
  state$serial<-state$serial+1L
  state$paths[s]<-cpl_checkpoint_sample(f,state$store$manifest$directory[s],paste0(tag,"_",state$serial))
  state
}
cpl_coupled_profile_pass <- function(state,samples=seq_along(state$paths),roots=NULL,steps=4L,q=32L) {
  for(s in samples) {
    f<-readRDS(state$paths[s]);Y<-cpl_sample_data(state$store,s);mu<-cpl_coupled_mean(f,state)
    ids<-state$map[f$gid];block<-if(is.null(roots))seq_len(ncol(f$H)) else which(ids%in%roots)
    if(!length(block))next
    for(k in block[order(f$area[block]*state$scale[ids[block]],decreasing=TRUE)]) {
      at<-which(f$H[,k]>0);if(length(at)<3L)next
      h<-as.numeric(f$H[at,k]);sp<-state$P[ids[k],]*f$area[k]*state$scale[ids[k]]
      W<-if(length(f$weights)==1L)f$weights else f$weights[at,,drop=FALSE]
      z<-slcn_log_profile(Y[at,,drop=FALSE],mu[at,,drop=FALSE],h,sp,
        f$kind,f$phi,f$theta,as.numeric(W),length(at),censor=f$censor,
        maxit=as.integer(steps),q=as.integer(q),nthreads=as.integer(f$ncores),mode=2L,
        time=f$aligned_time[at],exposure=f$exposure[at])
      hnew<-z$profile;scale<-sum(hnew*f$quadrature[at])
      if(!is.finite(scale)||scale<=0)next
      nextmu<-mu[at,,drop=FALSE]+(f$exposure[at]*(hnew-h))%o%sp
      delta<-slcn_loglik(Y[at,,drop=FALSE],nextmu,
        f$kind,f$phi,f$theta,as.numeric(W),as.integer(f$ncores),f$censor)-
        slcn_loglik(Y[at,,drop=FALSE],mu[at,,drop=FALSE],f$kind,f$phi,f$theta,as.numeric(W),as.integer(f$ncores),f$censor)
      if(is.finite(delta)&&delta>= -1e-7) {
        f$H[at,k]<-hnew/scale;f$area[k]<-f$area[k]*scale;mu[at,]<-nextmu
      }
    }
    state<-cpl_coupled_checkpoint(state,s,f,"profile")
  }
  state
}

cpl_coupled_spectrum_pass <- function(state,samples=seq_along(state$paths),roots=NULL) {
  if(is.null(roots))roots<-sort(unique(state$map))
  if(!length(roots))return(state)
  m<-ncol(state$P);gradient<-information<-matrix(0,length(roots),m)
  for(s in samples) {
    f<-readRDS(state$paths[s]);Y<-cpl_sample_data(state$store,s);mu<-cpl_coupled_mean(f,state)
    ids<-state$map[f$gid]
    for(k in which(ids%in%roots)) {
      at<-which(f$H[,k]>0);h<-as.numeric(f$H[at,k])*f$exposure[at]*f$area[k]*state$scale[ids[k]]
      W<-if(length(f$weights)==1L)matrix(f$weights,length(at),m) else f$weights[at,,drop=FALSE]
      z<-slcn_pooled_terms(Y[at,,drop=FALSE],mu[at,,drop=FALSE],h,
        f$kind,f$phi,f$theta,W,as.integer(f$ncores),f$censor)
      j<-match(ids[k],roots);gradient[j,]<-gradient[j,]+z$gradient;information[j,]<-information[j,]+z$information
    }
  }
  old<-cpl_coupled_score(state,samples,penalty=FALSE)
  floor<-pmax(apply(information,1,max)*1e-10,1e-20)
  direction<- -gradient/pmax(information,floor)
  for(half in 0:20) {
    nextstate<-state;u<-pmax(state$P[roots,,drop=FALSE]+direction/2^half,0);scale<-rowSums(u)
    if(any(scale<=0))next
    nextstate$P[roots,]<-u/scale;nextstate$scale[roots]<-state$scale[roots]*scale
    score<-cpl_coupled_score(nextstate,samples,penalty=FALSE)
    if(is.finite(score)&&score<=old+1e-7)return(nextstate)
  }
  state
}
cpl_coupled_background_pass <- function(state,samples=seq_along(state$paths)) {
  for(s in samples) {
    f<-readRDS(state$paths[s]);if(!ncol(f$background))next
    Y<-cpl_sample_data(state$store,s);mu<-cpl_coupled_mean(f,state)
    B<-f$background*f$raw_exposure;fixed<-pmax(mu-B%*%f$alpha,0)
    W<-if(length(f$weights)==1L)matrix(f$weights,nrow(Y),ncol(Y)) else f$weights
    z<-if(f$kind==0L)nuisance_gaussian_cpp(B,Y,fixed,W/f$phi,f$alpha,
      rep(1L,ncol(B)),as.integer(f$ncores)) else if(f$kind%in%1:2)
      nuisance_poisson_cpp(B,Y,fixed,W/f$phi,f$alpha,80L,1e-8,as.integer(f$ncores)) else NULL
    if(is.null(z)) {
      z<-slcn_spectral(Y,mu,as(Matrix(B,sparse=TRUE),"dgCMatrix"),t(f$alpha),
        f$kind,f$phi,f$theta,as.numeric(W),nrow(Y),0L,4L,as.integer(f$ncores),20L,1e-8,f$censor)
      z<-list(coefficients=t(z$S),background=z$mu-fixed)
    }
    if(cpl_coupled_ll(f,Y,fixed+z$background)>=cpl_coupled_ll(f,Y,mu)-1e-7)f$alpha<-z$coefficients
    state<-cpl_coupled_checkpoint(state,s,f,"background")
  }
  state
}
cpl_optimize_coupled <- function(state,cycles=4L,samples=seq_along(state$paths),roots=NULL,record=TRUE,
                             q=32L,joint_profiles=FALSE) {
  previous<-cpl_coupled_score(state,samples)
  for(cycle in seq_len(cycles)) {
    start<-proc.time()[["elapsed"]]
    state<-cpl_coupled_profile_pass(state,samples,if(joint_profiles)NULL else roots,q=q)
    for(i in 1:3)state<-cpl_coupled_spectrum_pass(state,samples,roots)
    state<-cpl_coupled_background_pass(state,samples)
    state<-cpl_noise_pass(state,samples)
    score<-cpl_coupled_score(state,samples)
    if(score>previous+1e-4)stop("Joint objective increased.")
    if(record)state$history<-rbind(state$history,data.frame(stage="coupled_cycle",cycle,
      objective=score,delta=score-previous,seconds=proc.time()[["elapsed"]]-start))
    previous<-score
  }
  state
}

cpl_coupled_occurrences <- function(state) {
  out<-list()
  for(s in seq_along(state$paths)) {
    f<-readRDS(state$paths[s]);ids<-state$map[f$gid]
    geom<-cpl_profile_anchors(f$H,state$P[ids,,drop=FALSE]*(f$area*state$scale[ids]),
      f$aligned_time,f$small_width)
    out[[s]]<-data.frame(sample=rep.int(s,length(ids)),local=seq_along(ids),gid=f$gid,root=ids,geom)
  }
  do.call(rbind,out)
}
cpl_coupled_proposals <- function(state,method="cosine",rt_tolerance=state$control$rt_tolerance) {
  occ<-cpl_coupled_occurrences(state);groups<-split(seq_len(nrow(occ)),occ$root);ids<-as.integer(names(groups))
  center<-vapply(groups,function(i)weighted.mean(occ$anchor[i],sqrt(occ$area[i])),0)
  proposals<-list()
  # RT bins bound comparisons; no all-occurrence pairwise similarity matrix.
  sorted<-order(center)
  for(aa in seq_along(sorted)) {
    a<-sorted[aa];later<-sorted[which(seq_along(sorted)>aa & center[sorted]<=center[a]+rt_tolerance)]
    for(b in later) {
      if(length(intersect(occ$sample[groups[[a]]],occ$sample[groups[[b]]])))next
      distance<-as.numeric(cpl_spectral_distance(state$P[ids[a],,drop=FALSE],state$P[ids[b],,drop=FALSE],method))
      if(distance > if(method=="cosine") .035 else .18)next
      proposals[[length(proposals)+1L]]<-data.frame(a=ids[a],b=ids[b],distance=distance,rt_gap=abs(center[a]-center[b]))
    }
  }
  if(!length(proposals))return(data.frame())
  out<-do.call(rbind,proposals);out[order(out$distance,out$rt_gap,out$a,out$b),]
}

# Algebraic clustering performs no GLM or NMF refits. One joint fit follows.
cpl_batch_spectral_groups <- function(state,method="cosine",ridge=TRUE) {
  begun<-proc.time()[["elapsed"]];before<-cpl_coupled_score(state)
  # Put shared scaling into occurrence amplitudes before pooling spectra.
  for(s in seq_along(state$paths)) {
    f<-readRDS(state$paths[s]);f$area<-f$area*state$scale[state$map[f$gid]]
    state<-cpl_coupled_checkpoint(state,s,f,"batch_scale")
  }
  state$scale[]<-1;rejected<-character();audit<-list()
  repeat {
    proposals<-cpl_coupled_proposals(state,method)
    if(!nrow(proposals))break
    keys<-paste(proposals$a,proposals$b,sep=";");proposals<-proposals[!keys%in%rejected,,drop=FALSE]
    if(!nrow(proposals))break
    p<-proposals[1L,];occ<-cpl_coupled_occurrences(state);hit<-which(occ$root%in%c(p$a,p$b))
    profiles<-times<-list();areas<-numeric(length(hit))
    for(j in seq_along(hit)) {
      v<-occ[hit[j],];f<-readRDS(state$paths[v$sample]);profiles[[j]]<-as.numeric(f$H[,v$local])
      times[[j]]<-f$aligned_time;areas[j]<-v$area
    }
    check<-if(ridge)cpl_ridge_compatibility(profiles,times,areas,max_shift=state$control$rt_tolerance)else list(pass=TRUE,coverage=1)
    audit[[length(audit)+1L]]<-data.frame(a=p$a,b=p$b,distance=p$distance,
      coverage=min(check$coverage),accepted=check$pass)
    if(check$pass) {
      weight<-vapply(c(p$a,p$b),function(g)sum(occ$area[occ$root==g]),0)
      state$P[p$a,]<-as.numeric(crossprod(weight,state$P[c(p$a,p$b),,drop=FALSE]))/sum(weight)
      state$map[state$map==p$b]<-p$a
    } else rejected<-c(rejected,paste(p$a,p$b,sep=";"))
  }
  state$batch_grouping<-if(length(audit))do.call(rbind,audit)else data.frame()
  after<-cpl_coupled_score(state)
  state$history<-rbind(state$history,data.frame(stage="algebraic_grouping",cycle=0L,
    objective=after,delta=after-before,seconds=proc.time()[["elapsed"]]-begun))
  state
}

cpl_compact_coupled <- function(state) {
  occ<-cpl_coupled_occurrences(state);roots<-sort(unique(state$map))
  center<-vapply(roots,function(g)weighted.mean(occ$anchor[occ$root==g],sqrt(occ$area[occ$root==g])),0)
  roots<-roots[order(center)];state$root_lookup<-data.frame(compound=seq_along(roots),source_root=roots)
  state$P<-state$P[roots,,drop=FALSE];state$scale<-state$scale[roots]
  state$map<-match(state$map,roots)
  if(is.null(colnames(state$P))) {
    meta<-readRDS(file.path(state$store$manifest$directory[1L],"metadata.rds"));colnames(state$P)<-meta$mass
  }
  rownames(state$P)<-if(nrow(state$P))paste0("compound",seq_len(nrow(state$P)))else character()
  state$compounds<-cpl_coupled_occurrences(state)
  state$compounds<-state$compounds[order(state$compounds$root,state$compounds$sample),]
  state$areas<-matrix(NA_real_,length(state$paths),nrow(state$P),
    dimnames=list(state$store$manifest$sample_id,rownames(state$P)))
  state$areas[cbind(state$compounds$sample,state$compounds$root)]<-state$compounds$area
  state$consensus_counts<-state$P*colSums(state$areas,na.rm=TRUE)
  class(state)<-"coupledNMF"
  state
}

cpl_profile_anchors <- function(H,S,time,small_width) {
  K<-ncol(H);out<-data.frame(component=seq_len(K),area=rowSums(S),start=numeric(K),apex=numeric(K),width=numeric(K),anchor=numeric(K))
  for(k in seq_len(K)) {
    h<-as.numeric(H[,k]);at<-which(h>max(h)*.01);peak<-which.max(h)
    # Linear interpolation of the 1%-height crossing on the measured grid.
    first<-min(at);left<-if(first>1L) first-1L else first
    out$start[k]<-if(left<first)time[left]+(time[first]-time[left])*(.01*max(h)-h[left])/(h[first]-h[left]) else time[first]
    out$apex[k]<-time[peak];out$width[k]<-diff(range(time[at]))
    out$anchor[k]<-out$start[k]+rep_len(small_width,K)[k]/2
  }
  out
}
cpl_spectral_distance <- function(A,B,method=c("cosine","hellinger")) {
  method<-match.arg(method)
  if(method=="cosine") return(1-tcrossprod(A,B)/sqrt(rowSums(A^2)%o%rowSums(B^2)))
  A<-A/rowSums(A);B<-B/rowSums(B)
  sqrt(pmax(1-tcrossprod(sqrt(A),sqrt(B)),0))
}
cpl_monotone_map <- function(x,y,weight,domain) {
  o<-order(x);x<-x[o];y<-y[o];weight<-weight[o]
  live<-!duplicated(x);x<-x[live];y<-y[live];weight<-weight[live]
  if(length(x)<4L) {
    shift<-if(length(x))median(y-x) else 0
    return(list(evaluate=function(t)t+shift,derivative=function(t)rep(1,length(t)),
      anchors=data.frame(x,y,weight),method="shift"))
  }
  robust<-rep(1,length(x));fit<-NULL
  for(it in 1:5) {
    fit<-smooth.spline(x,y,w=pmax(1e-5,weight*robust),df=min(5L,length(x)-1L))
    residual<-y-predict(fit,x)$y;scale<-max(1e-5,mad(residual))
    u<-residual/(4.685*scale);robust<-ifelse(abs(u)<1,(1-u^2)^2,0)
  }
  grid<-sort(unique(c(domain,x,seq(min(domain),max(domain),length.out=100L))))
  values<-predict(fit,grid)$y
  # Project onto a strictly increasing map; monotone Hermite interpolation
  # then retains positivity between knots, without resampling any observations.
  minimum_slope<-.25
  values<-isoreg(grid,values-minimum_slope*grid)$yf+minimum_slope*grid
  fun<-splinefun(grid,values,method="monoH.FC")
  list(evaluate=fun,
    anchors=data.frame(x,y,weight,robust),method="robust_monotone_spline")
}
cpl_estimate_alignment <- function(fits,metadata,method="cosine",small_width=NULL,alignment_window=.5) {
  geom<-lapply(seq_along(fits),function(s) {
    width<-if(is.null(small_width))2*sqrt(2*log(100))*1.6*median(diff(metadata[[s]]$time)) else small_width[[s]]
    cpl_profile_anchors(fits[[s]]$H,fits[[s]]$S,metadata[[s]]$time,width)
  })
  ref<-which.max(vapply(fits,function(f)ncol(f$H),0L));maps<-anchors<-vector("list",length(fits))
  for(s in seq_along(fits)) {
    if(s==ref) {maps[[s]]<-list(evaluate=identity,method="reference");next}
    if(!ncol(fits[[s]]$H)||!ncol(fits[[ref]]$H)) {
      maps[[s]]<-cpl_monotone_map(numeric(),numeric(),numeric(),range(metadata[[s]]$edges));next
    }
    d<-cpl_spectral_distance(fits[[s]]$S,fits[[ref]]$S,method)
    gap<-abs(outer(geom[[s]]$anchor,geom[[ref]]$anchor,"-"))
    d[gap>alignment_window]<-Inf
    j<-apply(d,1L,which.min);i<-apply(d,2L,which.min)
    good<-which(i[j]==seq_len(nrow(d)) & d[cbind(seq_len(nrow(d)),j)] < if(method=="cosine") .025 else .16)
    x<-geom[[s]]$anchor[good];y<-geom[[ref]]$anchor[j[good]]
    w<-sqrt(pmin(geom[[s]]$area[good],geom[[ref]]$area[j[good]]))
    if(length(w))w<-pmin(w/median(w),4)
    ratio<-log(geom[[s]]$area[good]/geom[[ref]]$area[j[good]])
    maps[[s]]<-cpl_area_adjusted_map(x,y,w,ratio,range(metadata[[s]]$edges))
    anchors[[s]]<-data.frame(sample=good,reference=j[good],x,y)
  }
  list(maps=maps,anchors=anchors,geometry=geom,reference=ref,method=method)
}
cpl_project_irregular_lc <- function(h,t,min_relative=1e-12) {
  eta<-log(pmax(h,max(h)*min_relative));hull<-integer()
  for(i in seq_along(t)) {
    hull<-c(hull,i)
    while(length(hull)>2L) {
      z<-tail(hull,3L)
      if((eta[z[2]]-eta[z[1]])/(t[z[2]]-t[z[1]]) >= (eta[z[3]]-eta[z[2]])/(t[z[3]]-t[z[2]]))break
      hull<-hull[-(length(hull)-1L)]
    }
  }
  exp(approx(t[hull],eta[hull],xout=t,rule=2)$y)
}

# Compatibility diagnostic for an area-ordered ridge. The smaller profile
# should mostly fit underneath the larger one after the sample-level warp.
# Separate concavity is enforced by the optimizer; joint 2D concavity is not assumed.
cpl_ridge_compatibility <- function(profiles,times,areas,tolerance=.7,max_shift=.12) {
  if(length(profiles)<2L)return(list(pass=TRUE,coverage=1))
  o<-order(areas);coverage<-raw_coverage<-shift<-numeric(length(o)-1L)
  for(j in seq_along(coverage)) {
    a<-o[j];b<-o[j+1L];grid<-sort(unique(c(times[[a]],times[[b]])))
    low<-approx(times[[a]],profiles[[a]]*areas[a],grid,yleft=0,yright=0)$y
    high<-approx(times[[b]],profiles[[b]]*areas[b],grid,yleft=0,yright=0)$y
    dt<-diff(grid);integral<-function(x)sum(dt*(head(x,-1L)+tail(x,-1L))/2)
    raw_coverage[j]<-integral(pmin(low,high))/max(1e-20,integral(low))
    # Concentration can move the ridge as well as widen it. A bounded
    # leading-edge translation is used only for this shape diagnostic.
    # Observations, fitted profiles and exposure grids are never moved here.
    start_a<-times[[a]][which(profiles[[a]]>max(profiles[[a]])*.01)[1L]]
    start_b<-times[[b]][which(profiles[[b]]>max(profiles[[b]])*.01)[1L]]
    shift[j]<-max(-max_shift,min(max_shift,start_a-start_b))
    aligned_high<-approx(times[[b]]+shift[j],profiles[[b]]*areas[b],grid,yleft=0,yright=0)$y
    coverage[j]<-max(raw_coverage[j],integral(pmin(low,aligned_high))/max(1e-20,integral(low)))
  }
  list(pass=all(coverage>=tolerance),coverage=coverage,raw_coverage=raw_coverage,shift=shift)
}
# Robust low-df drift spline with an area-ratio nuisance term. The latter
# explains concentration-related anchor motion and is omitted from the warp.
cpl_area_adjusted_map <- function(x,y,weight,log_ratio,domain) {
  n<-length(x)
  if(n<4L)return(cpl_monotone_map(x,y,weight,domain))
  df<-min(5L,max(2L,floor(n/5L)))
  basis<-splines::ns(x,df=df-1L,Boundary.knots=domain)
  B<-cbind(1,basis);D<-cbind(B,log_ratio);target<-y-x
  robust<-rep(1,n);coefs<-numeric(ncol(D))
  for(it in 1:8) {
    W<-pmax(1e-6,weight*robust);G<-crossprod(D,D*W)
    diag(G)<-diag(G)+1e-8*max(diag(G))
    coefs<-as.numeric(solve(G,crossprod(D,W*target)))
    r<-target-as.numeric(D%*%coefs);scale<-max(1e-5,mad(r))
    robust<-pmin(1,1.345*scale/pmax(abs(r),1e-20))
  }
  grid<-sort(unique(c(domain,x,seq(min(domain),max(domain),length.out=100L))))
  values<-grid+as.numeric(cbind(1,predict(basis,grid))%*%coefs[seq_len(ncol(B))])
  values<-isoreg(grid,values-.25*grid)$yf+.25*grid
  fun<-splinefun(grid,values,method="monoH.FC")
  list(evaluate=fun,
    anchors=data.frame(x,y,weight,robust,log_ratio),area_coefficient=tail(coefs,1L),
    method="robust_spline_with_area_nuisance",df=df)
}

cpl_write_matrix_file <- function(Y,path) {
  con<-file(path,"wb");on.exit(close(con));writeBin(as.double(Y),con,size=8L,endian="little")
  invisible(path)
}
cpl_read_matrix_file <- function(path,n,m,channels=NULL,rows=NULL) {
  if(is.null(channels)) channels<-seq_len(m)
  if(is.null(rows)) rows<-seq_len(n)
  if(any(channels<1|channels>m)||any(rows<1|rows>n))stop("Invalid matrix slice.")
  con<-file(path,"rb");on.exit(close(con))
  if(identical(as.integer(channels),seq_len(m))&&identical(as.integer(rows),seq_len(n))) {
    x<-readBin(con,"double",n=as.double(n)*m,size=8L,endian="little")
    if(length(x)!=as.double(n)*m)stop("Truncated sample matrix.")
    return(matrix(x,n,m))
  }
  Y<-matrix(0,length(rows),length(channels))
  contiguous<-length(rows)>0L && all(diff(rows)==1L)
  for(j in seq_along(channels)) {
    start<-if(contiguous)rows[1L]-1L else 0L
    count<-if(contiguous)length(rows) else n
    seek(con,((channels[j]-1)*as.double(n)+start)*8,origin="start")
    z<-readBin(con,"double",n=count,size=8L,endian="little")
    if(length(z)!=count)stop("Truncated sample matrix.")
    Y[,j]<-if(contiguous)z else z[rows]
  }
  Y
}
cpl_sample_data <- function(store,s,channels=NULL,rows=NULL,aggregate=FALSE) {
  d<-store$manifest$directory[s];meta<-readRDS(file.path(d,"metadata.rds"))
  m<-if(aggregate)meta$aggregate_channels else meta$m
  cpl_read_matrix_file(file.path(d,if(aggregate)"aggregate.bin" else "observed.bin"),meta$n,m,channels,rows)
}
cpl_checkpoint_sample <- function(fit,directory,tag) {
  # Immutable checkpoint names avoid replacing a file while it is being read.
  path<-tempfile(pattern=paste0("fit_",tag,"_"),tmpdir=directory,fileext=".rds")
  saveRDS(fit,path,compress=FALSE)
  check<-readRDS(path);stopifnot(identical(dim(check$H),dim(fit$H)))
  normalizePath(path,winslash="/",mustWork=TRUE)
}
