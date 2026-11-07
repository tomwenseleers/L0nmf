# Alternative local starts for nearly coincident component modes.
# Completed restarts are accepted using the original-channel likelihood.
slcn_profile_starts <- function(f,Y,cycles=16L,knots=128L) {
  begun<-proc.time()[["elapsed"]]
  f<-slcn_sync(f,Y);K<-ncol(f$H)
  if(K<2L)return(f)
  geometry<-slcn_profile_geometry(f)
  sorted<-order(geometry$apex);eligible<-rep(FALSE,K)
  for(a in seq_len(K-1L)) {
    k<-sorted[a];b<-a+1L
    while(b<=K && geometry$apex[sorted[b]]-geometry$apex[k]<=geometry$width[k]/4) {
      j<-sorted[b]
      if(geometry$apex[j]-geometry$apex[k]<=geometry$width[j]/4)eligible[c(k,j)]<-TRUE
      b<-b+1L
    }
  }
  if(!any(eligible))return(f)
  audit<-list();sources<-f$members
  criterion<-function(z)slcn_loglik(Y,z$mu,z$spec$kind,z$spec$dispersion,z$spec$theta,
    as.numeric(z$weights),as.integer(z$ncores),z$control$censor)
  deform<-function(h,left,right) {
    at<-which(h>0);apex<-which.max(h);grid<-seq_along(h)
    position<-apex+(grid-apex)/ifelse(grid<apex,left,right)
    h<-exp(stats::approx(at,log(h[at]),xout=position,yleft=-Inf,yright=-Inf)$y)
    h/sum(h)
  }
  for(k in which(eligible)[order(rowSums(f$S)[eligible],decreasing=TRUE)]) {
    current<-criterion(f);h<-as.numeric(f$H[,k])
    if(sum(h>0)<3L)next
    geometry<-slcn_profile_geometry(f)
    trials<-lapply(list(c(2,1),c(1,2),c(2,2)),function(scale)deform(h,scale[1],scale[2]))
    at<-which(Reduce(function(x,y)x+y,trials)+h>max(h)*1e-6)
    block<-which(geometry$left<=max(at)&geometry$right>=min(at))
    screen<-function(z)slcn_refine_block(z,Y,block=block,cycles=1L,profile_steps=2L,
      adaptive=FALSE,joint=TRUE,pad=2L,free=block,knots=16L)
    control<-screen(f);control_score<-if(ncol(control$H)==K)criterion(control)else -Inf
    best<-f;best_score<-current;choice<-"none"
    for(j in seq_along(trials)) {
      z<-f;z$H[,k]<-trials[[j]]
      z$mu<-as.matrix(z$H%*%z$S)+z$background%*%z$alpha+z$offset
      z<-screen(z)
      score<-if(ncol(z$H)==K)criterion(z)else -Inf
      if(is.finite(score)&&score>max(best_score,control_score)+1e-6) {
        best<-z;best_score<-score;choice<-c("left","right","both")[j]
      }
    }
    accepted<-FALSE;delta<-0
    if(choice!="none") {
      polish<-function(z,iterations=cycles)slcn_refine_block(z,Y,block=block,cycles=iterations,
        profile_steps=6L,adaptive=FALSE,joint=TRUE,pad=2L,free=block,knots=knots)
      best<-polish(best);control<-polish(f,cycles+1L)
      score<-if(ncol(best$H)==K)criterion(best)else -Inf
      reference<-if(ncol(control$H)==K)criterion(control)else -Inf
      if(is.finite(score)&&score>max(current,reference)+1e-6) {
        delta<- -2*(score-current);f<-best;accepted<-TRUE
      }
    }
    audit[[length(audit)+1L]]<-data.frame(component=k,start=choice,
      screen_delta=-2*(best_score-current),delta=delta,accepted=accepted)
  }
  f<-slcn_sync(f,Y)
  f$profile_starts<-list(audit=if(length(audit))do.call(rbind,audit)else data.frame(),
    seconds=proc.time()[["elapsed"]]-begun,source_components=sources)
  f
}
