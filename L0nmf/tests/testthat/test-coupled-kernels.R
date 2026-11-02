test_that("irregular profile and pooled spectrum updates use the original likelihood", {
  set.seed(88);n<-41L;m<-3L;t<-cumsum(runif(n,.7,1.3));e<-runif(n,.4,1.6)
  h<-dnorm(t,mean(t),4);s<-c(1600,2300,700);mu<-2+e*h%o%s
  Y<-matrix(rpois(length(mu),mu),n,m);W<-matrix(runif(n*m,.5,1.5),n,m)
  for(kind in 0:6) {
    phi<-if(kind==4L).3 else 3
    run<-function(nt)slcn_log_profile(Y,mu,h,s,kind,phi,10,as.numeric(W),n,
      q=24L,nthreads=nt,time=t,exposure=e)
    a<-run(1L)
    expect_identical(a,run(2L))
    expect_lt(max(diff(diff(log(a$profile))/diff(t))),1e-6)
    expect_true(all(a$history<=1e-7))
    z<-slcn_pooled_terms(Y,mu,e*h,kind,phi,10,W)
    expect_identical(z,slcn_pooled_terms(Y,mu,e*h,kind,phi,10,W,2L))
    eps<-1e-3;d<-e*h%o%c(0,eps,0)
    ll<-function(v)slcn_loglik(Y,v,kind,phi,10,as.numeric(W),1L,.5)
    expect_lt(abs(z$gradient[2]+(ll(mu+d)-ll(mu-d))/(2*eps)),1e-5)
  }
  a<-slcn_log_profile(Y,mu,h,s,1L,1,1,1,n,q=24L)
  b<-slcn_log_profile(Y,mu,h,s,1L,1,1,1,n,q=24L,time=seq_len(n),exposure=rep(1,n))
  expect_identical(a,b)
})
