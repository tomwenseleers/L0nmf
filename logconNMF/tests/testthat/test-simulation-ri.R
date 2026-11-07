test_that("RI placement, drift and concentration-dependent leading edges are coherent", {
  sim <- simulateGCMSSamples(samples=6L,n=500L,k=10L,presence=1,
    preset="time_area_acid",acid_fraction=.1,irregular=0,drift=.5,rt_area=0,
    retentionMap=data.frame(ri=c(650,1900),time=c(.65,7.6)),seed=17)
  expect_length(sim$samples[[1L]]$acid,1L)
  expect_equal(sim$samples[[1L]]$centers,
    splinefun(sim$retentionMap$ri,sim$retentionMap$time,method="monoH.FC")(
      sim$samples[[1L]]$retentionIndex))
  e <- new.env()
  data("simulatedSpectra",package="L0glm",envir=e)
  expect_equal(sim$samples[[1L]]$retentionIndex,
    unname(attr(e$simulatedSpectra,"RI")[sim$compound_ids]))
  onsets <- vapply(sim$samples,function(s)s$onset,numeric(10))
  expect_lt(max(apply(onsets,1,function(x)diff(range(x)))),1e-12)
  for (s in seq_len(5L)) {
    advance <- 60*(sim$samples[[s+1L]]$retentionCalibration$time-
      sim$samples[[s]]$retentionCalibration$time)
    expect_true(all(advance>0))
    expect_lt(max(advance),.7)
  }
  for (s in sim$samples)
    expect_equal(colSums(s$H_true*s$exposure),rep(1,10),tolerance=1e-12)
  acid <- sim$samples[[1L]]$acid
  widths <- vapply(sim$samples,function(s) {
    h<-s$H_true[,acid]
    diff(range(s$reference_time[h>=.01*max(h)]))*60
  },0)
  areas <- vapply(sim$samples,function(s)s$area_true[acid],0)
  expect_gt(cor(widths,areas),.7)
})

test_that("RI simulation rejects incoherent calibration metadata", {
  map <- data.frame(ri=c(650,1900),time=c(.6,7.8))
  expect_error(simulateGCMSSamples(n=500L,retentionMap=map[2:1,]),
    "increasing ri/time")
  expect_error(simulateGCMSSamples(n=500L,retentionMap=map,retentionIndex=1:3),
    "one finite RI per library row")
  expect_error(simulateGCMSSamples(n=500L,
    retentionMap=data.frame(ri=c(650,1900),time=c(20,30))),
    "inside the scan window")
})

test_that("supplied baseline intensities retain their scale and scan grids", {
  S<-matrix(c(1,0,0,1),2,2,dimnames=list(NULL,c("43","44")))
  b<-rbind(c(1,3),c(4,6))
  sim<-simulateGCMSSamples(samples=2L,n=c(80L,100L),k=2L,S=S,
    background=b,presence=0,irregular=0,seed=13)
  for(s in sim$samples) {
    expected<-cbind(1+3*s$reference_time*60/100,3+3*s$reference_time*60/100)
    expect_equal(s$background/s$exposure,expected)
    expect_equal(unname(s$mean),unname(s$background))
  }
  flat<-simulateGCMSSamples(samples=1L,n=80L,k=2L,S=S,background=c(2,5))
  expect_equal(flat$samples[[1L]]$background/flat$samples[[1L]]$exposure,
    matrix(c(2,5),80L,2L,byrow=TRUE))
  expect_error(simulateGCMSSamples(k=2L,S=S,background=-1),"nonnegative")
})
