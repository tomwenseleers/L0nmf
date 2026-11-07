test_that("local starts preserve native objectives, constraints and nuisance terms", {
  for(kind in 0:6) {
    d <- polish_example(kind)
    h <- as.numeric(d$f$H[,1L]);at <- which(h>0);mode <- which.max(h)
    position <- mode+(seq_along(h)-mode)*1.7
    h <- exp(approx(at,log(h[at]),xout=position,yleft=-Inf,yright=-Inf)$y)
    d$f$H[,1L] <- h/sum(h)
    d$f <- slcn_sync(d$f,d$Y)
    z <- slcn_profile_starts(d$f,d$Y,cycles=4L,knots=32L)
    expect_lte(z$objective,d$f$objective+1e-6)
    expect_equal(ncol(z$H),ncol(d$f$H))
    expect_equal(as.numeric(Matrix::colSums(z$H)),rep(1,ncol(z$H)),tolerance=1e-9)
    for(k in seq_len(ncol(z$H))) {
      h <- as.numeric(z$H[,k]);h <- h[h>max(h)*1e-6]
      expect_lte(max(diff(log(h),differences=2L)),1e-6)
    }
    expect_equal(z$mu,as.matrix(z$H%*%z$S)+z$background%*%z$alpha+z$offset,
      tolerance=1e-9)
    expect_identical(z$spec,d$f$spec)
    expect_identical(z$weights,d$f$weights)
    expect_identical(z$offset,d$f$offset)
    expect_identical(z$alpha,d$f$alpha)
  }
})

test_that("the restart pass is enabled by default and validates its control", {
  n <- 35L
  X <- convMatrix(convKernel("gaussian",9L,list(sigma=1.5)),n)
  H <- X[,c(12L,25L)]
  Y <- as.matrix(H%*%rbind(c(120,30,2),c(4,30,100)))+2
  expect_error(splogconNMF(X,Y,family=gaussian(),dispersion=1,
    control=list(profile_restarts=1),ncores=1L),"profile_restarts")
  z <- splogconNMF(X,Y,family=gaussian(),dispersion=1,background=matrix(1,n,1),
    correlation="none",ncores=1L,control=list(fusion_passes=0L))
  expect_true(z$control$profile_restarts)
  expect_true("local_profile_restarts"%in%z$history$stage)
  expect_lt(abs(z$ledger_error),1e-5)
  expect_true(all(z$history$delta[z$history$accepted]<=1e-6))
  z <- splogconNMF(X,Y,family=gaussian(),dispersion=1,background=matrix(1,n,1),
    correlation="none",ncores=1L,control=list(profile_restarts=FALSE))
  expect_false("local_profile_restarts"%in%z$history$stage)
})

test_that("many-channel screening retains the original nuisance model", {
  d <- polish_example(2L)
  f <- d$f
  index <- rep(seq_len(ncol(d$Y)),10L)
  Y <- d$Y[,index]
  f$S <- f$S[,index];f$alpha <- f$alpha[,index]
  f$offset <- f$offset[,index];f$weights <- 1
  f <- slcn_sync(f,Y)
  z <- slcn_profile_starts(f,Y,cycles=2L,knots=32L)
  expect_equal(dim(z$S),dim(f$S))
  expect_lte(z$objective,f$objective+1e-6)
  expect_identical(z$offset,f$offset)
  expect_identical(z$alpha,f$alpha)
  expect_identical(z$spec,f$spec)
  expect_identical(z$weights,f$weights)
  expect_equal(z$mu,as.matrix(z$H%*%z$S)+z$background%*%z$alpha+z$offset,
    tolerance=1e-9)
})
