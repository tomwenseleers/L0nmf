test_that("local fusion refreshes the entire support of every updated spectrum", {
  d <- polish_example()
  h <- numeric(120);at <- 90:120;h[at] <- dnorm(at,108,7);h <- h/sum(h)
  f <- d$f;f$H <- cbind(f$H,Matrix::Matrix(h,sparse=TRUE))
  f$S <- rbind(f$S,c(20000,500,15000,10000));f$members <- list(1L,2L,3L)
  Y <- as.matrix(f$H%*%f$S)+f$background%*%f$alpha+f$offset
  f <- slcn_sync(f,Y)
  Y[111:120,] <- Y[111:120,]*1.4
  z <- slcn_merge_anchor(f,Y,list(group=1:2,anchor=1L,shape="anchor"),
    cycles=0L,neighbor_mass=.01)
  expect_equal(z$mu,as.matrix(z$H%*%z$S)+z$background%*%z$alpha+z$offset,
    tolerance=1e-9)
  expect_gt(max(abs(z$S[1L,]-f$S[3L,])),1)
  expect_equal(dim(z$mu),dim(Y))
})

test_that("contained fusion respects all native losses and rejects distinct ion profiles", {
  for(kind in 0:6) {
    d <- polish_example(kind)
    z <- slcn_contained_broad(d$f,d$Y,psf_width=2L)
    expect_equal(dim(z$S),c(1L,4L),info=paste("kind",kind))
    expect_lte(z$objective,d$f$objective+1e-6)
    expect_equal(z$mu,as.matrix(z$H%*%z$S)+z$background%*%z$alpha+z$offset,
      tolerance=1e-9)
    expect_equal(as.numeric(Matrix::colSums(z$H)),1,tolerance=1e-9)
    h <- as.numeric(z$H);h <- h[h>0]
    expect_lte(max(diff(log(h),differences=2L)),1e-7)
    expect_identical(z$weights,d$f$weights)
    expect_identical(z$offset,d$f$offset)
    expect_identical(z$alpha,d$f$alpha)
    d <- polish_example(kind,distinct=TRUE,precision=1e5)
    z <- slcn_contained_broad(d$f,d$Y,psf_width=2L)
    expect_equal(ncol(z$H),2L,info=paste("distinct kind",kind))
    expect_lte(z$objective,d$f$objective+1e-6)
  }
})

test_that("initialization and final correlation-adjusted ledgers use their own criteria", {
  set.seed(347);n <- 70L;m <- 6L
  X <- convMatrix(convKernel("gaussian",11L,list(sigma=1.5)),n)
  H <- X[,c(20L,47L)];H <- H%*%Matrix::Diagonal(x=1/Matrix::colSums(H))
  S <- matrix(rexp(2*m,rate=1/2500),2,m)
  Y <- matrix(rpois(n*m,as.matrix(H%*%S)+2),n,m)
  z <- splogconNMF(X,Y,family=quasipoisson("identity"),background=matrix(1,n,1),
    dispersion=3,correlation=.2,precluster=2L,ncores=1L)
  expect_equal(z$lambda,m*1.4*presetLambda("gic",p=n,n=n))
  expect_equal(z$initialization$lambda,z$lambda/1.4)
  expect_true(all(z$history$delta[z$history$accepted]<=1e-6))
  expect_lt(abs(z$ledger_error),1e-5)
  expect_equal(z$mu,as.matrix(z$H%*%z$S)+z$background%*%z$alpha+z$offset,
    tolerance=1e-9)
  expect_equal(z$apex_scan,sort(z$apex_scan))
})
