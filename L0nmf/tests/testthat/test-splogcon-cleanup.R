test_that("original-channel cleanup deletes, adds and records complete model moves", {
  set.seed(812);n <- 100L;m <- 5L
  X <- convmatrix(dnorm(-6:6, 0, 1.5), n)
  H <- X[, c(25L, 72L)];H <- H %*% Matrix::Diagonal(x = 1/Matrix::colSums(H))
  S <- matrix(runif(2*m, 1500, 3000), 2L, m)
  Y <- matrix(rpois(n*m, as.matrix(H %*% S) + 2), n, m)
  # One genuine component is initially absent; one false component is present.
  start <- list(H = X[, c(25L, 48L)], S = rbind(S[1, ], rep(.01, m)), alpha = matrix(2, 1L, m))
  f <- splogconNMF(X, Y, background = matrix(1, n, 1), start = start,
    ncores = 1L, correlation = "none", control = list(profile_cycles = 0L,
      profile_fusion = FALSE, broad_repeats = 0L, full_sweeps = 0L,
      final_polishes = 0L, original_cleanup = FALSE))
  z <- slcn_cleanup(f, Y, X)
  expect_lt(z$objective, f$objective)
  expect_true(any(abs(z$apex_scan - 72L) <= 2L))
  expect_false(any(abs(z$apex_scan - 48L) <= 2L))
  expect_true(any(z$cleanup$moves$original_residual$accepted))
  expect_equal(dim(z$S), c(z$rank, m))
  expect_equal(unname(z$mu), unname(as.matrix(z$H %*% z$S) + z$background %*% z$alpha), tolerance = 1e-8)
  expect_lt(abs(z$ledger_error), 1e-4)
  expect_error(splogconNMF(X, Y, control = list(original_cleanup = NA)), "original_cleanup")
})

test_that("cached deletion scores equal complete rescoring with matrix weights", {
  set.seed(813);n <- 75L;m <- 4L
  X <- convmatrix(dnorm(-5:5, 0, 1.2), n)
  W <- matrix(runif(n*m, .2, 1.2), n, m);W[1:4, ] <- 0
  f <- splogconNMF(X, matrix(0, n, m), weights = W, offset = 1,
    correlation = "none", ncores = 1L, control = list(original_cleanup = FALSE))
  f$H <- X[, c(20L, 23L, 55L)];f$H <- f$H %*% Matrix::Diagonal(x = 1/Matrix::colSums(f$H))
  f$S <- matrix(.02, 3L, m);f$mu <- as.matrix(f$H %*% f$S) + 1
  f$members <- as.list(1:3)
  z <- slcn_cleanup_prune(f, matrix(0, n, m))
  expect_equal(ncol(z$fit$H), 0L)
  expect_equal(z$fit$mu, matrix(1, n, m))
  expected <- slcn_cleanup_delta(f, matrix(0,n,m), f$mu, z$fit$mu, seq_len(n)) - 3*f$lambda
  expect_equal(sum(z$audit$delta), expected, tolerance = 1e-10)
})
