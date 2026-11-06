test_that("PSF generators and shifts provide normalized nonnegative kernels", {
  for (type in c("gaussian", "bigaussian", "emg", "lorentzian", "laplace",
                 "exponential", "gamma", "gammaexp")) {
    kernel <- convKernel(type, 51L)
    expect_length(kernel, 51L)
    expect_true(all(is.finite(kernel) & kernel >= 0))
    expect_equal(sum(kernel), 1, tolerance = 1e-12)
    expect_equal(shiftKernel(kernel, perc = 0), kernel, tolerance = 1e-12)
  }
})

test_that("convolution matrices reproduce direct full and same convolutions", {
  set.seed(192)
  x <- rnorm(40)
  kernels <- list(convKernel("gaussian", 15L),
    convKernel("bigaussian", 15L, list(sigma1 = 1, sigma2 = 3)),
    convKernel("emg", 15L, list(sigma = 1.2, tau = 3)))
  for (kernel in kernels) for (shape in c("same", "full")) {
    design <- convMatrix(kernel, length(x), shape = shape)
    expect_s4_class(design, "dgCMatrix")
    expect_equal(as.numeric(design %*% x),
      as.numeric(gsignal::conv(x, kernel, shape = shape)), tolerance = 1e-12)
  }
})

test_that("asymmetric PSF estimation recovers a bigaussian from a Gaussian start", {
  set.seed(12)
  n <- 500L; L <- 51L
  kernel <- convKernel("bigaussian", L, list(sigma1 = 1.25, sigma2 = 5), perc = "mean")
  A <- matrix(0, n, 3L); A[c(60, 180, 320, 440), ] <- 30000
  Y <- matrix(rpois(n * 3L, as.matrix(convMatrix(kernel, n) %*% A) + 10), n, 3L)
  initial <- convKernel("gaussian", L, list(sigma = 4))
  fit <- estimatePSF(Y, family = poisson("identity"), standard = NULL,
    method = "blind", precluster = 1L, offset = 10, ncores = 2L,
    control = list(kernel_est = initial, maxiter = 6L))
  cosine <- function(a, b) sum(a * b) / sqrt(sum(a * a) * sum(b * b))
  expect_gt(cosine(fit$g_hat, kernel), .99)
  expect_gt(cosine(fit$g_hat, kernel), cosine(initial, kernel))
  expect_equal(sum(fit$g_hat), 1, tolerance = 1e-12)
})
