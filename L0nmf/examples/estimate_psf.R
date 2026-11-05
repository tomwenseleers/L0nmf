library(L0glm)
set.seed(12)
n <- 500L
time <- seq_len(n) / 60
L <- 51L
kernel <- convkernel("bigaussian", L,
  list(sigma1 = 1.25, sigma2 = 5), perc = "mean")
X <- convmatrix(kernel, n)
A <- matrix(0, n, 3L)
A[c(60, 180, 320, 440), ] <- 30000
Y <- matrix(rpois(n * 3L, as.matrix(X %*% A) + 10), n, 3L)
rownames(Y) <- time

# Calibrate from isolated peaks with the same asymmetric instrumental response.
standard_spikes <- numeric(n)
standard_spikes[c(80, 250, 420)] <- c(90000, 120000, 100000)
standard <- list(Y = rpois(n, as.numeric(X %*% standard_spikes) + 20), time = time)
calibration <- estimate_psf(Y, time = time, family = poisson("identity"),
  standard = standard, method = "standard", nwindows = 1L, ncores = 2L,
  control = list(L = L, baseline_df = 4L, n_peaks = 3L))
dim(calibration$X)

# Blind estimation starts from a symmetric Gaussian, using summed channels.
initial <- convkernel("gaussian", L, list(sigma = 4))
blind <- estimate_psf(Y, time = time, family = poisson("identity"),
  standard = NULL, method = "blind", precluster = 1L, offset = 10,
  ncores = 2L, control = list(kernel_est = initial, maxiter = 6L))
blind$width_scans

lag <- (seq_len(L) - (L + 1L) / 2) / 60
matplot(lag, cbind(kernel, initial, calibration$g_hat, blind$g_hat),
  type = "l", lty = c(1, 2, 1, 1), lwd = 2,
  col = c("black", "grey60", "#008b9a", "#cc3344"),
  xlab = "Relative retention time (min)", ylab = "Unit-area PSF")
legend("topright", c("True bigaussian", "Initial Gaussian", "Calibrated", "Blind"),
  lty = c(1, 2, 1, 1), lwd = 2,
  col = c("black", "grey60", "#008b9a", "#cc3344"), bty = "n")

\dontrun{
data(GCMS_data2)
Y <- GCMS_data2
# The alkane TIC and its time axis are taken from attr(Y, "alkanes_TIC").
psf <- estimate_psf(Y, nwindows = 3L)
B <- splines::bs(as.numeric(rownames(Y)), df = 4L, intercept = TRUE)
fit <- splogconNMF(psf$X, Y, background = B,
  family = quasipoisson("identity"), precluster = 32L)

# Calibration alone, or a single common PSF, use the same interface.
calibrated <- estimate_psf(Y, method = "standard", nwindows = 3L)
constant <- estimate_psf(Y, nwindows = 1L)
}
