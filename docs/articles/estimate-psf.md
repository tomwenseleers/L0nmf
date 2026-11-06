# Estimate PSFs

A short impulse is often recorded as a wider response: an optical blur,
an instrument response or a chromatographic peak. Its **point spread
function (PSF)** describes that response. With a known time-invariant
PSF, shifted copies form a convolution design matrix, and
[`l0glmFit()`](https://tomwenseleers.github.io/L0glm/reference/l0glmFit.html)
can estimate a sparse set of impulse locations and amplitudes directly.

When the PSF is unknown, alternate two conditional fits. Estimate the
spike train with the current kernel, then estimate the kernel with those
spikes held fixed.
[`estimatePSF()`](https://tomwenseleers.github.io/logconNMF/reference/estimatePSF.md)
can also calibrate the kernel from isolated standards and estimate how
it changes across overlapping windows.

The resulting sparse matrix `X` serves two purposes. Use it directly in
[`l0glmFit()`](https://tomwenseleers.github.io/L0glm/reference/l0glmFit.html)
when a common response shape at each position is adequate. Use it to
initialize [log-concave
NMF](https://tomwenseleers.github.io/logconNMF/articles/log-concave-nmf.md)
when components at the same position may have different widths or tails.

## Deconvolve with a known, constant PSF

This multichannel signal has eight impulses, an asymmetric bigaussian
PSF with left and right widths of 1.25 and 5 scans, and Poisson noise.
Channels share impulse locations but have different amplitudes.

``` r

library(L0glm)
library(logconNMF)
set.seed(611)
n <- 600L; m <- 6L
time <- seq_len(n) / 50
g <- convKernel("bigaussian", 51L, list(sigma1=1.25, sigma2=5),
  perc="mean")
X <- convMatrix(g, n)
A <- matrix(0, n, m)
A[c(65,125,190,255,335,405,475,540),] <- runif(8L*m, 1500, 3500)
Y <- matrix(rpois(n*m, as.matrix(X %*% A)+3), n, m)

known <- l0glmFit(X, Y, family=poisson("identity"), offset=3,
  nonnegative=TRUE, grouped=TRUE, relaxed=TRUE)
```

Selected rows in `coef(known)` locate the impulses. The [fixed-PSF
tutorial](https://tomwenseleers.github.io/L0glm/articles/deconvolution.html)
develops this regression model.

## Estimate an unknown kernel from overlapping signals

``` r

initial <- convKernel("gaussian", 51L, list(sigma=4))
blind <- estimatePSF(Y, time=time, family=poisson("identity"),
  standard=NULL, method="blind", nwindows=1L, precluster=1L,
  offset=3, control=list(kernel_est=initial, maxiter=6L))
blind$g_hat
blind$width_scans
```

`nwindows=1L` estimates one common PSF. The symmetric Gaussian start is
refined into an asymmetric response. The estimator starts with a
positive normalized kernel, alternates spike and kernel fits, and
retains the completed fit with the best requested information criterion.
Unit-sum normalization places scale in spike amplitudes; centroid
alignment fixes the shift convention.

The default `precluster=1L` uses the total across channels. `FALSE` uses
the original channels; an integer requests at most that many
intensity-adjusted groups. The returned `blind$X` is a sparse
scan-by-scan design ready for
[`l0glmFit()`](https://tomwenseleers.github.io/L0glm/reference/l0glmFit.html)
or
[`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md).

## Calibrate from isolated standard peaks

``` r

standard_spikes <- numeric(n)
standard_spikes[c(90,300,510)] <- c(12000,18000,15000)
standard <- list(Y=rpois(n, as.numeric(X %*% standard_spikes)+5), time=time)

calibrated <- estimatePSF(Y, time=time, standard=standard,
  family=poisson("identity"), method="standard", nwindows=1L,
  control=list(L=51L, baseline_df=4L, n_peaks=3L))
```

Calibration estimates a low-quantile background, identifies informative
isolated peaks and fits their log-concave profiles using the observation
likelihood. Normalized profiles are aligned by centroid and placed on
the sample scan interval. `method="standard"` retains this calibration;
`method="auto"` uses it to initialize blind refinement on the sample.

![Asymmetric bigaussian ground truth, symmetric initialization, blind
and calibrated PSFs, and window-specific
estimates](figures/psf-estimation.png)

Asymmetric bigaussian ground truth, symmetric initialization, blind and
calibrated PSFs, and window-specific estimates

## Allow the PSF to vary through the signal

``` r

varying <- estimatePSF(Y, time=time, standard=standard,
  family=poisson("identity"), method="auto", nwindows=3L,
  offset=3, control=list(L=51L, baseline_df=4L, n_peaks=3L, maxiter=6L))
varying$psf_timepoints
varying$psf_kernels
dim(varying$X)
```

Local kernels are estimated in overlapping windows and joined smoothly
on their log-profile scale. Each column of `varying$X` represents the
response appropriate to that impulse position. This synthetic signal has
a constant true kernel, so its local estimates should be similar; the
same construction captures gradual changes in instrumental response.

For an experimental example, bundled matrices include an alkane
calibration signal in `alkanes_TIC`, whose named time coordinate is
recognized automatically:

``` r

data(gcmsData2)
psf <- estimatePSF(gcmsData2, nwindows=3L, method="standard")
psf$method
psf$aggregation
```

## Reuse one standard across samples

Fit the total-ion chromatogram of an isolated-peak standard once, after
excluding the solvent region. Pass times explicitly for an unnamed
vector. Then use its calibration to build a sparse convolution design on
each sample’s scan grid:

``` r

keep <- standard_time >= 3  # minutes; choose the end of the solvent interval
calibration <- estimatePSF(standard_tic[keep], time=standard_time[keep],
  method="standard", nwindows=10L)
psf <- estimatePSF(sample_signal, time=sample_time, standard=calibration,
  method="standard", nwindows=10L)
X <- psf$X
```

The standard’s local peak shapes are reused. Each sample gets a design
on its own time grid and can have a different number of scans. Use
`method="auto"` if the sample also needs blind refinement from those
starting shapes.

Gaussian, Poisson, quasipoisson, NB2, NB1 and the two Gamma
parameterizations are supported with identity links. Noise scales are
shared across channels. Aggregation diagnostics record the variance
matching used in reduced models.

A position-dependent PSF describes a common instrumental response.
Components whose shapes also depend on their identity or intensity can
be modeled using [sparse log-concave
NMF](https://tomwenseleers.github.io/logconNMF/articles/log-concave-nmf.md).
The [mathematical
article](https://tomwenseleers.github.io/logconNMF/articles/algorithm.html#psf-estimation)
gives the objectives, initialization and pseudocode. File-based
calibration and automated reports are available in
[GCMSprime](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html).
