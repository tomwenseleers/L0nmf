# Package index

## Estimate point spread functions and deconvolve signals

- [`convkernel()`](https://tomwenseleers.github.io/logconNMF/reference/convkernel.md)
  [`shift_kernel()`](https://tomwenseleers.github.io/logconNMF/reference/convkernel.md)
  : Generate normalised 1-D convolution kernels (point-spread functions)
- [`convmatrix()`](https://tomwenseleers.github.io/logconNMF/reference/convmatrix.md)
  : Make banded or rectangular convolution design matrices
- [`estimate_psf()`](https://tomwenseleers.github.io/logconNMF/reference/estimate_psf.md)
  : Estimate constant or time-varying point spread functions
- [`estimate_psf_from_std()`](https://tomwenseleers.github.io/logconNMF/reference/estimate_psf_from_std.md)
  : Estimate a time-varying PSF matrix from an alkane ladder TIC
- [`L0glm_blind_deco()`](https://tomwenseleers.github.io/logconNMF/reference/L0glm_blind_deco.md)
  : Alternating sparse L0 glm powered blind deconvolution
- [`L0glm_sparse_deco()`](https://tomwenseleers.github.io/logconNMF/reference/L0glm_sparse_deco.md)
  : Sparse deconvolution with a fixed PSF via L0-penalised GLMs
- [`L0glm_deco_estkernel()`](https://tomwenseleers.github.io/logconNMF/reference/L0glm_deco_estkernel.md)
  : Estimate the PSF given a spike train estimate via L0-penalised GLMs

## Sparse and coupled factorization

- [`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
  [`coef(`*`<logconNMF_fit>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
  [`fitted(`*`<logconNMF_fit>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
  : Fit a log-concave factorization with automatic rank selection
- [`residuals(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`splogconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`fitted(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`nobs(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`coef(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`logLik(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`AIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`BIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`GIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`EBIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`MBIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`print(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  : Sparse log-concave nonnegative matrix factorization
- [`coupledNMF()`](https://tomwenseleers.github.io/logconNMF/reference/coupledNMF.md)
  [`fitted(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/coupledNMF.md)
  [`coef(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/coupledNMF.md)
  [`logLik(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/coupledNMF.md)
  [`nobs(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/coupledNMF.md)
  [`AIC(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/coupledNMF.md)
  [`BIC(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/coupledNMF.md)
  [`GIC(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/coupledNMF.md)
  [`EBIC(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/coupledNMF.md)
  [`MBIC(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/coupledNMF.md)
  [`print(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/coupledNMF.md)
  [`plot(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/coupledNMF.md)
  : Coupled sparse log-concave factorization of multiple samples
- [`sampleNMF()`](https://tomwenseleers.github.io/logconNMF/reference/sampleNMF.md)
  : Extract one sample's factorization from a coupled fit
- [`channel_groups()`](https://tomwenseleers.github.io/logconNMF/reference/channel_groups.md)
  [`aggregate_channels()`](https://tomwenseleers.github.io/logconNMF/reference/channel_groups.md)
  : Learn and reuse a common channel aggregation

## Noise calibration and profile constraints

- [`estimate_signal_noise()`](https://tomwenseleers.github.io/logconNMF/reference/estimate_signal_noise.md)
  : Estimate a shared noise model for ordered multichannel signals
- [`log_conc_smooth()`](https://tomwenseleers.github.io/logconNMF/reference/log_conc_smooth.md)
  : Log-concave spline smoothing of a discrete convolution kernel
  estimate

## Simulations and experimental data

- [`simulate_GCMS()`](https://tomwenseleers.github.io/logconNMF/reference/simulate_GCMS.md)
  : Simulate synthetic GC-MS chromatograms with ground-truth signal +
  background
- [`simulate_GCMS_samples()`](https://tomwenseleers.github.io/logconNMF/reference/simulate_GCMS_samples.md)
  : Simulate related chromatograms with shared compound spectra
- [`GCMS_data1`](https://tomwenseleers.github.io/logconNMF/reference/GCMS_data1.md)
  : GC/MS Chromatogram Section - Cuticular Hydrocarbon Profile of a
  Cataglyphis iberica Ant Queen
- [`GCMS_data2`](https://tomwenseleers.github.io/logconNMF/reference/GCMS_data2.md)
  : GC/MS Chromatogram Section - Beer Headspace Volatiles

## Shared tools and package overview

- [`reexports`](https://tomwenseleers.github.io/logconNMF/reference/reexports.md)
  [`GIC`](https://tomwenseleers.github.io/logconNMF/reference/reexports.md)
  [`EBIC`](https://tomwenseleers.github.io/logconNMF/reference/reexports.md)
  [`MBIC`](https://tomwenseleers.github.io/logconNMF/reference/reexports.md)
  : Objects exported from other packages
- [`logconNMF-package`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF-package.md)
  : Log-Concave Nonnegative Matrix Factorization with Automatic Rank
  Selection
