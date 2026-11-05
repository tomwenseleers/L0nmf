# Package index

## Estimate point spread functions and deconvolve signals

- [`convkernel()`](https://tomwenseleers.github.io/L0nmf/reference/convkernel.md)
  [`shift_kernel()`](https://tomwenseleers.github.io/L0nmf/reference/convkernel.md)
  : Generate normalised 1-D convolution kernels (point-spread functions)
- [`convmatrix()`](https://tomwenseleers.github.io/L0nmf/reference/convmatrix.md)
  : Make banded or rectangular convolution design matrices
- [`estimate_psf()`](https://tomwenseleers.github.io/L0nmf/reference/estimate_psf.md)
  : Estimate constant or time-varying point spread functions
- [`estimate_psf_from_std()`](https://tomwenseleers.github.io/L0nmf/reference/estimate_psf_from_std.md)
  : Estimate a time-varying PSF matrix from an alkane ladder TIC
- [`L0glm_blind_deco()`](https://tomwenseleers.github.io/L0nmf/reference/L0glm_blind_deco.md)
  : Alternating sparse L0 glm powered blind deconvolution
- [`L0glm_sparse_deco()`](https://tomwenseleers.github.io/L0nmf/reference/L0glm_sparse_deco.md)
  : Sparse deconvolution with a fixed PSF via L0-penalised GLMs
- [`L0glm_deco_estkernel()`](https://tomwenseleers.github.io/L0nmf/reference/L0glm_deco_estkernel.md)
  : Estimate the PSF given a spike train estimate via L0-penalised GLMs

## Sparse and coupled factorization

- [`residuals(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
  [`splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
  [`fitted(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
  [`nobs(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
  [`coef(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
  [`logLik(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
  [`AIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
  [`BIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
  [`GIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
  [`EBIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
  [`MBIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
  [`print(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
  : Sparse log-concave nonnegative matrix factorization
- [`coupledNMF()`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md)
  [`fitted(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md)
  [`coef(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md)
  [`logLik(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md)
  [`nobs(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md)
  [`AIC(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md)
  [`BIC(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md)
  [`GIC(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md)
  [`EBIC(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md)
  [`MBIC(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md)
  [`print(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md)
  [`plot(`*`<coupledNMF>`*`)`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md)
  : Coupled sparse log-concave factorization of multiple samples
- [`sampleNMF()`](https://tomwenseleers.github.io/L0nmf/reference/sampleNMF.md)
  : Extract one sample's factorization from a coupled fit

## Noise calibration and profile constraints

- [`estimate_signal_noise()`](https://tomwenseleers.github.io/L0nmf/reference/estimate_signal_noise.md)
  : Estimate a shared noise model for ordered multichannel signals
- [`log_conc_smooth()`](https://tomwenseleers.github.io/L0nmf/reference/log_conc_smooth.md)
  : Log-concave spline smoothing of a discrete convolution kernel
  estimate

## Simulations and experimental data

- [`simulate_GCMS()`](https://tomwenseleers.github.io/L0nmf/reference/simulate_GCMS.md)
  : Simulate synthetic GC-MS chromatograms with ground-truth signal +
  background
- [`simulate_GCMS_samples()`](https://tomwenseleers.github.io/L0nmf/reference/simulate_GCMS_samples.md)
  : Simulate related chromatograms with shared compound spectra
- [`GCMS_data1`](https://tomwenseleers.github.io/L0nmf/reference/GCMS_data1.md)
  : GC/MS Chromatogram Section - Cuticular Hydrocarbon Profile of a
  Cataglyphis iberica Ant Queen
- [`GCMS_data2`](https://tomwenseleers.github.io/L0nmf/reference/GCMS_data2.md)
  : GC/MS Chromatogram Section - Beer Headspace Volatiles

## Shared tools and package overview

- [`reexports`](https://tomwenseleers.github.io/L0nmf/reference/reexports.md)
  [`GIC`](https://tomwenseleers.github.io/L0nmf/reference/reexports.md)
  [`EBIC`](https://tomwenseleers.github.io/L0nmf/reference/reexports.md)
  [`MBIC`](https://tomwenseleers.github.io/L0nmf/reference/reexports.md)
  : Objects exported from other packages
- [`L0nmf`](https://tomwenseleers.github.io/L0nmf/reference/L0nmf-package.md)
  [`L0nmf-package`](https://tomwenseleers.github.io/L0nmf/reference/L0nmf-package.md)
  : Scalable Sparse Nonnegative Matrix Factorization with Automatic Rank
  Selection and Flexible Noise Models
