# logconNMF 0.1.0.9000

## PSF estimation and deconvolution

* `estimatePSF()` estimates constant or smoothly time-varying point spread
  functions from standards or overlapping signals and returns a sparse design.
* Alternating spike-train and kernel estimation supports single or multiple
  channels, smooth backgrounds and information-criterion checkpoints.

## Sparse log-concave matrix factorization

* `logconNMF()` estimates nonnegative profiles and loadings with hard
  log-concavity constraints and automatic best rank tuning.
* Local profile, spectrum and background updates, algebraic candidate fusion,
  conditional deletion and residual additions use original-channel criteria.
* Gaussian, Poisson, quasipoisson, NB2, NB1, censored Gamma and common-scale
  Gamma models, with pooled noise calibration and optional temporal adjustment.
* Optional channel aggregation accelerates initialization; final refinement
  retains every original outcome channel.
* Permissive initialization and native fine-scale profile/spectrum refinement
  preserve candidate components for final correlation-aware rank selection.
* Contained broad-peak fragments are compared with equally refined unmerged
  controls. Complete local supports, reconstructed original-channel means
  and separate initialization/final-criterion ledgers make each decision auditable.

## Coupled multisample models

* `logconNMF()` with a sample list shares spectra across samples with sample-specific profiles,
  areas and backgrounds, allowing compounds to be absent from individual samples.
* Monotone alignment, area-ordered profile compatibility and exposure-aware
  likelihoods retain each sample's original irregular scan grid.
* Disk-backed observations, batched native spectral scores, shared row inputs,
  sparse profiles and OpenMP support efficient fitting.

## Simulation and interpretation

* Single- and multisample simulations include drifting retention times,
  concentration-dependent width/tailing, broad profiles and missing compounds.
* Experimental multichannel signals, profile and spectral plots, factor
  extraction, likelihood and information-criterion methods.
* Tutorials follow PSF estimation, single-sample factorization and coupled
  factorization; a separate mathematical article gives objectives and pseudocode.
