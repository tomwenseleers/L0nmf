# Simulate related chromatograms with shared compound spectra

Generate exponentially modified Gaussian and optional triangular acid
peaks, smooth backgrounds, irregular scan grids, sample drift and
concentration- dependent peak widths and tails. Missing compounds can be
simulated explicitly.

## Usage

``` r
simulate_GCMS_samples(
  samples = 3L,
  n = 600L,
  k = 26L,
  S = NULL,
  preset = c("time_area_acid", "same", "time", "time_area"),
  presence = 1,
  absent = NULL,
  irregular = 0.12,
  drift = 2,
  area = 2e+05,
  seed = 811,
  width_time = 0.5,
  width_area = 0.45,
  tail_time = 1,
  tail_area = 1.5,
  rt_area = 0.25,
  acid_fraction = 0.08,
  family = c("poisson", "quasipoisson", "nb", "gaussian", "Gamma", "nb1", "Gamma_scale"),
  dispersion = 1,
  theta = 20
)
```

## Arguments

- samples:

  Number of samples.

- n:

  Number of scans, recycled across samples.

- k:

  Number of distinct compounds.

- S:

  Optional spectral library with one compound per row.

- preset:

  Peak-shape scenario: identical shapes, time-dependent shapes, time and
  concentration-dependent shapes, or the latter with acid peaks.

- presence:

  Probability that each compound occurs in each sample.

- absent:

  Optional logical compounds-by-samples absence matrix.

- irregular:

  Relative scan-width variation, from zero to 0.5.

- drift:

  Drift per sample, in seconds.

- area:

  Typical integrated compound intensity.

- seed:

  Random-number seed.

- width_time, width_area, tail_time, tail_area, rt_area:

  Effects of time and concentration on peak width, tail and retention
  time.

- acid_fraction:

  Fraction of broad triangular acid profiles.

- family:

  Signal distribution name supported by
  [`coupledNMF()`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md).

- dispersion, theta:

  Noise parameters.

## Value

Simulation truth and a list of samples with observations, original bin
edges and times (minutes), integrated-count means, true spectra,
profiles, areas, presence indicators and acquisition dates. A simulated
observation is integrated over its bin; `raw_exposure` is the original
bin duration in seconds. For
[`coupledNMF()`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md),
supply it as exposure.

## See also

[`coupledNMF()`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md),
[`simulate_GCMS()`](https://tomwenseleers.github.io/L0nmf/reference/simulate_GCMS.md)

## Examples

``` r
sim <- simulate_GCMS_samples(samples=3L,n=c(180L,190L,200L),k=8L,
  presence=0.5,seed=3)
sim$present
```
