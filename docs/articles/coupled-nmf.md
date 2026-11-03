# Coupled multisample log-concave NMF with automatic best rank tuning

[`coupledNMF()`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md)
combines related multichannel samples into a common spectral library
while fitting each sample’s own compound areas, peak shapes and smooth
background. Samples retain their measured time grids. The model
accommodates retention-time drift, concentration-dependent widths and
tails, and compounds that occur in only some samples.

``` r

fit <- coupledNMF(samples, family=quasipoisson("identity"), precluster=32L)
fit$areas                         # samples x consensus compounds
coef(fit, normalize="l1")          # common spectra
fitted(fit, sample=1L)             # original-grid fitted signal
sampleNMF(fit, sample=1L)          # standard per-sample NMF object
GIC(fit)
```

Each sample description contains `y`, `time`, and either a candidate PSF
matrix `X` or an independent
[`splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
fit. Matrices can be supplied through RDS paths. Observations, factors,
weights and offsets are stored per sample, and the optimization loop
reads one sample at a time. Common spectra and compact alignment
summaries are shared in memory.

## Three samples, shared spectra and different peak shapes

The simulator provides exponentially modified Gaussian peaks, broad
triangular acid profiles, smooth backgrounds and irregular scan
intervals. Peak width, tailing and retention time can depend on
concentration and position in the chromatogram. `presence=0.5` makes
compounds absent independently in half of the samples on average.

``` r

sim <- simulate_GCMS_samples(samples=3L, n=c(580L,600L,620L), k=26L,
  preset="time_area_acid", presence=0.5, seed=811)

samples <- lapply(sim$samples, function(s) list(
  y=s$Y, time=s$time, edges=s$edges, exposure=s$raw_exposure,
  sample_id=s$sample_id,
  X=convmatrix(dnorm(-8:8, sd=1.6), nrow(s$Y)),
  background=splines::bs(s$time, df=4L, intercept=TRUE)))

fit <- coupledNMF(samples, family=poisson("identity"), ncores=8L,
  control=list(initial=list(profile_cycles=4L, fusion_passes=2L)))
fit$areas
plot(fit, type="profiles", compounds=c(1L,12L))
```

The three-sample validation retains all **258 mass channels**. With
every compound present, it recovers **26 shared compounds and all 78
occurrences**. With `presence=0.5`, the realized data contain **23
shared compounds and 37 occurrences**, all recovered. Missing
occurrences remain `NA` in the area table.

| Simulation | True shared compounds | Recovered | True occurrences | Recovered | Total time |
|:---|---:|---:|---:|---:|---:|
| All compounds present | 26 | 26 | 78 | 78 | 29.69 s |
| 50% presence probability | 23 | 23 | 37 | 37 | 21.92 s |

These elapsed times include independent initialization, alignment,
grouping and coupled refinement on an eight-core Intel i9-11950H
workstation. Both examples include different scan grids and
concentration-dependent peak shapes.

![Aligned component profiles across the three
samples](figures/coupled-profiles.png)

Aligned component profiles across the three samples

![Recovered occurrence areas and absent
compounds](figures/coupled-areas.png)

Recovered occurrence areas and absent compounds

Shared spectra are estimated from all detected occurrences of a
compound:

![Shared spectra recovered from all sample
occurrences](figures/coupled-spectra.png)

Shared spectra recovered from all sample occurrences

## Original scan grids, shared spectra and missing compounds

Each sample retains its measured time grid, background and noise
parameters. Monotone alignment changes profile coordinates and bin
widths; observations are evaluated directly without resampling. Each
detected occurrence has its own log-concave shape and amplitude, while
its spectrum is shared.

The joint criterion counts one unit-sum spectrum per consensus component
and one amplitude per detected occurrence. Profile and background
complexity are unpenalized. Original-channel likelihood checks decide
whether a complete grouped/refined candidate is retained.

``` r

fit$areas
fit$compounds
fit$alignment
single <- sampleNMF(fit, sample=1L)
dim(single$H)
dim(coef(single))
plot(fit, type="areas")
```

An absent occurrence stays `NA`; this model does not fill undetected
compounds. The [coupled model and
algorithm](https://tomwenseleers.github.io/L0nmf/articles/algorithm.html#coupled-model)
give the exact likelihood, Jacobian, normalization, pooled spectral
updates and rank charge.

## Noise, exposure and disk-backed fitting

Gaussian, Poisson, quasipoisson, NB2, NB1, censored Gamma and
common-scale Gamma are supported. Noise parameters are shared over
channels within each sample. Use positive `dispersion`/`theta` to fix
known parameters; negative NB2 `theta` enables bounded
likelihood-improving updates from its absolute starting value.

Known observation exposure can be given as a scalar or scan vector in
each sample description. It multiplies the mean; the identity-link
`offset` adds known signal. The default unit exposure preserves the
units of supplied data. Alignment’s change in bin widths preserves
integrated occurrence areas.

Observation matrices and factor checkpoints are stored per sample on
disk. Only one sample’s observations are needed for a local update.
Shared spectra, areas and compact matching summaries stay in memory.
Supply `y` or independent `fit` components as RDS paths for file-backed
input; use a new `path=` for each fit’s persistent store. Keep that
store alongside saved fitted objects.

The native kernel pools spectral scores over sparse profile supports and
uses OpenMP over channels. On 14,949 scans, 518 channels and 325
profiles, one pooled score pass took 0.19 s; the [implementation
details](https://tomwenseleers.github.io/L0nmf/articles/algorithm.html#coupled-algorithm)
describe the update and threading structure.

For sample-group directories, acquisition ordering, alkane standards and
Excel/MSP reports, use the workflow in
[GCMSprime](https://tomwenseleers.github.io/GCMSprime/articles/multisample.html).
