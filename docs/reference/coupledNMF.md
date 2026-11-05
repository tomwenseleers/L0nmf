# Coupled sparse log-concave factorization of multiple samples

Align related multichannel samples and estimate common spectra with
sample-specific nonnegative areas, log-concave profiles and backgrounds.
Observations remain on their original time grids. Samples are read from
a disk-backed store in bounded blocks during fitting. Independent
samples are processed concurrently within a single logical-core budget.

## Usage

``` r
coupledNMF(
  samples,
  family = poisson("identity"),
  lambda = "gic",
  path = tempfile("coupledNMF_"),
  precluster = 32L,
  dispersion = NULL,
  theta = NULL,
  method = c("cosine", "hellinger"),
  align = TRUE,
  ncores = NULL,
  control = list(),
  verbose = FALSE
)

# S3 method for class 'coupledNMF'
fitted(object, sample = 1L, ...)

# S3 method for class 'coupledNMF'
coef(object, sample = NULL, normalize = c("none", "l1"), ...)

# S3 method for class 'coupledNMF'
logLik(object, ...)

# S3 method for class 'coupledNMF'
nobs(object, ...)

# S3 method for class 'coupledNMF'
AIC(object, ..., k = 2)

# S3 method for class 'coupledNMF'
BIC(object, ...)

# S3 method for class 'coupledNMF'
GIC(object, ...)

# S3 method for class 'coupledNMF'
EBIC(object, ...)

# S3 method for class 'coupledNMF'
MBIC(object, ...)

# S3 method for class 'coupledNMF'
print(x, ...)

# S3 method for class 'coupledNMF'
plot(
  x,
  type = c("areas", "spectra", "profiles"),
  sample = NULL,
  compounds = seq_len(min(6L, nrow(x$P))),
  ...
)
```

## Arguments

- samples:

  List of sample descriptions. Each contains `y` (a numeric
  scans-by-channels matrix, or an RDS file containing one), increasing
  `time`, and optionally `X`, `background`, `fit`, `weights`, `offset`,
  `exposure`, `sample_id`, `group`, `small_width` and `edges`. Large
  `weights` and `offset` matrices can also be supplied as RDS paths.
  `fit` is a
  [`splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
  initializer or its RDS path. Without a fit, `X` supplies candidate
  profiles; otherwise
  [`estimate_psf()`](https://tomwenseleers.github.io/L0nmf/reference/estimate_psf.md)
  estimates them. Channel names and ordering must agree across samples.
  Times and `small_width` use the same units. Shared scan-wise weights
  and offsets remain vectors throughout coupled fitting; each mass
  channel reuses the same values.

- family:

  Identity-link signal family supported by
  [`splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md):
  Gaussian, Poisson, quasipoisson, NB2, NB1, censored Gamma or
  common-scale Gamma.

- lambda:

  Information-criterion preset or nonnegative penalty per selection
  degree of freedom. The selection count is
  `(channels - 1) * shared_compounds + detected_occurrences`, plus
  estimated noise parameters. Background and profile complexity are
  unpenalized.

- path:

  New directory for the sample store and fitted factors. Keep this
  directory with the saved R object; its paths are used by
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and plots.

- precluster:

  Maximum common outcome groups, default `32L`; `NULL` or `FALSE`
  disables aggregation. A saved
  [`channel_groups()`](https://tomwenseleers.github.io/L0nmf/reference/channel_groups.md)
  map can be reused. The default map is learned from samples summed in
  common time bins. Matching and Poisson/quasipoisson profile proposals
  use this common representation; likelihood checks and spectral
  estimation retain all original channels.

- dispersion, theta:

  Shared-across-channels noise parameters, as in
  [`splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md).
  They may be supplied per sample in its description. Resolved
  dispersions are frozen; a negative NB2 theta enables bounded
  likelihood-improving shape updates in the coupled iterations.

- method:

  Spectral distance for anchor matching and grouping.

- align:

  Estimate monotone sample-level time warps.

- ncores:

  Total logical-CPU budget, defaulting to physical cores. Up to
  `min(length(samples), ncores)` sample workers are used. The remaining
  budget determines native threads per worker; BLAS uses one thread
  during fitting. Caller BLAS settings are restored afterwards.

- control:

  Named controls: `cycles=10L`, `initial_cycles=5L`, `rounds=3L`,
  `q=64L`, `rt_tolerance=NULL`, `alignment_window=NULL`, `ridge=TRUE`,
  `censor=0.5`, `correlation="none"`, `workers=NULL`, `cache=TRUE`, and
  `initial=list()` for controls passed to independent
  [`splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
  fits. Default matching windows are derived from the supplied or
  estimated narrow-peak widths. `workers` limits concurrent samples
  within `ncores`. `cache=TRUE` holds only aggregated observations in
  RAM during fitting; `FALSE` reloads those blocks from disk.

- verbose:

  Print initialization and coupled objective progress.

- object, x:

  A fitted coupled model.

- sample:

  One sample index or name; `NULL` in
  [`coef()`](https://rdrr.io/r/stats/coef.html) returns consensus
  spectra, and in
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) selects all
  samples.

- ...:

  Additional method arguments.

- normalize:

  Return unit-sum spectra or amplitude-scaled spectra.

- k:

  Information-criterion multiplier for
  [`AIC()`](https://rdrr.io/r/stats/AIC.html).

- type:

  Plot sample areas, common spectra or aligned compound profiles.

- compounds:

  Compound indices to display.

## Value

A `coupledNMF` object with consensus spectra `P`, integrated `areas`,
occurrence table `compounds`, monotone `alignment`, disk-backed factor
`paths` and sample `store`, objective history, likelihood and timings,
the reusable `channel_groups` map and `parallel` thread allocation.
`fitted(object, sample=...)` reconstructs one sample on its original
grid.

## Details

Let `u = w_s(t)` be a monotone alignment map and `J` its binwise
Jacobian. With supplied exposure `e`, the original-grid mean is
\$\$\mu\_{sij}=e\_{si}\\J\_{si}\sum_g a\_{sg}h\_{sg}(u\_{si})p\_{gj}
+(B_s A_s)\_{ij}\\+o\_{sij}.\$\$ Spectra satisfy \\\sum_j p\_{gj}=1\\;
profiles satisfy \\\sum_i\Delta u\_{si}h\_{sg}(u\_{si})=1\\ and have
concave logarithms in actual aligned time. Thus areas integrate the
fitted rate over original time. `exposure=1` models the supplied
intensities. For integrated counts, known integration durations can
instead be supplied as exposure. Scan spacing, per-mass dwell time and
vendor abundance units are distinct quantities. No exposure is inferred
from the number of nonzero mass traces.

Reciprocal spectral matches supply leading-edge anchors: the 1%-height
start plus half the local narrow-peak width. A robust low-dimensional
drift spline includes log area ratio as a nuisance covariate; only drift
enters the warp. Positive-slope projection gives a monotone map.
Grouping pools compatible spectra algebraically, with at most one
occurrence per sample and an area-ordered profile-overlap check. It
performs no model refits. Subsequent local profile, pooled spectral and
background updates are backtracked against the full original-channel
likelihood. A grouped proposal is retained only when the complete
penalized objective improves.

The grouping criterion counts shared spectral parameters and detected
sample-specific amplitudes; hard shape constraints control profile
geometry. It is a rank-selection criterion with unpenalized profile and
baseline complexity. Quasipoisson uses the scaled likelihood kernel
described in
[`splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md).
Noise scales are held fixed when comparing groupings. Undetected
sample-compound combinations have `NA` area and are not filled in.

Each sample retains its own sparse PSF design and original scan grid.
The shared cache contains scans-by-aggregate matrices, while
original-channel observations remain on disk. Shared spectral and
background updates read at most 32 original channels per sample block.
Poisson/quasipoisson profile proposals use sums of channels and are
accepted against the original-channel likelihood in the affected time
window. Other noise families use their original-channel profile
objective. Parallel sample updates hold spectra fixed, and spectral
contributions are reduced in acquisition order. Only the assigned sample
and compact factors cross worker sockets. The cache and live worker
connections are released before returning the fitted object.

## See also

[`splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md),
[`estimate_psf()`](https://tomwenseleers.github.io/L0nmf/reference/estimate_psf.md),
[`simulate_GCMS_samples()`](https://tomwenseleers.github.io/L0nmf/reference/simulate_GCMS_samples.md)

## Examples

``` r
set.seed(19)
n <- 100L; m <- 6L
X <- convmatrix(dnorm(-6:6, sd=1.5), n)
H <- X[, c(25, 65)]; H <- H %*% Matrix::Diagonal(x=1/colSums(H))
S <- matrix(rexp(2*m, 1/2000), 2, m)
samples <- lapply(1:3, function(s) list(
  y=matrix(rpois(n*m, as.matrix(H %*% (S*s))+1), n, m),
  time=seq_len(n), X=X, background=matrix(1,n,1), small_width=9))
fit <- coupledNMF(samples, ncores=1L,
  control=list(cycles=2L, initial_cycles=1L, rounds=1L,
    initial=list(profile_cycles=0L, original_cleanup=FALSE)))
fit$areas
GIC(fit)
```
