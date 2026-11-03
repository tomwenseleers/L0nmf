# Sparse log-concave NMF with automatic best rank tuning

A convolution model gives every impulse the same local response shape.
Real components can be broader, more asymmetric or more strongly tailed.
Sparse log-concave NMF lets each component acquire its own positive
profile and loading vector, while selecting how many components the
observations support.

[`splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
starts from sparse multichannel deconvolution on a candidate PSF design,
fuses compatible components, refines their profiles and loadings, and
checks residual evidence for additions or deletions. Local updates use
the specified observation likelihood and hard log-concavity constraints.

## Fit profiles, loadings and background

The signal below contains six co-eluting components. Rows index ordered
scans; columns index mass channels. The formulation also applies to
other ordered multichannel measurements.

``` r

library(L0glm)
library(L0nmf)
data(GCMS_data1)
Y <- GCMS_data1
X <- convmatrix(convkernel("gaussian", 65L, list(sigma=8)), nrow(Y))
B <- matrix(1, nrow(Y), 1L)
fit <- splogconNMF(X, Y, background=B,
  family=quasipoisson("identity"), precluster=16L, correlation="auto")
fit$rank
dim(fit$H)
dim(coef(fit))
GIC(fit)
```

The mean is `H %*% S + B %*% alpha + offset`. Each sparse column of `H`
has contiguous support, unit sum and a concave log profile. Rows of `S`
contain loadings on every original channel; their sums give the
component’s total discrete integrated intensity. `B` is a nonnegative
background basis. For unknown kernels, [estimate the
PSF](https://tomwenseleers.github.io/L0nmf/articles/estimate-psf.md)
first and pass its `X`.

## Inspect the recovered components

``` r

matplot(as.matrix(fit$H), type="l", lty=1, xlab="Scan", ylab="Profile")
matplot(t(coef(fit)), type="h", lty=1, xlab="Channel", ylab="Loading")
plotmat(t(fitted(fit)), xlab="Scan", ylab="Channel")
patterns <- coef(fit) / rowSums(coef(fit))
fitted_signal <- fitted(fit)
fit$seconds
```

![Recovered component profiles and full-channel loading
patterns](../reference/figures/gcms-data1-profiles-spectra.png)

Recovered component profiles and full-channel loading patterns

All 671 channels are retained, and factors are returned in increasing
profile-apex order. The [GCMSprime
tutorial](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html)
develops the chromatographic interpretation, file input and reports.

## A broad asymmetric component

The beer-volatile window combines narrow overlapping peaks with a much
broader acid contribution. A low-dimensional spline accounts for the
smooth background:

``` r

data(GCMS_data2)
Y2 <- GCMS_data2
t2 <- as.numeric(rownames(Y2))
X2 <- convmatrix(convkernel("gaussian", 65L, list(sigma=6)), nrow(Y2))
B2 <- splines::bs(t2, df=4L, intercept=TRUE)
fit2 <- splogconNMF(X2, Y2, background=B2,
  family=quasipoisson("identity"), precluster=32L, correlation="auto")
fit2$rank
fit2$contained_fusion$broad
```

![Log-concave profiles and full-channel spectra, including the broad
component](../reference/figures/gcms-data2-profiles-spectra.png)

Log-concave profiles and full-channel spectra, including the broad
component

The corresponding GCMSprime CDF workflows use the same full-channel
backend and the default temporal information adjustment:

| Window         | Channels | Components | Quasi-GIC | Adjusted weighted R² |
|:---------------|---------:|-----------:|----------:|---------------------:|
| Hydrocarbons   |      671 |          6 |   246,473 |              0.98999 |
| Beer volatiles |      518 |         38 | 1,467,419 |              0.99498 |

The accompanying [CDF
workflow](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html)
reports measured timings and complete-criterion tables for both windows,
with downloadable traces, heatmaps, complete profile/spectrum panels and
MSP libraries.

## Automatic best rank tuning

The initialization keeps a permissive pool of components while native
constrained likelihood updates refine their profiles and spectra. Final
rank selection applies the requested temporal information adjustment.
Small fragments contained within a broad profile are compared with an
equally refined unmerged local model on all original channels.

The default `lambda="gic"` compares complete fits by likelihood and
component count. After dividing the likelihood by the original channel
count, each retained component incurs one rank charge. On the unscaled
likelihood the charge is `m * lambda_IC`, with an optional temporal
information adjustment. Profile complexity and background coefficients
have no extra rank penalty.

``` r

fit$objective
fit$history[, c("stage", "delta", "accepted", "compounds")]
fit$cleanup
fit$initialization
fit$contained_fusion
c(AIC=AIC(fit), BIC=BIC(fit), GIC=GIC(fit))
```

Reported criteria count `m * rank` plus estimated shared noise
parameters. Other presets are `"aic"`, `"bic"`, `"ebic"` and `"mbic"`. A
nonnegative numeric `lambda` supplies the total component penalty before
temporal adjustment. The [rank
objective](https://tomwenseleers.github.io/L0nmf/articles/algorithm.html#rank-selection)
gives the full convention.

## Choose a noise model

All seven models use identity-link nonnegative means and the original
observation scale. Gaussian measurements can have negative errors around
that mean.

| Model | R specification | Variance at mean `mu` |
|:---|:---|:---|
| Gaussian | `gaussian("identity")` | `dispersion` |
| Poisson | `poisson("identity")` | `mu` |
| Quasipoisson | `quasipoisson("identity")` | `dispersion * mu` |
| NB2 | `mgcv::nb(link="identity")` | `mu + mu^2/theta` |
| NB1 | `nb1(dispersion=2, link="identity")` | `dispersion * mu` |
| Censored Gamma | `Gamma_censored(dispersion=.2, link="identity")` | `dispersion * mu^2` before censoring |
| Common-scale Gamma | `Gamma_scale(scale=2, link="identity")` | `scale * mu` before censoring |

Positive `dispersion` or `theta` arguments fix the parameter. Omitted
shared parameters are estimated from temporal differences. A negative
NB2 `theta` initializes at its absolute value and enables
likelihood-improving updates. Quasipoisson uses a scaled Poisson kernel
and corresponding quasi-criteria.

Zero Gamma intensities represent left censoring at a detection limit,
set by `control=list(censor=0.5)`; their contribution is the probability
below that limit. See the [noise
models](https://tomwenseleers.github.io/L0nmf/articles/algorithm.html#noise-models)
for scores and parameterizations.

## Faster initialization and smooth backgrounds

`precluster=32L` initializes in at most 32 size-adjusted outcome groups,
then returns to every original channel. `precluster=FALSE` uses the full
matrix directly. `control$linkage` selects Ward Poisson or Bregman
grouping. The original channel count sets the rank charge in both
workflows; final scoring uses the original observations and specified
noise model.

For a longer signal, replace the constant background with a
low-dimensional nonnegative spline basis:

``` r

B <- splines::bs(seq_len(nrow(Y)), df=4L, intercept=TRUE)
```

Background coefficients are fitted conditionally on the component signal
and remain unpenalized in rank selection. Known additive signals use
`offset`; weights can be shared across channels or channel-specific.

The default `correlation="auto"` estimates a lag-one working information
adjustment. `"none"` uses independent-observation rank charges; a
supplied nonnegative correlation fixes the adjustment. Marginal noise
variance stays in the likelihood. The [mathematical
article](https://tomwenseleers.github.io/L0nmf/articles/algorithm.md)
develops local solvers, hard shape constraints, aggregation and complete
pseudocode. Continue with [coupled multisample
NMF](https://tomwenseleers.github.io/L0nmf/articles/coupled-nmf.md) for
shared loadings and sample-specific profiles.
