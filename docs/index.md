# logconNMF

SPARSE NONNEGATIVE MATRIX FACTORIZATION · LOG-CONCAVE PROFILES

Log-concave NMF.\
Automatic rank selection.
=========================

Separate overlapping signals into log-concave profiles and nonnegative
spectra. Use the same interface for one sample or a coupled fit across
an experiment.

[Get started
↗](https://tomwenseleers.github.io/logconNMF/articles/log-concave-nmf.md)[Explore
the
methods](https://tomwenseleers.github.io/logconNMF/reference/index.md)

R · Local likelihood updates · Full spectra · Constant & varying PSFs

![logconNMF hex logo with log-concave profiles and corresponding
spectra](logo-dark.svg)

[**Automatic rank**Estimate how many components the signal
supports](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
[**Peak shape constraints**Log-concave profiles and constant or
time-varying
PSFs](https://tomwenseleers.github.io/logconNMF/reference/estimate_psf.md)
[**7 noise families**Gaussian, Poisson, quasi-Poisson, NB2, NB1 and two
Gamma
models](https://tomwenseleers.github.io/logconNMF/articles/log-concave-nmf.md)

**logconNMF estimates the number of components from the data.** It
separates overlapping signals into nonnegative channel loadings and
log-concave profiles. Fit one matrix or share spectra across multiple
samples, with a separate peak shape and area for each occurrence. Seven
noise models cover count and continuous signals, including overdispersed
counts.

Gas chromatography–mass spectrometry (GC–MS) is a motivating
application. Several compounds can elute at the same time, so each scan
contains a mixture of their mass spectra. The fit separates the elution
profiles and spectra and estimates how many distinct components the
scans support. The same model applies to other ordered multichannel
signals with nonnegative, single-peaked components.

## Automatic rank selection

For a scan-by-channel matrix `Y`, the fitted signal is
`H %*% S + B %*% A`. Each column of `H` is a profile across scans and
its corresponding row of `S` contains the channel loadings. The optional
term `B %*% A` represents a smooth baseline. A known additive signal can
be supplied as an offset.

Both modes approximately minimize **−2 × log-likelihood + κ × λ × d**.
The weighted log-likelihood sums over the original channels and, for a
coupled fit, over samples. Quasi-Poisson divides the Poisson log score
by the estimated or supplied dispersion. The default `lambda = "gic"`
sets the penalty to `log(p) × log(log(n))` per selection degree of
freedom. For one sample, `p` is the PSF dictionary width and `n` is the
number of scans. For a coupled fit, `p` and `n` both equal the total
scan count. The optional temporal adjustment κ is 1 under independence.

| Fit | Selection count `d` for the components |
|:---|:---|
| One sample | `m × K` for `K` components in `m` channels |
| Coupled samples | `(m − 1) × G + L` for `G` shared spectra and `L` detected sample–compound occurrences |

Each shared spectrum has unit sum, giving `m − 1` spectral coordinates,
and each occurrence contributes one amplitude. In one sample, each
spectrum has one occurrence, so `G = L = K` and the counts agree.
Profile shapes and background coefficients carry no separate rank
charge. The reported information criteria also count estimated noise
parameters.

The search refines candidate profiles and compares alternative ranks
using the criterion. For example, overlapping components can fuse when
their combined profile is log-concave and the criterion improves. The
resulting `fit$rank` is the estimated number of components or shared
spectra. See [the mathematical
objective](https://tomwenseleers.github.io/logconNMF/articles/algorithm.html#common-rank-objective).

## Installation

``` r

install.packages("remotes")
remotes::install_github("tomwenseleers/logconNMF", subdir="logconNMF")
```

## Recover overlapping profiles and spectra

This short experimental window contains six co-eluting hydrocarbons. The
initial dictionary `X` places a narrow Gaussian point spread function
(PSF) at each possible peak location. The factorization can move away
from that shape and gives each retained component its own log-concave
profile.

``` r

library(logconNMF)
data(GCMS_data1)
Y <- GCMS_data1
g <- convkernel("gaussian", 65L, list(sigma = 8))
X <- convmatrix(g, nrow(Y))
fit <- logconNMF(Y, X = X, background = matrix(1, nrow(Y), 1L),
  family = quasipoisson("identity"), precluster = 16L)

fit$rank                         # estimated number of components
profiles <- fit$H                 # scans x components
loadings <- coef(fit)             # components x original channels
spectra <- coef(fit, normalize = "l1")
GIC(fit)
```

The columns of `profiles` sum to one, so each row sum in `loadings`
gives the component’s total fitted signal across scans and channels.
Grouping similar channels can speed up the initial search. The reported
spectra use every original channel.

![Observed signal, fitted log-concave profiles and recovered spectra for
six hydrocarbons](reference/figures/gcms-data1-decomposition.png)

Observed signal, fitted log-concave profiles and recovered spectra for
six hydrocarbons

Each spectrum is displayed on a square-root relative-intensity scale.

## Choose a workflow

| Data and objective | Method |
|:---|:---|
| Estimate a common response shape from a standard or overlapping signals | [Estimate PSFs](https://tomwenseleers.github.io/logconNMF/articles/estimate-psf.md) with [`estimate_psf()`](https://tomwenseleers.github.io/logconNMF/reference/estimate_psf.md). The returned sparse design can use a constant or time-varying PSF. |
| Resolve components with individual shapes in one matrix | [Sparse log-concave NMF](https://tomwenseleers.github.io/logconNMF/articles/log-concave-nmf.md) with [`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md). The tutorial covers rank selection, backgrounds and seven noise models. |
| Share spectra across related samples | [Coupled log-concave NMF](https://tomwenseleers.github.io/logconNMF/articles/coupled-nmf.md) with [`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md). Each sample retains its own scan grid, peak shapes and areas. Components may be absent from individual samples. |

For a known and common PSF,
[L0glm](https://tomwenseleers.github.io/L0glm/articles/deconvolution.html)
fits the sparse convolution model directly.
[GCMSprime](https://tomwenseleers.github.io/GCMSprime/) uses logconNMF
for chromatogram files, retention-index calibration and library search.
The [algorithm and
theory](https://tomwenseleers.github.io/logconNMF/articles/algorithm.md)
give the objectives, shape constraints and mathematical pseudocode.

## Share spectra across samples

In the [16-sample
benchmark](https://tomwenseleers.github.io/logconNMF/articles/coupled-nmf.html#sixteen-samples-with-parallel-fitting),
the coupled model recovers all 48 shared compounds and 380 occurrences,
with no false calls. The CDF workflow took **63.91 s** with 16 workers.
Each sample has 1,200 scans and 258 ion channels, with retention-time
drift and concentration-dependent peak shapes. The tutorial shows the
recovered areas and timings by stage.

## Citation

An article describing the methods implemented in logconNMF is in
preparation. If you use the package before the article is formally
published, we would appreciate a citation to the package:

Wenseleers T, Vanderaa C (2026). *logconNMF: Log-Concave Nonnegative
Matrix Factorization with Automatic Rank Selection*. R package version
0.1.0, <https://tomwenseleers.github.io/logconNMF/>.

The [Authors and Citation
page](https://tomwenseleers.github.io/logconNMF/authors.html#citation)
also provides the citation in BibTeX format. Please also cite the
methodological sources linked in the help pages where relevant.

## Acknowledgements

OpenAI Codex, using GPT-6 Astra and GPT-6.1 Sol, was used to optimize
and refactor code and help draft the documentation. The key ideas and
methodological decisions are the authors’ own, and the authors take full
responsibility for the work presented here.
