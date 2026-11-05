# L0nmf

<div class="l0-hero">
<div class="l0-hero-copy">
<p class="l0-eyebrow">SPARSE NONNEGATIVE MATRIX FACTORIZATION · LOG-CONCAVE PROFILES</p>
<h1>Sparse factors.<br><span>Interpretable signals.</span></h1>
<p class="l0-lead">Separate overlapping signals into log-concave profiles and nonnegative spectra. Estimate the number of components in one sample or share spectra across a whole experiment.</p>
<div class="l0-actions"><a class="btn btn-primary" href="articles/log-concave-nmf.html">Get started <span aria-hidden="true">↗</span></a><a class="btn btn-outline-light" href="reference/index.html">Explore the methods</a></div>
<p class="l0-hero-note">R · Local likelihood updates · Full spectra · Constant &amp; varying PSFs</p>
</div>
<img class="l0-hero-logo" src="logo-dark.svg" alt="L0nmf hex logo showing a sparse matrix factorized into profiles and spectral loadings" width="240" height="280">
</div>

<div class="l0-stat-grid">
<a href="reference/splogconNMF.html"><strong>Sparse rank</strong><span>Recover a parsimonious set of nonnegative components</span></a>
<a href="reference/estimate_psf.html"><strong>Peak shape constraints</strong><span>Log-concave profiles and constant or time-varying PSFs</span></a>
<a href="articles/log-concave-nmf.html"><strong>7 noise families</strong><span>Gaussian, Poisson, quasi-Poisson, NB2, NB1 and two Gamma models</span></a>
</div>

In a GC–MS chromatogram, several compounds can elute at the same time.
Each scan records the sum of their mass spectra, with contributions that
change as the compounds pass through the instrument. L0nmf separates this
mixture into an elution profile and a spectrum for each component, while
estimating how many components the signal supports. The same model applies
to ordered multichannel signals with nonnegative, single-peaked components.

## What the model estimates

For a scan-by-channel matrix `Y`, the fitted signal is `H %*% S + B %*% A`.
Column `k` of `H` gives the shape of component `k` across scans. Row `k`
of `S` gives its contributions to the original channels. A low-dimensional
basis `B` can describe a smooth baseline, with coefficients in `A`. A known
additive signal can also be supplied as an offset.

On its positive interval, a **log-concave profile** has non-increasing slopes on
the log scale. This admits symmetric peaks, asymmetric tails and broad
triangular peaks. It gives the fit a way to separate overlapping components
while allowing their shapes to differ. A rank penalty charges for each
retained component, and local likelihood updates refine the profiles and
spectra under the chosen noise model.

## Installation

```r
install.packages("remotes")
remotes::install_github("tomwenseleers/L0nmf", subdir="L0nmf")
```

## Recover overlapping profiles and spectra

This short experimental window contains six co-eluting hydrocarbons. The
initial dictionary `X` places a narrow Gaussian point spread function (PSF)
at each possible peak location. The factorization can move away from that
shape and gives each retained component its own log-concave profile.

```r
library(L0nmf)
data(GCMS_data1)
Y <- GCMS_data1
g <- convkernel("gaussian", 65L, list(sigma = 8))
X <- convmatrix(g, nrow(Y))
fit <- splogconNMF(X, Y, background = matrix(1, nrow(Y), 1L),
  family = quasipoisson("identity"), precluster = 16L)

fit$rank                         # estimated number of components
profiles <- fit$H                 # scans x components
loadings <- coef(fit)             # components x original channels
spectra <- loadings / rowSums(loadings)
GIC(fit)
```

The columns of `profiles` sum to one, so each row sum in `loadings`
gives the component's total fitted signal across scans and channels.
Grouping similar channels can speed up the initial search. The reported
spectra use every original channel.

![Observed signal, fitted log-concave profiles and recovered spectra for six hydrocarbons](reference/figures/gcms-data1-decomposition.png)

Each spectrum is displayed on a square-root relative-intensity scale.

## Choose a workflow

| Data and objective | Method |
|:--|:--|
| Estimate a common response shape from a standard or overlapping signals | [Estimate PSFs](articles/estimate-psf.html) with `estimate_psf()`. The returned sparse design can use a constant or time-varying PSF. |
| Resolve components with individual shapes in one matrix | [Sparse log-concave NMF](articles/log-concave-nmf.html) with `splogconNMF()`. The tutorial covers rank selection, backgrounds and seven noise models. |
| Share spectra across related samples | [Coupled log-concave NMF](articles/coupled-nmf.html) with `coupledNMF()`. Each sample retains its own scan grid, peak shapes and areas. Components may be absent from individual samples. |

For a known and common PSF, [L0glm](https://tomwenseleers.github.io/L0glm/articles/deconvolution.html)
fits the sparse convolution model directly. [GCMSprime](https://tomwenseleers.github.io/GCMSprime/)
uses L0nmf for chromatogram files, retention-index calibration and library
search. The [algorithm and theory](articles/algorithm.html) give the
objectives, shape constraints and mathematical pseudocode.

## Share spectra across samples

In the [16-sample benchmark](articles/coupled-nmf.html#sixteen-samples-with-parallel-fitting),
the coupled model recovers all 48 shared compounds and 380 occurrences,
with no false calls. The CDF workflow took **63.91 s** with 16 workers. Each sample has 1,200 scans and 258 ion channels,
with retention-time drift and concentration-dependent peak shapes.
The tutorial shows the recovered areas and timings by stage.

## Citation

An article describing the methods implemented in L0nmf is in preparation.
If you use the package before the article is formally published, we would
appreciate a citation to the package:

Wenseleers T, Vanderaa C (2026).
*L0nmf: Scalable Sparse Nonnegative Matrix Factorization with Automatic Rank Selection and Flexible Noise Models*.
R package version 0.1.0,
<https://tomwenseleers.github.io/L0nmf/>.

The [Authors and Citation page](authors.html#citation) also provides the
citation in BibTeX format. Please also cite the methodological sources linked
in the help pages where relevant.

## Acknowledgements

OpenAI Codex, using GPT-6 Astra and GPT-6.1 Sol, was used to optimize and
refactor code and help draft the documentation. The key ideas and methodological
decisions are the authors' own, and the authors take full responsibility for the
work presented here.
