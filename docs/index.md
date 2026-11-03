# L0nmf

SPARSE NONNEGATIVE MATRIX FACTORIZATION · LOG-CONCAVE PROFILES

Sparse factors.\
Interpretable signals.
======================

Sparse nonnegative matrix factorization, log-concave profiles and
flexible noise models.

[Get started
↗](https://tomwenseleers.github.io/L0nmf/articles/log-concave-nmf.md)[Explore
the methods](https://tomwenseleers.github.io/L0nmf/reference/index.md)

R · Local likelihood updates · Full spectra · Constant & varying PSFs

![L0nmf hex logo showing a sparse matrix factorized into profiles and
spectral loadings](logo-dark.svg)

[**Sparse rank**Recover a parsimonious set of nonnegative
components](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
[**Peak shapes**Log-concave profiles and constant or time-varying
PSFs](https://tomwenseleers.github.io/L0nmf/reference/estimate_psf.md)
[**7 families**Gaussian, Poisson, quasi-Poisson, NB2, NB1 and two Gamma
models](https://tomwenseleers.github.io/L0nmf/articles/log-concave-nmf.md)

L0nmf is part of a focused three-package toolkit:
[L0glm](https://tomwenseleers.github.io/L0glm/) supplies fast sparse GLM
selection, [L0nmf](https://tomwenseleers.github.io/L0nmf/) supplies
sparse factorization and PSF estimation, and
[GCMSprime](https://tomwenseleers.github.io/GCMSprime/) turns those
methods into a chromatographic workflow.

## Installation

``` r

install.packages("remotes")
remotes::install_github("tomwenseleers/L0nmf/L0nmf")
```

## Recover profiles and spectra

``` r

library(L0glm)
library(L0nmf)
data(GCMS_data1)
X <- convmatrix(convkernel("gaussian", 65L, list(sigma=8)), nrow(GCMS_data1))
fit <- splogconNMF(X, GCMS_data1, background=matrix(1,nrow(X),1),
                  family=quasipoisson("identity"), precluster=16L,
                  correlation="none")
fit
GIC(fit)
matplot(as.matrix(fit$H), type="l", lty=1, xlab="Scan", ylab="Elution profile")
```

The factors retain every original outcome channel. Profiles have unit
area, nonnegative spectra carry their intensities, and a smooth
unpenalized background captures baseline variation. Selection and fusion
use a parsimonious rank criterion with local likelihood updates.

![Recovered log-concave profiles and
spectra](reference/figures/gcms-data1-profiles-spectra.png)

Recovered log-concave profiles and spectra

The [worked
example](https://tomwenseleers.github.io/L0nmf/articles/log-concave-nmf.md)
shows the model, likelihood and local updates in detail.

The final original-channel pass checks component deletions, local
fusions and missed residual peaks, accepting each change only when the
complete likelihood-and-compound criterion improves.

## Methods and applications

- [Sparse log-concave factorization and its mathematical
  model](https://tomwenseleers.github.io/L0nmf/articles/log-concave-nmf.md)
- [Constant and time-varying
  PSFs](https://tomwenseleers.github.io/L0nmf/reference/estimate_psf.md)
- [Blind
  deconvolution](https://tomwenseleers.github.io/L0nmf/reference/L0glm_blind_deco.md)
- [Simulated
  chromatograms](https://tomwenseleers.github.io/L0nmf/reference/simulate_GCMS.md)
- [GC/MS input, output and reports in
  GCMSprime](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html)

Gaussian, Poisson, quasi-Poisson, negative-binomial and Gamma
observation models are supported for identity-link nonnegative signals.
The help pages describe aggregation, temporal dispersion, censored
observations and the full objective.

[`L0BestRankGNMF()`](https://tomwenseleers.github.io/L0nmf/reference/L0BestRankGNMF.md)
provides an experimental general factorization method. Its algorithm is
being developed separately from the validated local log-concave
workflow.

## Share spectra across samples

[`coupledNMF()`](https://tomwenseleers.github.io/L0nmf/reference/coupledNMF.md)
aligns related samples and estimates common spectra with sample-specific
log-concave profiles, integrated areas and smooth backgrounds. It
supports all seven signal families, irregular scan grids and disk-backed
observations. In a three-sample simulation with all 258 channels, it
recovers all 78 compound occurrences in 29.69 seconds. With compounds
absent from half the samples, all 37 present occurrences are recovered
in 21.92 seconds.

``` r

fit <- coupledNMF(samples, family=quasipoisson("identity"), precluster=32L)
fit$areas
coef(fit, normalize="l1")
```

See [Coupled multisample
factorization](https://tomwenseleers.github.io/L0nmf/articles/coupled-nmf.md)
for simulations, plots and mathematical details, and
[GCMSprime](https://tomwenseleers.github.io/GCMSprime/articles/multisample.html)
for the directory-based chromatogram workflow.

## Citation

An article describing the methods implemented in L0nmf is in
preparation. If you use the package before the article is formally
published, we would appreciate a citation to the package:

Wenseleers T, Vanderaa C (2026). *L0nmf: Scalable Sparse Nonnegative
Matrix Factorization with Automatic Rank Selection and Flexible Noise
Models*. R package version 0.1.0,
<https://tomwenseleers.github.io/L0nmf/>.

The [Authors and Citation
page](https://tomwenseleers.github.io/L0nmf/authors.html#citation) also
provides the citation in BibTeX format. Please also cite the
methodological sources linked in the help pages where relevant.

## Acknowledgements

OpenAI Codex, using GPT-6 Astra and GPT-6.1 Sol, was used to optimize
and refactor code and help draft the documentation. The key ideas and
methodological decisions are the authors’ own, and the authors take full
responsibility for the work presented here.
