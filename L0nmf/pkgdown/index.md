# L0nmf

<div class="l0-hero">
<div class="l0-hero-copy">
<p class="l0-eyebrow">SPARSE NONNEGATIVE MATRIX FACTORIZATION · LOG-CONCAVE PROFILES</p>
<h1>Sparse factors.<br><span>Interpretable signals.</span></h1>
<p class="l0-lead">Sparse nonnegative matrix factorization, log-concave profiles and flexible noise models.</p>
<div class="l0-actions"><a class="btn btn-primary" href="articles/log-concave-nmf.html">Get started <span aria-hidden="true">↗</span></a><a class="btn btn-outline-light" href="reference/index.html">Explore the methods</a></div>
<p class="l0-hero-note">R · Local likelihood updates · Full spectra · Constant &amp; varying PSFs</p>
</div>
<img class="l0-hero-logo" src="logo-dark.svg" alt="L0nmf hex logo showing a sparse matrix factorized into profiles and spectral loadings" width="240" height="280">
</div>

<div class="l0-stat-grid">
<a href="reference/splogconNMF.html"><strong>Sparse rank</strong><span>Recover a parsimonious set of nonnegative components</span></a>
<a href="reference/estimate_psf.html"><strong>Peak shapes</strong><span>Log-concave profiles and constant or time-varying PSFs</span></a>
<a href="articles/log-concave-nmf.html"><strong>5 families</strong><span>Gaussian, Poisson, quasi-Poisson, negative binomial and Gamma</span></a>
</div>

L0nmf is part of a focused three-package toolkit: [L0glm](https://tomwenseleers.github.io/L0glm/) supplies fast sparse GLM selection, [L0nmf](https://tomwenseleers.github.io/L0nmf/) supplies sparse factorization and PSF estimation, and [GCMSprime](https://tomwenseleers.github.io/GCMSprime/) turns those methods into a chromatographic workflow.

## Installation

```r
install.packages("remotes")
remotes::install_github("tomwenseleers/L0nmf/L0nmf")
```

## Recover profiles and spectra

```r
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

The factors retain every original outcome channel. Profiles have unit area, nonnegative spectra carry their intensities, and a smooth unpenalized background captures baseline variation. Selection and fusion use a parsimonious rank criterion with local likelihood updates.


![Recovered log-concave profiles and spectra](reference/figures/gcms-data1-profiles-spectra.png)

The [worked example](articles/log-concave-nmf.html) shows the model, likelihood and local updates in detail.

## Methods and applications

* [Sparse log-concave factorization and its mathematical model](articles/log-concave-nmf.html)
* [Constant and time-varying PSFs](reference/estimate_psf.html)
* [Blind deconvolution](reference/L0glm_blind_deco.html)
* [Simulated chromatograms](reference/simulate_GCMS.html)
* [GC/MS input, output and reports in GCMSprime](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html)

Gaussian, Poisson, quasi-Poisson, negative-binomial and Gamma observation models are supported for identity-link nonnegative signals. The help pages describe aggregation, temporal dispersion, censored observations and the full objective.

`L0BestRankGNMF()` provides an experimental general factorization method. Its algorithm is being developed separately from the validated local log-concave workflow.

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
