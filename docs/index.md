# L0nmf

SPARSE NONNEGATIVE MATRIX FACTORIZATION · LOG-CONCAVE PROFILES

Sparse factors.\
Interpretable signals.
======================

Scalable sparse nonnegative matrix factorization with log-concave
profile shape constraints, flexible noise models and automatic best rank
tuning for single or multiple samples.

[Get started
↗](https://tomwenseleers.github.io/L0nmf/articles/estimate-psf.md)[Explore
the methods](https://tomwenseleers.github.io/L0nmf/reference/index.md)

R · Local likelihood updates · Full spectra · Constant & varying PSFs

![L0nmf hex logo showing a sparse matrix factorized into profiles and
spectral loadings](logo-dark.svg)

[**Sparse rank**Recover a parsimonious set of nonnegative
components](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
[**Peak shape constraints**Log-concave profiles and constant or
time-varying
PSFs](https://tomwenseleers.github.io/L0nmf/reference/estimate_psf.md)
[**7 noise families**Gaussian, Poisson, quasi-Poisson, NB2, NB1 and two
Gamma
models](https://tomwenseleers.github.io/L0nmf/articles/log-concave-nmf.md)

Scalable sparse nonnegative matrix factorization with log-concave
profile shape constraints, flexible noise models and automatic best rank
tuning for single or multiple samples.

Start with a calibrated or estimated point spread function, then let
each component acquire its own log-concave profile. Coupled models share
spectra across samples while retaining sample-specific shapes, areas and
observation grids.

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
                  correlation="auto")
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
shows fitting and diagnostics. The [algorithm and
theory](https://tomwenseleers.github.io/L0nmf/articles/algorithm.md)
give the likelihood, constraints and mathematical pseudocode.

The final original-channel pass checks component deletions, local
fusions and missed residual peaks, accepting each change only when the
complete likelihood-and-compound criterion improves. Native log-concave
profile updates and joint spectral fits start from the current factors.
Small fragments contained in broad profiles are compared with equally
refined unmerged local models.

## Learn the methods

1.  [Estimate
    PSFs](https://tomwenseleers.github.io/L0nmf/articles/estimate-psf.md):
    known convolution, calibration and blind kernel estimation.
2.  [Sparse log-concave NMF with automatic best rank
    tuning](https://tomwenseleers.github.io/L0nmf/articles/log-concave-nmf.md).
3.  [Coupled multisample log-concave NMF with automatic best rank
    tuning](https://tomwenseleers.github.io/L0nmf/articles/coupled-nmf.md).

The [mathematical
article](https://tomwenseleers.github.io/L0nmf/articles/algorithm.md)
develops the objectives and algorithms.
[L0glm](https://tomwenseleers.github.io/L0glm/) supplies sparse GLM
selection; [GCMSprime](https://tomwenseleers.github.io/GCMSprime/)
applies these models to chromatogram files, quantification and spectral
identification.

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

See [Coupled multisample log-concave
NMF](https://tomwenseleers.github.io/L0nmf/articles/coupled-nmf.md) for
simulations, plots and mathematical details, and
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
