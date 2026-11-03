# Algorithms and theory for PSFs and sparse log-concave factorization

L0nmf connects convolution and blind system identification to
nonnegative factorization with flexible noise models. The ordered factor
has a positive, log-concave profile; the other factor gives its
multichannel loadings. Rank selection, backgrounds and noise calibration
are part of the model.

This article follows the three tutorials: [PSF
estimation](https://tomwenseleers.github.io/L0nmf/articles/estimate-psf.md),
[single-sample
factorization](https://tomwenseleers.github.io/L0nmf/articles/log-concave-nmf.md)
and [coupled multisample
factorization](https://tomwenseleers.github.io/L0nmf/articles/coupled-nmf.md).
It specifies the actual objectives, coordinate conventions,
initialization and acceptance rules.

## Convolution and PSF estimation

Let Y\in\mathbb R^{n\times m} contain ordered multichannel observations,
C\in\mathbb R\_+^{p\times m} a spike-amplitude matrix and g a kernel.
With a nonnegative background basis B, coefficients A\geq0 and known
offset O, the fixed-kernel model is

M=X(g)C+BA+O,\qquad C,A\geq0.

For a constant kernel, X(g)\_{ir}=g\_{i-r}, with the specified centering
and boundary convention. Sparse regression selects individual entries or
whole rows of C. Its target combines the observation log likelihood with
an entrywise or row-support penalty. L0glm supplies the adaptive-ridge
IRLS, normalization, compulsory-background treatment and relaxed refit;
its [algorithm
article](https://tomwenseleers.github.io/L0glm/articles/algorithm.html)
gives those updates in detail. The sparse pulse-train formulation builds
on [de Rooi & Eilers (2011)](https://doi.org/10.1016/j.aca.2011.05.030).

### Calibration with isolated peaks

For a local isolated standard peak, write

\mu\_{ij}=b\_{ij}+o\_{ij}+h_i s_j, \qquad h_i\>0,\quad s_j\geq0,\quad
D^2\log h\leq0.

A low-quantile spline estimates background. Signal-to-noise, separation
and profile geometry identify informative peak windows. A native
constrained log-profile update fits each window on the original
observation scale. Normalize its profile, align by centroid, and
resample onto the sample’s median scan interval. A high quantile of
fitted profile support widths supplies the common support when
`control$L` is omitted. `control$L` fixes an odd length.

**Algorithm P1 — standard PSF calibration**

1.  Resolve the shared noise parameters, weights and additive offsets.
    Optionally aggregate channels; carry the corresponding sum-model
    parameters.
2.  Estimate a lower-quantile background and detect informative isolated
    peaks.
3.  For each retained local window, initialize s_j from positive
    integrated signal and h from its positive total trace projected to
    log concavity.
4.  Update \eta=\log h under D^2\eta\leq0, using the likelihood score
    and positive working curvature. Backtrack against the original
    likelihood.
5.  Set h\leftarrow h/\sum_i h_i, align its centroid, and place it on
    the common lag grid. Infer support from fitted peak mass or use the
    supplied length.
6.  Interpolate the calibrated log kernels at the requested PSF anchor
    positions, normalize and construct the nonnegative sparse
    convolution design.

### Alternating blind estimation

For a fixed channel spike train c_j, convolution is linear in the
kernel:

X(g)c_j=T(c_j)g,

where T(c_j) contains shifted copies of the spike train. This gives a
second nonnegative GLM. Its ridge center is the current kernel. The
implementation normalizes each informative channel’s spike train for
conditioning and rescales the prior center consistently. Positive starts
retain kernel tails. Its default kernel update keeps alternating lag
positions as compulsory terms; the remaining positions use adaptive
selection. A log-concave smoother and centroid alignment produce the
kernel proposal, which is normalized to unit sum.

Completed proposals are compared after re-estimating their spike trains
and backgrounds. A full-tail proposal and a pedestal-removed proposal
can both be considered. The returned checkpoint is the synchronized
spike/background fit with the lowest requested criterion among these
completed candidates.

**Algorithm P2 — blind constant-PSF estimation**

1.  Set g^{(0)} to the supplied/calibrated kernel, or a normalized
    Gaussian on the requested lag grid. Resolve shared noise parameters
    once in
    [`estimate_psf()`](https://tomwenseleers.github.io/L0nmf/reference/estimate_psf.md).
2.  For t=0,\ldots,T-1:
    1.  Construct X(g^{(t)}). Fit \[X(g^{(t)}),B\] to Y, with compulsory
        background and the chosen sparse/grouped GLM preset.
    2.  Store this spike/background fit and its criterion.
    3.  For every informative channel, construct T(c_j) and fit its
        kernel using background plus known signal as offset, the current
        kernel as ridge center, and the same observation weights/noise
        model.
    4.  Smooth, recenter and normalize the channel kernels; combine
        informative channel estimates. Construct the optional
        concentrated alternative.
    5.  Refit spike/background models for the competing proposals and
        select the better completed candidate.
    6.  Stop if \\g^{(t+1)}-g^{(t)}\\\_1/\\g^{(t)}\\\_1\<\epsilon.
3.  Synchronize the last kernel with a spike/background refit. Return
    the best scored checkpoint, its sparse design and its selection
    trace.

`control$maxiter` bounds this loop and `control$tol` is \epsilon.
`method="auto"` runs P1 followed by P2; `"standard"` runs P1 and
`"blind"` runs P2. A supplied kernel with `maxiter=0L` constructs its
design directly.

### Smoothly time-varying kernels

With `nwindows>1L`, run local estimation in overlapping windows centered
at \tau_1\<\cdots\<\tau_R. Use a common lag grid and interpolate
adjacent log kernels:

\eta_t=(1-a)\log g_r+a\log g\_{r+1},\quad a=u^2(3-2u),\quad
u=\frac{t-\tau_r}{\tau\_{r+1}-\tau_r}, \qquad
g_t=\frac{\exp\eta_t}{\sum\_\ell\exp\eta\_{t\ell}}.

A convex combination of concave log profiles is concave; normalization
changes only its intercept. The smoothstep weights give zero
interpolation slope at anchor joins. Each design column uses its
position’s kernel, is projected to contiguous log-concave support, and
is normalized after boundary truncation. The measured observations stay
on their scan grid; kernels are indexed in scans.

**Algorithm P3 — windowed PSF estimation**

1.  Choose evenly spaced sample anchor positions and overlapping local
    windows.
2.  Run P1/P2 for each anchor, retaining a common lag grid and its
    diagnostics.
3.  Interpolate log kernels using the weights above, with endpoint
    continuation.
4.  Place each local kernel at its candidate impulse position, truncate
    at the observation boundary, enforce log concavity and normalize the
    column.
5.  Return sparse X, anchor kernels, times, widths and convergence
    diagnostics.

These PSFs describe a common response at each position. The
factorization below subsequently estimates component-specific shapes.

## Single-sample model and rank criterion

For scans i, channels j and compounds r,

\mu\_{ij}=\sum\_{r=1}^K H\_{ir}S\_{rj}+(BA)\_{ij}+O\_{ij},\qquad
H,S,A\geq0.

The columns of H have unit sum, contiguous support and \Delta^2\log
H\_{\cdot r}\leq0. Spectra S therefore have the interpretation of
integrated channel intensities. The objective is

J=-2\sum\_{ij}w\_{ij}\ell(y\_{ij};\mu\_{ij},\nu)+\lambda_R K, \qquad
\lambda_R=m\kappa\lambda\_{IC}.

Here m is the original channel count, \nu contains pooled noise
parameters and \kappa=1 for independent observations. Numeric `lambda`
specifies the total compound penalty before the optional covariance
adjustment. Background coefficients and profile complexity have no
additional penalty. The original channel count is retained in this
penalty when the observations are aggregated, so changing the number of
outcome groups does not change the compound penalty. Dividing by the
channel count gives the equivalent criterion
-2\ell/m+\kappa\lambda\_{IC}K, with one rank penalty per compound. The
reported IC methods use the compound-count convention mK, plus any
estimated shared noise parameters.

## Generalized matrix factorization

The broad statistical setting is **generalized matrix factorization
(GMF)**: an exponential-family mean is linked to a low-rank product,
optionally with observed covariate effects and offsets. See [Kidziński
et al. (2022)](https://www.jmlr.org/papers/v23/20-1104.html) and
[Castiglione et
al. (2026)](https://doi.org/10.1093/biostatistics/kxag010). For example,
the latter mean-model structure is

g(\mu)=X\_{\rm row}B\_{\rm row}^{\top} +B\_{\rm col}X\_{\rm
col}^{\top}+UV^{\top}+O.

Our model uses the identity link, U=H, V^{\top}=S, and X\_{\rm row}=B,
B\_{\rm row}^{\top}=A for the smooth background; the column-covariate
term is omitted. It adds nonnegativity, unit-area log-concave profiles
and a compound-count penalty to select the rank. Thus **sparse
log-concave nonnegative GMF** describes the model, while the function
retains the familiar
[`splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
name.

## Deviance matrix factorization and Bregman NMF

[`splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
is a **rank-penalized, nonnegative deviance matrix factorization with
log-concave profiles**. It belongs to the deviance-based framework of
[Wang & Carvalho (2023)](https://doi.org/10.1214/23-EJS2174), using an
identity link for the additive mean \mu=HS+BA+O. The loss also connects
directly to [NMF with Bregman
divergences](https://doi.org/10.3390/s23167229).

For a regular exponential-dispersion model with fixed noise parameters,

\ell(y;\theta,\phi)=\frac{y\theta-b(\theta)}{\phi}+c(y,\phi), \qquad
\mu=b'(\theta),

let \Psi=b^\* be the convex conjugate of b. Its Bregman divergence is

D\_\Psi(y,\mu)=\Psi(y)-\Psi(\mu)-\Psi'(\mu)(y-\mu).

Since \Psi'(\mu)=\theta(\mu) and
\Psi(\mu)=\mu\theta(\mu)-b\\\theta(\mu)\\, substitution gives

\begin{aligned} D\_\Psi(y,\mu)
&=y\theta(y)-b\\\theta(y)\\-y\theta(\mu)+b\\\theta(\mu)\\,\\
2\\\ell(y;y,\phi)-\ell(y;\mu,\phi)\\ &=\frac{2}{\phi}D\_\Psi(y,\mu).
\end{aligned}

Here the likelihood on the second line is written in mean coordinates.
Boundary observations use the corresponding limits. Thus, conditional on
the noise parameters, fitting is equivalent to minimizing

2\sum\_{ij}\frac{w\_{ij}}{\phi}D\_\Psi(y\_{ij},(HS+BA+O)\_{ij})
+\lambda_RK

under the nonnegativity, profile-normalization and log-concavity
constraints. The saturated likelihood contributes a constant in these
conditional model comparisons. The corresponding divergences include:

| Noise model | Generator \Psi(u) | Half unit deviance D\_\Psi(y,\mu) |
|:---|:---|:---|
| Gaussian | u^2/2 | (y-\mu)^2/2 |
| Poisson | u\log u-u | y\log(y/\mu)-y+\mu |
| Gamma, uncensored y\>0 | -\log u | y/\mu-\log(y/\mu)-1 |
| NB2, fixed shape \theta | u\log u-(u+\theta)\log(u+\theta) | y\log(y/\mu)-(y+\theta)\log\\(y+\theta)/(\mu+\theta)\\ |

For NB2, \phi=1 and the shape \theta enters the generator. Gaussian,
Poisson and Gamma therefore give squared-error, generalized KL and
Itakura–Saito NMF losses. Quasipoisson gives the KL kernel scaled by its
dispersion. Censored Gamma and the mean-dependent-shape extensions use
their full likelihood kernels described below; estimated noise
parameters are compared using the full likelihood as well.

## Noise models and likelihood kernels

| Noise model | Variance | Negative log-likelihood score in the mean |
|:---|:---|:---|
| Gaussian | \phi | (\mu-y)/\phi |
| Poisson | \mu | 1-y/\mu |
| Quasipoisson | \phi\mu | (1-y/\mu)/\phi |
| NB2 | \mu+\mu^2/\theta | \theta(\mu-y)/\\\mu(\mu+\theta)\\ |
| Gamma | \phi\mu^2 | (\mu-y)/(\phi\mu^2) |
| NB1 | \phi\mu | \\\psi(a)-\psi(a+y)+\log\phi\\/(\phi-1), a=\mu/(\phi-1) |
| Common-scale Gamma | b_0\mu | \\\psi(\mu/b_0)+\log b_0-\log y\\/b_0 |

The final two models use the exact distributions documented in
[`signal_families`](https://tomwenseleers.github.io/L0glm/reference/signal_families.html).
Their observed positive curvatures are differences of trigamma functions
for NB1 and \psi_1(\mu/b_0)/b_0^2 for common-scale Gamma. The latter is
convex in the mean; the NB1 curvature is nonnegative because trigamma
decreases on the positive real line.

## Censored Gamma intensities

Recorded zeros may represent intensities below a detection limit d,
whose contribution is \log P(Y\<d;\mu). Positive observations contribute
their ordinary density. For standard Gamma, set a=1/\phi, z=d/(\phi\mu);
its negative log-CDF derivative is

g(\mu)=\frac{z^a e^{-z}}{\Gamma(a)P(a,z)\mu}.

For common-scale Gamma, z=d/b_0 stays fixed and a=\mu/b_0 changes. The
convergent incomplete-Gamma series

P(a,z)=\frac{e^{-z}z^a}{\Gamma(a+1)}
\left\\1+\sum\_{k\geq1}\frac{z^k}{(a+1)\cdots(a+k)}\right\\

provides its shape derivative when z is small. The general branch uses
the Gamma CDF and a numerical shape derivative. Both branches are
checked against the actual censored likelihood. The default detection
limit is 0.5 intensity units and can be changed through
`control$censor`.

## Local quadratic updates

For a coefficient b\geq0, keep the rest of the signal in an offset r_i
and write \mu_i=r_i+x_i b. With negative log-likelihood score g_i and
positive working curvature I_i, propose

b^+=\max\left\\0,b-\frac{\sum_i w_i x_i g_i} {\sum_i w_i x_i^2
I_i}\right\\.

Actual-loss backtracking and a zero-boundary check complete each scalar
solve. The Poisson branch uses a faster bracketed likelihood solve. Its
score \sum_iw_ix_i(1-y_i/\mu_i) is monotone, and its nonnegative root is
bounded above by \sum_iw_iy_i/\sum_iw_ix_i. Channels are independent
conditional on the profiles, allowing deterministic OpenMP parallelism.

Local profile updates use the original observation likelihood. For
\eta_i=\log h_i, \mu\_{ij}=r\_{ij}+h_i s_j, the chain-rule score and
working information are

g\_{\eta,i}=h_i\sum_jw\_{ij}s_jg\_{\mu,ij},\qquad
I\_{\eta,i}=h_i^2\sum_jw\_{ij}s_j^2I\_{\mu,ij}.

The working target is \eta_i-g\_{\eta,i}/I\_{\eta,i}, with tiny
stabilization where the information vanishes. This quadratic
approximation updates the **log profile parameters**; the observed
intensities remain on their original scale and acceptance uses their
actual likelihood.

A smooth concave log profile can be written as

\eta(t)=a+bt-\frac12\sum_j c_j \\(t-u_j)\_+^2-(t-v_j)\_+^2\\,\qquad
c_j\geq0,

where the intervals \[u_j,v_j) partition the local time window. Its
second derivative is -\sum_jc_j\mathbf1\_{\[u_j,v_j)}(t)\leq0, giving
hard log-concavity. The intercept and slope are unrestricted, and a
small nonnegative quadratic solve estimates the curvature coefficients.
For more detailed local adjustments, piecewise-linear log-profile
increments preserve the existing fine-scale curvature; a
concave-envelope projection removes numerical constraint residuals.
Backtracking between feasible log profiles preserves concavity and
accepts only likelihood improvements.

Fusion starts with spectral compatibility and profile overlap.
Broad-anchor and contained-tail proposals are refined only within the
affected windows; other compounds remain in an offset. The full
original-channel criterion decides whether each proposed rank change is
retained. A final conservative pass can reunite broad profiles with
complementary ion sets when their profile cosine is at least 0.95.
Neighbouring modes are guarded during its local refit.
`control=list(profile_cycles=0L)` disables the main shape-refinement
stage; `control=list(profile_fusion=FALSE)` disables the final
conservative pass.

### Complete the fit on the original ion channels

The final pass checks conditional deletions, fusions of overlapping
profiles, and PSF candidates suggested by the original-channel
residuals. Nearly proportional profiles can be considered even when
neither is a broad peak. Candidate construction uses matrix products and
performs no refitting. Accepted proposals use native local likelihood
fits with all other compounds and the fitted background held in the
offset.

For J=-2\ell+\lambda_RK, a deletion has change \Delta J=\Delta
D-\lambda_R, a fusion of r factors has \Delta J=\Delta D-\lambda_R(r-1),
and an addition has \Delta J=\Delta D+\lambda_R. The complete
original-channel criterion must improve. Noise parameters, offsets,
weights and the numerical compound penalty stay fixed. Cached deletion
scores need updating only for profiles whose temporal supports overlap
the deleted factor.

The pass uses bounded scan windows to screen missed peaks, then refines
the accepted profiles locally under hard log concavity. The same
procedure is available for every supported noise family; working-score
screening proposes candidates and each family’s likelihood decides
acceptance. The returned `fit$cleanup` records the moves and full
objective change. `control=list(original_cleanup=FALSE)` skips this
pass.

## Fusion, aggregation and rank selection

1.  Obtain a variance-matched grouped deconvolution start, retaining the
    spline background as compulsory unpenalized terms.
2.  Compare normalized spectra by Jensen–Shannon information loss.
    Propose fusions when their elution profiles overlap and their
    combined profile can be made log-concave. Clustering itself performs
    no regression fits.
3.  Refine each proposal locally, holding all outside contributions
    fixed. Overlapping spectra and the candidate profile can adjust
    together.
4.  Accept a rank change when
    -2(\ell\_{new}-\ell\_{old})+\lambda_R(K\_{new}-K\_{old})\leq0.
5.  Update broad profiles, smooth background coefficients and spectra,
    followed by final boundary checks and rank pruning.

Algebraic fusion uses a least concave majorant of the log profile,
followed by normalization. Local likelihood updates impose hard
concavity constraints in their quadratic solves. Every complete accepted
stage decreases the actual criterion, within numerical tolerance.
`history` and `ledger_error` make this property directly checkable.

With `precluster=32L`, the initial workflow uses at most 32 aggregated
outcome groups. Size factors `colSums(pmax(Y,0))` account for different
total ion intensities. Final refinement returns to every original
channel. Independent NB1 variables sharing p=1/\phi are closed under
sums, because their probability-generating functions multiply as

\prod_j\left\\\frac{p}{1-(1-p)t}\right\\^{\mu_j/(\phi-1)}
=\left\\\frac{p}{1-(1-p)t}\right\\^{\sum_j\mu_j/(\phi-1)}.

The same property holds for common-scale Gamma latent signals, through
\prod_j(1+b_0t)^{-\mu_j/b_0}=(1+b_0t)^{-\sum_j\mu_j/b_0}. Standard NB2
and Gamma use moment matching in the reduced problem. Censoring or
different channel weights also give working aggregate approximations;
final scoring always uses the original responses and their specified
model.

## Temporal variation

High-order differences suppress the slowly varying mean. With
a_k=(-1)^k{r\choose k} and stationary covariance c_h,

\frac{\operatorname{Var}(\Delta^r e)}{{2r\choose r}}
=c_0+2\sum\_{h=1}^r(-1)^h \frac{{2r\choose r+h}}{{2r\choose r}}c_h.

The identity follows by expanding the variance and using \sum_k a_k
a\_{k+h}=(-1)^h{2r\choose r+h}. With lag-one covariance this reduces to
c_0-2r c_1/(r+1). Fourth and eighth differences estimate the marginal
variance and lag-one covariance. The optional penalty multiplier
\kappa=1+2\widehat\rho_1 approximates the local information loss for
smooth directions, while the marginal scale stays in the component
likelihood. This is a scalar working covariance adjustment.
Composite-likelihood theory relates the broader correction to
\operatorname{tr}(JH^{-1}), where H is sensitivity and J is score
covariance (Varin, Reid & Firth, 2011).

## Single-sample algorithm

Let L=-\sum\_{ij}w\_{ij}\ell(y\_{ij};\mu\_{ij},\nu). Resolve the pooled
noise parameters \nu and information multiplier \kappa before rank
comparisons. The default initializer uses a positive smoothed mean
\widetilde\mu and variance-matched Poisson working weights

w^{\rm pilot}\_{ij}=w\_{ij}\\
\frac{\widetilde\mu\_{ij}}{\operatorname{Var}(Y\_{ij}\mid\widetilde\mu\_{ij},\nu)}.

It fits grouped nonnegative
[`L0glm.fit()`](https://tomwenseleers.github.io/L0glm/reference/L0glm.fit.html)
on the PSF design with compulsory background, using \lambda_R/m. This
supplies candidate components; subsequent fits and acceptance use the
requested family. A supplied `start` with profiles, spectra and optional
background coefficients bypasses screening.

**Algorithm N1 — sparse log-concave NMF with rank selection**

1.  Normalize candidate PSF columns. Resolve the noise scale, optional
    lag-one information adjustment, weights, offsets and background
    basis.
2.  If preclustering is requested:
    1.  Form size-adjusted outcome groups and sum their
        observations/offsets.
    2.  Construct the corresponding exact or moment-matched aggregate
        noise model.
    3.  Initialize and fuse components in that reduced space, retaining
        the original channel count in \lambda_R.
    4.  Lift loadings/background with group intensity fractions; return
        to every original channel for likelihood refinement.
3.  Otherwise obtain the variance-matched grouped deconvolution start
    directly. Normalize each profile and transfer its scale into its
    loadings.
4.  Run the bounded local refinement stages:
    1.  Hold other components and background fixed in a local offset.
        Update each spectrum by nonnegative scalar likelihood solves,
        including zero-boundary checks; update each concave log profile
        with the quadratic subroutine below.
    2.  Update nonnegative background coefficients conditionally on the
        component signal. Preserve the specified family likelihood.
    3.  Generate spectral/overlap-compatible fusion candidates
        algebraically, starting with stronger components. Refine only
        the affected scan windows.
    4.  Accept whole candidates only if their complete J decreases
        within numerical tolerance. Check conditional deletions and zero
        spectra.
    5.  If NB2 shape updating is enabled, propose bounded steps in
        \log\theta and backtrack until the original likelihood improves.
5.  Refine full-channel shapes, broad/contained-tail fusions and shared
    local backgrounds. Finally check original-channel deletions,
    compatible profile fusions and additions suggested by residual
    working scores.
6.  Return temporally ordered H,S, background, family/noise parameters,
    likelihood, rank criterion, stage acceptance ledger and timings.

`profile_cycles`, `fusion_passes` and the other stage controls bound the
work; local solvers also stop when their score/step change is
sufficiently small. For fixed noise parameters, every accepted stage
uses the same original-channel objective. The reduced-space search
supplies an initializer, followed by original-channel rank checks.

**Subroutine Q — local likelihood and shape update**

1.  Subtract the focus component from the current mean to obtain the
    fixed offset.
2.  Compute mean-scale likelihood scores and positive working
    information.
3.  For a loading b, propose a nonnegative scalar scoring step or use
    the bracketed Poisson root. For a profile, set \eta=\log h, apply
    the chain-rule scores and solve a weighted quadratic under
    nonincreasing slopes.
4.  Halve the step until the original likelihood improves and
    constraints hold. A convex combination of feasible log profiles
    remains concave.
5.  Normalize the profile, transfer its scale to loadings, and update
    only the affected mean block. Outside-window observations are
    unchanged.

For a rank change \Delta K, acceptance is based on

\Delta J=2(L\_{\rm new}-L\_{\rm old})+\lambda_R\Delta K.

Because outside-window means are unchanged, their likelihood
contributions cancel exactly. Local computation therefore gives the
complete objective difference without rescoring unaffected observations.

## Coupled model on original scan grids

For sample s, let u=w_s(t) be its monotone alignment map. Native bin
edges t\_{si-1/2},t\_{si+1/2} define

\Delta u\_{si}=w_s(t\_{si+1/2})-w_s(t\_{si-1/2}),\qquad
J\_{si}=\frac{\Delta u\_{si}}{\Delta t\_{si}}.

With known exposure e\_{si}\>0, the mean model is

\mu\_{sij}=e\_{si}\left\\J\_{si}\sum_g a\_{sg}h\_{sg}(u\_{si})p\_{gj}
+(B_sA_s)\_{ij}\right\\+o\_{sij}.

Spectra have p\_{gj}\geq0 and \sum_jp\_{gj}=1. Each detected profile
satisfies

h\_{sg}\geq0,\quad \sum_i\Delta u\_{si}h\_{sg}(u\_{si})=1,\quad
\frac{\log h\_{sg,i+1}-\log h\_{sg,i}}{u\_{s,i+1}-u\_{si}} \text{ is
nonincreasing}.

The Jacobian preserves integrated area under alignment. No interpolation
or resampling is applied to observations. Local likelihood evaluations
use their original channels, weights, offsets and sample-specific noise
parameters.

The rank-selection objective is

\mathcal J=-2\sum\_{sij}w\_{sij}\ell_s(y\_{sij};\mu\_{sij})
+\kappa\lambda\_{IC}\\(m-1)G+L+d\_{\rm noise}\\.

Here G is the number of common spectra and L the number of detected
sample-compound occurrences. A unit-sum spectrum contributes m-1 free
coordinates and each occurrence contributes an amplitude. Profile and
baseline complexity are unpenalized shape-constrained nuisance terms in
this rank criterion. For temporal-correlation adjustment, \kappa is the
largest resolved sample information factor. Quasipoisson uses its
frozen, scaled Poisson kernel. The same
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html) and
[`GIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html)
interfaces apply.

## Alignment and local updates

1.  Fit samples independently, optionally initializing in 32 aggregated
    channel groups and returning to every original channel.
2.  Match high-confidence spectra to a reference sample. Use the
    1%-height peak start plus half the narrow-peak width as the timing
    anchor.
3.  Fit a robust low-dimensional drift spline, with log area ratio as a
    nuisance term for concentration-related anchor motion. Retain its
    drift component and enforce positive slope.
4.  Group spectra algebraically when their timing and area-ordered
    profile overlap are compatible. Each group has at most one
    occurrence per sample. This clustering stage performs no NMF or GLM
    refits.
5.  Update local concave log profiles, pool spectral scores over
    samples, and update each smooth background. Backtrack against the
    original likelihood.
6.  Retain a grouped fit when its complete penalized objective improves.
    Return temporally ordered consensus spectra and sample-specific
    areas.

For a common spectral coefficient p\_{gj}, with design value
x\_{sgi}=e\_{si}J\_{si}a\_{sg}h\_{sg}(u\_{si}), the pooled
negative-log-likelihood score and Fisher information are

U\_{gj}=\sum\_{si:g\in s}w\_{sij}x\_{sgi}g\_{\mu,sij},\qquad
I\_{gj}=\sum\_{si:g\in s}w\_{sij}x\_{sgi}^{2}I\_{\mu,sij}.

A nonnegative scoring step is followed by full-likelihood backtracking.
Renormalizing a spectrum transfers its scale into the occurrence
amplitudes and preserves the fitted signal. Profile updates use the
corresponding chain-rule score in log coordinates and linear slope
constraints on actual aligned times. Native OpenMP updates give
identical results across thread counts in the validation examples.

The native spectral kernel evaluates each cell’s likelihood derivatives
once and pools them over all sparse profiles in a single threaded pass.
On a beer chromatogram with 14,949 scans, 518 channels and 325 profiles,
this reduced one pooled-score calculation from **5.45 s to 0.19 s**, a
**28.7-fold speedup**, while matching the profile-wise calculation to
numerical precision. Cached peak geometry also avoids rescanning every
profile after each algebraic merge. Shared row-wise weights and offsets
retain one vector throughout the coupled fit, with channel-specific
matrices used when required.

All seven [signal
families](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.md)
are supported: Gaussian, Poisson, quasipoisson, NB2, NB1, censored Gamma
and common-scale Gamma. Dispersion parameters are resolved from the
original channels and frozen for group comparisons. NB2 additionally
supports likelihood-improving shape updates when initialized with a
negative `theta`.

## Coupled algorithm and normalization

For coupled models use n\_\*=p\_\*=\max_s n_s in the IC preset, and the
largest sample information multiplier as \kappa. Noise parameters are
shared across channels within each sample. The joint rank count is
(m-1)G+L+d\_{\rm noise}, reducing to mK+d\_{\rm noise} when one sample
has G=L=K. Common unit-sum spectra contribute m-1 coordinates each; each
observed sample-compound occurrence contributes one amplitude.

If \sum_i\Delta u\_{si}h\_{sg}(u\_{si})=1, then

\sum_i\Delta t\_{si}J\_{si}h\_{sg}(u\_{si}) =\sum_i\Delta
u\_{si}h\_{sg}(u\_{si})=1.

Thus an occurrence amplitude represents integrated intensity in the
chosen time unit, before its known exposure is applied.
[`sampleNMF()`](https://tomwenseleers.github.io/L0nmf/reference/sampleNMF.md)
also returns unit-sum original-grid profiles and their compatible
discrete loadings.

**Algorithm N2 — coupled multisample log-concave NMF**

1.  Fit samples independently by N1, or read supplied fitted factors.
    Persist observations and factor checkpoints sample by sample.
2.  Match high-confidence spectra to a reference. Obtain
    onset-plus-half-width anchors and fit a robust monotone drift map,
    accounting for log area ratio as a nuisance predictor of anchor
    motion. Retain the drift component.
3.  Map each original bin’s edges to aligned time. Compute \Delta u,J,
    initialize log-concave profiles there, normalize by \Delta u, and
    put their scale into occurrence amplitudes. Initialize each
    occurrence with its own spectrum.
4.  Perform `initial_cycles` of local profile, pooled spectrum and
    background updates; retain this ungrouped-spectrum reference
    checkpoint.
5.  For up to `rounds` grouping rounds:
    1.  Restrict spectrum comparisons to nearby aligned time bins.
        Compare cosine or Hellinger distance, with at most one
        occurrence per sample per group.
    2.  Check area-ordered log-concave-ridge compatibility. Pool
        compatible spectra by occurrence area. Keep profiles fixed
        during this algebraic grouping.
    3.  Run up to `cycles` joint local-update cycles: constrained sample
        profiles; batched, pooled spectral scores; nonnegative sample
        backgrounds; optional NB2 shape updates. Backtrack each update
        on original-channel likelihood.
    4.  Stop grouping when a subsequent round makes no new joins.
6.  Compare the complete grouped/refined candidate against the reference
    checkpoint using the joint criterion. Retain the lower-criterion
    model.
7.  Compact consensus IDs in temporal order. Return shared spectra,
    sample factors, alignment, diagnostics and an area table with
    missing occurrences `NA`.

Clustering constructs a complete candidate before joint refinement. The
final comparison includes every original channel and all detected
occurrences.

### Area-ordered ridge compatibility

Put matched profiles on a small comparison grid for geometry only;
observed signals remain on their native grids. Order profiles by total
area and allow bounded local alignment. Their thresholded level sets
must overlap and support a common log-concave ridge across the area
ordering. This screens concentration- dependent broadening and tailing
before sharing a spectrum. The subsequent original-data likelihood
decides whether the combined proposal is retained.

### Streaming and parallel computation

Sparse profile supports restrict work to affected scan windows. The
coupled spectral kernel evaluates each observation’s likelihood
derivatives once and pools them across its active profiles. Channels can
be processed independently in OpenMP; sample contributions are
accumulated in a fixed order. The main optimization reads one sample at
a time, with shared spectra and compact geometry summaries in memory.
Shared row weights/offsets retain one vector. RT bins avoid an
all-occurrence pairwise similarity matrix.

## Methodological foundations and implementation

Sparse pulse-train regression builds on de Rooi & Eilers, and the
adaptive ridge initialization uses L0glm’s extensions of reweighted
quadratic selection. Eilers’ work motivates positive shape-constrained
profiles; generalized and deviance matrix factorization provide the
likelihood framework. Bregman and size-adjusted Poisson clustering
supply useful spectral and channel geometries.

The likelihood-based hard shape solves, complete original-channel rank
moves, windowed kernel checkpoints and coupled normalization/acceptance
rules above describe the L0nmf algorithms. The coupled criterion counts
shared spectra and sample-specific amplitudes explicitly. The general
references establish the underlying frameworks; the objectives and
pseudocode here specify their adaptation to these constrained signal
models.

## References

Kidziński, Ł., Hui, F. K. C., Warton, D. I. & Hastie, T. J. (2022).
[Generalized matrix factorization: efficient algorithms for fitting
generalized linear latent variable models to large data
arrays.](https://www.jmlr.org/papers/v23/20-1104.html)

Castiglione, C., Segers, A., Clement, L. & Risso, D. (2026). [Stochastic
gradient descent estimation of generalized matrix factorization models
with application to single-cell RNA sequencing
data.](https://doi.org/10.1093/biostatistics/kxag010)

Wang, L. & Carvalho, L. E. (2023). [Deviance matrix
factorization.](https://doi.org/10.1214/23-EJS2174)

Li, X. & Tyagi, A. (2023). [Block-active ADMM to minimize NMF with
Bregman divergences.](https://doi.org/10.3390/s23167229)

Banerjee, A. et al. (2005). [Clustering with Bregman
divergences.](https://jmlr.org/papers/v6/banerjee05b.html)

Witten, D. M. (2011). [Classification and clustering of sequencing data
using a Poisson model.](https://doi.org/10.1214/11-AOAS493)

Lee, D. D. & Seung, H. S. (2000). [Algorithms for non-negative matrix
factorization.](https://proceedings.neurips.cc/paper_files/paper/2000/file/f9d1152547c0bde01830b7e8bd60024c-Paper.pdf)

Eilers, P. H. C. (2005). [Unimodal
smoothing.](https://doi.org/10.1002/cem.935)

Varin, C., Reid, N. & Firth, D. (2011). [An overview of composite
likelihood
methods.](https://www3.stat.sinica.edu.tw/statistica/oldpdf/A21n11.pdf)

The likelihood distributions follow the R
[negative-binomial](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/NegBinomial.html)
and
[Gamma](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/GammaDist.html)
parameterizations. The compound-count criterion, local log-concave
fusion and their integration with L0glm’s grouped initialization are
developments of this package.
