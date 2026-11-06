# Estimate a shared noise model for ordered multichannel signals

Estimate a shared noise model for ordered multichannel signals

## Usage

``` r
estimateSignalNoise(
  y,
  family = quasipoisson("identity"),
  dispersion = NULL,
  theta = NULL,
  weights = NULL,
  correlation = "none",
  censor = 0.5,
  ncores = NULL
)
```

## Arguments

- y:

  Observations-by-channels matrix.

- family:

  Identity-link family accepted by \[splogconNMF()\].

- dispersion, theta:

  Fixed noise parameters, or NULL to estimate them.

- weights:

  Scalar, row vector or matrix of observation precisions.

- correlation:

  Correlation convention used by \[splogconNMF()\].

- censor:

  Gamma detection limit.

- ncores:

  Number of physical worker cores; NULL detects them.

## Value

A list containing the fitted noise specification, pilot means and
equivalent Poisson working precision for PSF diagnostics.

## Details

Pooled marginal variances are estimated from temporal differences. A
shared information factor accounts for estimated lag-one correlation.
The same calibration can be reused during PSF estimation and
factorization.
