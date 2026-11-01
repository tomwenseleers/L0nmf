# Log-concave spline smoothing of a discrete convolution kernel estimate

Fits a concave quadratic B-spline (via
[`cobs`](https://rdrr.io/pkg/cobs/man/cobs.html)) to \\\log\_{10}(y +
\varepsilon)\\ and back-transforms and normalises the fit so the result
is a \*\*log-concave, non-negative, unit L1 norm\*\* smoothed version of
the original convolution kernel estimate \\y\\. A small offset is added
internally to avoid taking the logarithm of zero; the final vector is
renormalised to unit L1 norm.

## Usage

``` r
log_conc_smooth(y, tau = 0.5, nknots = NULL)
```

## Arguments

- y:

  numeric vector (length \\n \ge 3\\) with the raw, non-negative data
  sampled at the integer points `1:length(y)`.

- tau:

  quantile to be fitted by `cobs`; defaults to `0.5` (median fit).

- nknots:

  integer; number of interior knots for the spline. If `NULL` (default)
  it is set to `round(length(y)/5)`.

## Value

A numeric vector of length `length(y)`:

- strictly non-negative,

- log-concave over its support,

- rescaled to unit L1 norm.

## Details

The weights passed to `cobs()` are proportional to `y + offs`, where
`offs = max(y) * 1e-6`; in practice that gives the most stable fits for
sparse or strongly peaked inputs.

## See also

[`cobs`](https://rdrr.io/pkg/cobs/man/cobs.html) for the underlying
constrained B-spline regression.

## Author

Tom Wenseleers <tom.wenseleers@kuleuven.be> \\ Laboratory of
Socioecology and Social Evolution, KU Leuven \\
<https://bio.kuleuven.be/eeb/tw/research> \\ ORCID:
[0000-0002-1434-861X](https://orcid.org/0000-0002-1434-861X)

## Examples

``` r
library(L0glm)
if (requireNamespace("cobs", quietly = TRUE)) {
set.seed(1)
x      <- seq(-5, 5, length.out = 101)

## underlying "true" density on the grid
f_true <- dnorm(x)
f_true <- f_true / sum(f_true)          # ensure it integrates to 1

N      <- 100000                        # total expected counts
lambda <- N * f_true                    # Poisson means in each bin

y_cnt  <- rpois(length(x), lambda)      # Poisson-noisy counts
y      <- y_cnt / sum(y_cnt)            # rescale back to a density

g <- log_conc_smooth(y)

matplot(cbind(y, f_true, g), type = "l", lty = 1,
        col = c(2, 1, 4), ylab = "density",
        main = "Raw (Poisson) vs. log-concave smoothed curve + ground truth")
legend("topleft", legend = c("raw (Poisson)", "smoothed", "ground truth"),
       col = c(2, 4, 1), lty = 1, bty = "n")
}
```
