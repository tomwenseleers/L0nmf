#' Fit a log-concave factorization with automatic rank selection
#'
#' Estimate ordered, nonnegative log-concave profiles and channel loadings from
#' one signal or several related samples. The information criterion selects the
#' number of components. Coupled fits share spectra while retaining the original
#' scan grids and estimating a separate profile and area for each occurrence.
#'
#' @param y A numeric scans-by-channels matrix or vector for one sample, or a
#'   list of at least two sample descriptions for a coupled fit. Each description
#'   contains `y` and `time` and may include its own `X`, `background`,
#'   `weights`, `offset`, `exposure` or existing `fit`. Matrices and initial
#'   fits in a description can be supplied as RDS paths. Samples must share
#'   channel names and order. List names supply sample IDs when a description
#'   omits `sample_id`. See [coupledNMF()] for all descriptor fields.
#' @param X Candidate profile matrix for a single sample, such as the sparse
#'   design from [estimate_psf()]. In a coupled fit, supply it in each sample
#'   description or let the coupled initializer estimate it.
#' @param lambda An information criterion (`"gic"`, `"aic"`, `"bic"`,
#'   `"ebic"` or `"mbic"`) or a nonnegative numeric penalty per selection
#'   degree of freedom, before any temporal information adjustment. With
#'   `m` channels, a single component incurs `m * lambda`; in a coupled fit
#'   a shared spectrum incurs `(m - 1) * lambda` and an occurrence incurs
#'   `lambda`.
#' @param precluster Maximum channel groups or a saved [channel_groups()] map.
#'   `"auto"` uses the original channels for a single sample and up to 32
#'   common groups for a coupled fit. `FALSE` or `NULL` disables aggregation.
#'   Final fits use every original channel.
#' @param correlation Lag-one working information adjustment. `NULL` uses the
#'   fitter's default (`"auto"` for one sample and `"none"` for coupled fits),
#'   or `control$correlation` if supplied for a coupled fit.
#' @param control Named controls for the selected fitter. See [splogconNMF()]
#'   and [coupledNMF()] for their respective controls.
#' @param ncores Total core budget. `NULL` uses the detected physical core count.
#' @param object A fitted model returned by `logconNMF()`.
#' @param sample Sample index or ID. For a single matrix use `1L` or
#'   `"sample1"`. In `coef()`, `NULL` returns the common spectra for a
#'   coupled fit.
#' @param normalize `"none"` returns amplitude-scaled loadings. `"l1"`
#'   returns unit-sum rows.
#' @param ... Further arguments passed to the underlying extractor.
#' @inheritParams splogconNMF
#' @inheritParams coupledNMF
#' @details
#' For a coupled fit, put `background`, `weights`, `offset` and any initial
#' `fit` in the sample descriptions. `path = NULL` uses a new temporary
#' sample store. Supply a persistent path when saving a coupled fit for later
#' use, and retain the store with the fitted object.
#'
#' Numeric penalties in `logconNMF()` have the same scale in both modes.
#' To reproduce a single-sample [splogconNMF()] call using a numeric compound
#' penalty, divide that penalty by the number of channels.
#' @return A fitted object with estimated `rank`, the number of single-sample
#'   components or shared components. Use [coef()], [fitted()] and [GIC()] in
#'   either mode. `fitted()` returns one sample on its original grid. A coupled
#'   fit also provides its sample-by-compound `areas` and [sampleNMF()] views.
#' @examples
#' X <- convmatrix(convkernel("gaussian", 9L), 50L)
#' signal <- cbind(700 * X[, 16L] + 350 * X[, 35L],
#'                 150 * X[, 16L] + 500 * X[, 35L])
#' set.seed(12)
#' Y <- matrix(rpois(length(signal), signal), nrow(signal))
#' fit <- logconNMF(Y, X = X, ncores = 1L)
#' fit$rank
#' coef(fit, normalize = "l1")
#'
#' second <- matrix(rpois(length(signal), 1.4 * signal), nrow(signal))
#' samples <- list(first = list(y = Y, time = seq_len(nrow(Y)), X = X),
#'   second = list(y = second, time = seq_len(nrow(Y)), X = X))
#' joint <- logconNMF(samples, ncores = 1L)
#' joint$rank
#' joint$areas
#' fitted(joint, sample = "second")[1:3, ]
#' @md
#' @export
logconNMF <- function(y, X = NULL, family = poisson("identity"), lambda = "gic",
                      background = NULL, dispersion = NULL, theta = NULL,
                      precluster = "auto", weights = NULL, offset = NULL,
                      correlation = NULL, ncores = NULL, control = list(),
                      start = NULL, path = NULL, align = TRUE,
                      method = c("cosine", "hellinger"), verbose = FALSE) {
  if (is.numeric(lambda) && (length(lambda) != 1L || is.complex(lambda) ||
      !is.finite(lambda) || lambda < 0)) {
    stop("Use one finite, nonnegative penalty or an information criterion.", call. = FALSE)
  }
  if (!is.list(control)) stop("Use a named control list.", call. = FALSE)
  coupled <- is.list(y) && !is.data.frame(y)
  if (identical(precluster, "auto")) precluster <- if (coupled) 32L else NULL
  if (coupled) {
    if (length(y) < 2L) stop("Supply a matrix for one sample or at least two sample descriptions.", call. = FALSE)
    if (!is.null(X) || !is.null(background) || !is.null(weights) ||
        !is.null(offset) || !is.null(start)) {
      stop("Supply X, background, weights, offset and initial fits in the sample descriptions.",
        call. = FALSE)
    }
    samples <- lapply(seq_along(y), function(i) {
      sample <- y[[i]]
      if (!is.list(sample) || is.null(sample$y) || is.null(sample$time)) {
        stop("Each sample description needs y and time.", call. = FALSE)
      }
      name <- names(y)[i]
      if (is.null(sample$sample_id) && length(name) && !is.na(name) && nzchar(name)) {
        sample$sample_id <- name
      }
      sample
    })
    if (!is.null(correlation)) control$correlation <- correlation
    if (is.null(path)) path <- tempfile("logconNMF_")
    fit <- coupledNMF(samples, family = family, lambda = lambda, path = path,
      precluster = precluster, dispersion = dispersion, theta = theta,
      method = match.arg(method), align = align, ncores = ncores, control = control,
      verbose = verbose)
  } else {
    if (is.null(X)) stop("Supply X for a single-sample fit.", call. = FALSE)
    if (!is.null(path) || !missing(align) || !missing(method)) {
      stop("path, align and method apply to coupled fits.", call. = FALSE)
    }
    if (is.numeric(lambda)) lambda <- lambda * NCOL(y)
    if (is.null(correlation)) correlation <- "auto"
    fit <- splogconNMF(X, y, family = family, background = background, lambda = lambda,
      dispersion = dispersion, theta = theta, precluster = precluster, weights = weights,
      offset = offset, correlation = correlation, ncores = ncores, control = control,
      start = start, verbose = verbose)
  }
  class(fit) <- c("logconNMF_fit", class(fit))
  fit
}

logconNMF_check_sample <- function(sample) {
  if (!is.null(sample) && !isTRUE(sample %in% c("1", "sample1"))) {
    stop("Select sample 1 or 'sample1' for a single-sample fit.", call. = FALSE)
  }
}

#' @rdname logconNMF
#' @export
coef.logconNMF_fit <- function(object, sample = NULL, normalize = c("none", "l1"), ...) {
  normalize <- match.arg(normalize)
  if (inherits(object, "coupledNMF")) {
    return(coef.coupledNMF(object, sample = sample, normalize = normalize, ...))
  }
  logconNMF_check_sample(sample)
  loadings <- coef.splogconNMF(object, ...)
  if (normalize == "l1" && nrow(loadings)) {
    totals <- rowSums(loadings)
    totals[totals == 0] <- 1
    loadings <- loadings / totals
  }
  loadings
}

#' @rdname logconNMF
#' @export
fitted.logconNMF_fit <- function(object, sample = 1L, ...) {
  if (inherits(object, "coupledNMF")) return(fitted.coupledNMF(object, sample = sample, ...))
  logconNMF_check_sample(sample)
  fitted.splogconNMF(object, ...)
}
