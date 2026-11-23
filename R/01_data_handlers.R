# ── Data Ingestion: Adapter Pattern ───────────────────────────────────────────
# ensure_computable() normalizes any time-series input to a clean numeric vector
# plus metadata. Results stay numeric; the metadata (tsp / idx) is kept in
# `meta` so plots and future as.ts()/as.xts() methods can map them back.

# ── Normalizer ───────────────────────────────────────────────────────────────

#' Normalize a time-series input to a clean numeric vector
#'
#' Detects the class of `x`, extracts values and temporal metadata, and
#' validates that the series contains no `NA` or `Inf` values.
#'
#' Dispatch uses `inherits()`, so subclasses work: `zooreg` (what
#' `zoo::as.zoo()` returns for a `ts`) is handled as `zoo`, and `xts` is
#' tested first because it inherits from `zoo`. Filters are univariate: any
#' input with more than one column (`mts`, a multi-column `xts`/`zoo`, a
#' matrix) is an error instead of being flattened into one long vector. A
#' single-column matrix is accepted as a vector.
#'
#' @param x A `numeric`, `ts`, `xts`, or `zoo` object.
#' @return A list with components:
#'   \describe{
#'     \item{y}{Clean `double` vector of values.}
#'     \item{idx}{Time index vector or `NULL` (for `numeric`/`ts`).}
#'     \item{tsp}{`tsp` attribute or `NULL` (only for `ts`).}
#'     \item{class}{Character: the original class name.}
#'   }
#' @keywords internal
#' @noRd
ensure_computable <- function(x) {
  if (!is.null(dim(x)) && NCOL(x) > 1L) {
    stop(
      sprintf(
        "`x` has %d columns; filters are univariate. Pass one column (e.g. x[, 1]).",
        NCOL(x)
      ),
      call. = FALSE
    )
  }

  handler <- if (inherits(x, "xts")) {
    .handle_xts(x)                       # before zoo: xts inherits from zoo
  } else if (inherits(x, "zoo")) {
    .handle_zoo(x)                       # includes zooreg
  } else if (inherits(x, "ts")) {
    .handle_ts(x)
  } else if (is.numeric(x) && is.null(attr(x, "class"))) {
    .handle_numeric(x)                   # numeric, integer, 1-column matrix
  } else {
    # data.frame / data.table / character / anything else -> error
    stop(
      sprintf(
        "Unsupported input class '%s'. Pass a numeric vector, ts, xts, or zoo object.\n  Hint: for data.frames use df$column.",
        class(x)[1L]
      ),
      call. = FALSE
    )
  }

  .validate_values(handler$y)
  handler
}

# ── Class-specific handlers (internal) ───────────────────────────────────────

.handle_numeric <- function(x) {
  list(
    y        = as.double(x),
    idx      = NULL,
    tsp      = NULL,
    class    = "numeric"
  )
}

.handle_ts <- function(x) {
  list(
    y        = as.double(x),
    idx      = NULL,
    tsp      = stats::tsp(x),
    class    = "ts"
  )
}

.handle_xts <- function(x) {
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("Package 'zoo' is required to handle xts objects.", call. = FALSE)
  }
  list(
    y        = as.double(zoo::coredata(x)),
    idx      = zoo::index(x),
    tsp      = NULL,
    class    = "xts"
  )
}

.handle_zoo <- function(x) {
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("Package 'zoo' is required to handle zoo objects.", call. = FALSE)
  }
  list(
    y        = as.double(zoo::coredata(x)),
    idx      = zoo::index(x),
    tsp      = NULL,
    class    = "zoo"
  )
}

# ── Value validation ─────────────────────────────────────────────────────────

.validate_values <- function(y) {
  n_nan <- sum(is.nan(y))
  n_na  <- sum(is.na(y)) - n_nan   # is.na() counts NaN too
  n_inf <- sum(is.infinite(y))

  issues <- character(0L)
  if (n_na > 0L)  issues <- c(issues, sprintf("%d NA", n_na))
  if (n_nan > 0L) issues <- c(issues, sprintf("%d NaN", n_nan))
  if (n_inf > 0L) issues <- c(issues, sprintf("%d Inf/-Inf", n_inf))

  if (length(issues) > 0L) {
    stop(
      sprintf(
        "Input contains non-finite values (%s). Clean the series before filtering.",
        paste(issues, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}
