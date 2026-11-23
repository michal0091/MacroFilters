# ── MacroFilters S3 Class System ──────────────────────────────────────────────
# Low-level constructor, validator, and print method for the `macrofilter` class.

# ── Constructor (low-level) ──────────────────────────────────────────────────

#' Create a macrofilter object (low-level)
#'
#' The single place where filter results are assembled. The cycle is derived
#' here as `data - trend`, so the additive identity holds by construction and
#' needs no runtime check. The result is validated before it is returned.
#'
#' @param trend Numeric vector. The trend component (`NA` over the lead-in
#'   window declared by `meta$na_lead`).
#' @param data  Numeric vector. The original series.
#' @param meta  A named list with at least `method` (character) and `na_lead`
#'   (number of leading positions with no estimate: `h + p - 1` for Hamilton,
#'   `0` otherwise).
#' @param bands `NULL`, or a list with numeric `lower` and `upper` vectors
#'   (the output of `.boot_engine()`, padded to full length).
#' @return An object of class `c("macrofilter", "list")`.
#' @keywords internal
#' @noRd
new_macrofilter <- function(trend, data, meta, bands = NULL) {
  trend <- as.double(trend)
  data  <- as.double(data)
  stopifnot(length(trend) == length(data), is.list(meta),
            is.character(meta$method))

  out <- list(trend = trend, cycle = data - trend, data = data, meta = meta)
  if (!is.null(bands)) {
    out$trend_lower <- as.double(bands$lower)
    out$trend_upper <- as.double(bands$upper)
  }
  class(out) <- c("macrofilter", "list")
  validate_macrofilter(out)
}

# ── Validator ────────────────────────────────────────────────────────────────

#' Validate a macrofilter object
#'
#' Checks what can actually go wrong, not the additive identity (which
#' `new_macrofilter()` guarantees by construction):
#' \itemize{
#'   \item `trend`, `cycle`, `data` (and the bands, if any) have equal length;
#'   \item `trend` is `NA` exactly over the first `meta$na_lead` positions and
#'     finite everywhere else;
#'   \item bands, if present, come as a pair, share that `NA` pattern, and
#'     satisfy `lower <= upper`.
#' }
#' `lower <= trend <= upper` is deliberately not checked: with symmetric
#' normal bands it holds by construction.
#'
#' @param x A `macrofilter` object.
#' @return `x`, invisibly, if valid.
#' @keywords internal
#' @noRd
validate_macrofilter <- function(x) {
  n <- length(x$data)
  lens <- c(trend = length(x$trend), cycle = length(x$cycle))
  if (any(lens != n)) {
    stop(sprintf("Length mismatch: data has %d obs, trend has %d, cycle has %d.",
                 n, lens[["trend"]], lens[["cycle"]]), call. = FALSE)
  }

  na_lead <- x$meta$na_lead
  if (length(na_lead) != 1L || is.na(na_lead) || na_lead < 0L || na_lead >= n) {
    stop("`meta$na_lead` must be a single integer in [0, n).", call. = FALSE)
  }
  lead <- seq_len(na_lead)
  body <- setdiff(seq_len(n), lead)

  .check_window <- function(v, what) {
    if (!all(is.na(v[lead]))) {
      stop(sprintf("`%s` must be NA over the first %d positions (meta$na_lead).",
                   what, na_lead), call. = FALSE)
    }
    if (!all(is.finite(v[body]))) {
      stop(sprintf("`%s` has %d non-finite values outside the NA lead-in.",
                   what, sum(!is.finite(v[body]))), call. = FALSE)
    }
  }
  .check_window(x$trend, "trend")

  has_lo <- !is.null(x$trend_lower)
  has_up <- !is.null(x$trend_upper)
  if (has_lo != has_up) {
    stop("`trend_lower` and `trend_upper` must be both present or both absent.",
         call. = FALSE)
  }
  if (has_lo) {
    if (length(x$trend_lower) != n || length(x$trend_upper) != n) {
      stop("Bootstrap bands must have the same length as the data.", call. = FALSE)
    }
    .check_window(x$trend_lower, "trend_lower")
    .check_window(x$trend_upper, "trend_upper")
    if (any(x$trend_lower[body] > x$trend_upper[body])) {
      stop("Bootstrap bands are inverted: `trend_lower > trend_upper`.",
           call. = FALSE)
    }
  }

  invisible(x)
}

# ── Print method ─────────────────────────────────────────────────────────────

#' @export
#' @importFrom stats sd
print.macrofilter <- function(x, ...) {
  m <- x$meta
  n <- length(x$data)

  cat(sprintf("-- MacroFilter [%s] --\n", m$method))
  cat(sprintf("   Observations : %d\n", n))


  # Show relevant scalar parameters per method
  # (exclude bookkeeping + temporal-identity fields, not user-facing params)
  params <- setdiff(names(m), c("method", "compute_time", "na_lead",
                                "ts_class", "tsp", "idx"))
  params <- Filter(function(p) length(m[[p]]) == 1L, params)
  if (length(params) > 0L) {
    vals <- vapply(params, function(p) format(m[[p]], digits = 4), character(1))
    cat(sprintf("   Parameters   : %s\n", paste(params, vals, sep = " = ", collapse = ", ")))
  }

  # Cycle summary (skip NAs)
  cy <- x$cycle[!is.na(x$cycle)]
  if (length(cy) > 0L) {
    cat(sprintf(
      "   Cycle range  : [%s, %s]  sd = %s\n",
      format(min(cy), digits = 4),
      format(max(cy), digits = 4),
      format(sd(cy),  digits = 4)
    ))
  }

  if (!is.null(m$compute_time)) {
    cat(sprintf("   Compute time : %.3f s\n", as.numeric(m$compute_time)))
  }

  invisible(x)
}