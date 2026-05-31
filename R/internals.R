# ── Internal utilities (not exported) ─────────────────────────────────────────

#' Fast HP cycle for internal calibration
#'
#' Solves (I + lambda * D'D) t = y via sparse Cholesky and returns only the
#' cycle (y - trend). No input validation, no class wrapping, no timing.
#' `stats::frequency(x)` drives the Ravn-Uhlig lambda so the threshold is
#' always on the correct residual scale.
#'
#' @param x Numeric vector or `ts` object.
#' @return Numeric vector of cyclical residuals, same length as `x`.
#' @keywords internal
#' @noRd
.hp_fast <- function(x) {
  n <- length(x)
  if (n < 3L) return(rep(0, n))          # bandSparse needs n-2 >= 1 rows

  freq   <- stats::frequency(x)          # 1 for numeric/annual, 4/12 for ts
  y      <- as.double(x)
  lambda <- 6.25 * freq^4                # Ravn-Uhlig rule

  D <- Matrix::bandSparse(
    n - 2L, n,
    k         = c(0L, 1L, 2L),
    diagonals = list(rep(1, n - 2L), rep(-2, n - 2L), rep(1, n - 2L))
  )
  Q <- Matrix::Diagonal(n) + lambda * Matrix::crossprod(D)
  y - as.numeric(Matrix::solve(Q, y))
}


#' Stripped MBH trend for bootstrap iterations
#'
#' Runs mboost with a reduced iteration budget so that applying it hundreds
#' of times in `.boot_engine` is tractable. Uses the same spline geometry as
#' the base fit (same `knots`, `df`, `boundary.knots`) so bootstrap trends
#' live in the same function space.
#'
#' @param y      Numeric vector (synthetic bootstrap draw).
#' @param d_val  Calibrated Huber threshold from the base fit.
#' @param knots  Clamped knot count from the base fit.
#' @param mstop_boot Reduced iteration budget (typically `max(50, mstop %/% 5)`).
#' @param nu     Learning rate from the base fit.
#' @param df     P-spline degrees of freedom from the base fit.
#' @param boundary.knots Boundary knots from the base fit (or `NULL`).
#' @return Numeric trend vector, same length as `y`.
#' @keywords internal
#' @noRd
.mbh_fast_trend <- function(y, d_val, knots, mstop_boot, nu, df,
                             boundary.knots) {
  n        <- length(y)
  time_idx <- seq_len(n)

  bl_linear <- mboost::bols(time_idx, intercept = TRUE)
  bl_smooth <- mboost::bbs(time_idx, knots = knots, degree = 3L,
                            differences = 2L, df = df,
                            boundary.knots = boundary.knots)

  # Native silencing via trace = FALSE, avoids the per-iteration text
  # connection that capture.output() opens and tears down on every call.
  ctrl <- mboost::boost_control(mstop = mstop_boot, nu = nu, trace = FALSE)

  suppressWarnings({
    mod <- mboost::mboost(y ~ bl_linear + bl_smooth,
                          data = data.frame(y = y, time_idx = time_idx),
                          family = mboost::Huber(d = d_val),
                          control = ctrl)
  })
  as.numeric(stats::fitted(mod))
}


#' Universal block-bootstrap engine (Circular Block Bootstrap)
#'
#' Builds a long `data.table` of synthetic series (one per iteration) by
#' circular-block resampling of the pseudo-residuals, applies `filter_func`
#' to each block via `data.table`'s `by` grouping, and returns 2.5%/97.5%
#' quantile bands by time point.
#'
#' Uses the Circular Block Bootstrap (Politis & Romano, 1992): block starts
#' are drawn uniformly from `1:n` and blocks wrap around the series end with
#' modular arithmetic. Every observation therefore appears in exactly
#' `block_size` candidate blocks, eliminating the end-point under-representation
#' of the plain moving-block bootstrap — critical because the trend's
#' end-of-sample uncertainty is precisely what we want to quantify.
#'
#' Index generation uses pure vector algebra:
#'   1. Sample `blk_n` block starts per iteration (all at once, from `1:n`).
#'   2. Expand each start with `rep()` + offset arithmetic.
#'   3. Wrap overshoot with `(idx - 1) %% n + 1` (circular, preserves variance).
#'   4. Flatten into the long `data.table` via positional arithmetic.
#'
#' @param filter_func Function `f(y_numeric) -> numeric trend`, same length.
#' @param trend_base  Numeric trend from the base fit (length `n`).
#' @param cycle_base  Numeric cycle (pseudo-residuals) from the base fit.
#' @param boot_iter   Positive integer number of iterations.
#' @param block_size  Block length (already resolved and validated by caller).
#' @return Named list with `lower` and `upper` numeric vectors (length `n`).
#' @keywords internal
#' @noRd
.boot_engine <- function(filter_func, trend_base, cycle_base,
                          boot_iter, block_size) {
  n       <- length(trend_base)
  blk_n   <- ceiling(n / block_size)
  raw_len <- blk_n * block_size            # >= n

  # Step 1-3: All indices for all iterations
  all_starts <- sample.int(n,
                            size    = boot_iter * blk_n,
                            replace = TRUE)

  offsets <- 0L:(block_size - 1L)
  raw_idx <- rep(all_starts, each = block_size) +
             rep(offsets, times  = boot_iter * blk_n)
  raw_idx <- (raw_idx - 1L) %% n + 1L      # circular wrap (no variance collapse)

  # Step 4: Map (iter, time) pairs → positions in raw_idx
  iter_id  <- rep(seq_len(boot_iter), each = n)
  time_pos <- rep(seq_len(n), times = boot_iter)
  flat_pos <- (iter_id - 1L) * raw_len + time_pos

  # Step 5: Build long data.table
  DT <- data.table::data.table(
    iter     = iter_id,
    time_idx = time_pos,
    y_boot   = trend_base[time_pos] + cycle_base[raw_idx[flat_pos]]
  )

  # Step 6: Apply filter and compute quantiles by reference
  DT[, boot_trend := filter_func(y_boot), by = iter]

  q_DT <- DT[, .(
    lower = stats::quantile(boot_trend, 0.025, names = FALSE),
    upper = stats::quantile(boot_trend, 0.975, names = FALSE)
  ), by = time_idx]
  data.table::setorder(q_DT, time_idx)

  list(lower = q_DT$lower, upper = q_DT$upper)
}
