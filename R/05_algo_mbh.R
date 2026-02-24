# ── MacroBoost Hybrid Filter (MBH) ────────────────────────────────────────────
# The flagship algorithm. Combines P-Splines with Component-wise Gradient
# Boosting and Huber loss to achieve robustness against structural breaks.

#' MacroBoost Hybrid (MBH) Filter
#'
#' Decomposes a time series into trend and cycle using a robust boosting
#' algorithm. Unlike the HP filter, MBH uses the Huber loss function to
#' automatically downweight outliers (like the COVID-19 shock), preventing
#' them from distorting the trend.
#'
#' @param x Numeric vector, `ts`, `xts`, or `zoo` object.
#' @param knots Integer.
#'   Number of interior knots for the P-Spline.
#'   If `NULL` (default), it is calculated as `max(20, floor(n / 2))`.
#'   High knot density is required for the trend to be flexible enough.
#' @param mstop Integer.
#'   Number of boosting iterations (default 500).
#'   Robust loss functions require more iterations to converge than L2 loss.
#' @param d Numeric.
#'   The delta parameter for Huber loss (default 0.01).
#'   Values larger than the shock size (in log-diff) will treat the shock as
#'   normal.  Values smaller than the shock size will treat it as an outlier.
#' @param nu Numeric.
#'   The learning rate (shrinkage) for boosting (default 0.2).
#'
#' @details
#' The model estimated is an additive model:
#' \deqn{y_t = \text{Linear}(t) + \text{Smooth}(t) + \epsilon_t}
#'
#' It is fitted using [mboost::mboost()] with:
#' \itemize{
#'   \item **Base Learners:** A linear time trend ([mboost::bols()]) to capture
#'     the global path, plus a B-spline ([mboost::bbs()]) to capture local
#'     curvature.
#'   \item **Loss Function:** Huber loss ([mboost::Huber()]) with parameter
#'     `d`.  This is the key to robustness.
#' }
#'
#' The default parameters (`knots = n/2`, `mstop = 500`) are calibrated to
#' mimic the flexibility of a standard HP filter while retaining the robustness
#' of the Huber loss.
#'
#' @return A `macrofilter` object with `trend`, `cycle`, `data`, and `meta`
#'   components.  The `meta` list contains `method = "MBH"`, `knots`, `d`,
#'   `mstop`, `nu`, and `compute_time`.
#'
#' @export
#' @importFrom mboost mboost bbs bols Huber boost_control
#' @importFrom stats fitted
#'
#' @examples
#' # Quarterly GDP-like series
#' set.seed(42)
#' y <- ts(cumsum(rnorm(200)), start = c(2000, 1), frequency = 4)
#' result <- mbh_filter(y)
#' print(result)
mbh_filter <- function(x, knots = NULL, mstop = 500L, d = 0.01, nu = 0.2) {

  # 1. Ingest ----------------------------------------------------------------
  inputs <- ensure_computable(x)
  y <- inputs$y
  n <- length(y)

  # 2. Length validation ------------------------------------------------------
  if (n < 6L) {
    stop("Series must have at least 6 observations for MBH filter.", call. = FALSE)
  }

  # 3. Knots heuristic (aggressive) ------------------------------------------
  # High flexibility (~1 knot every 2 obs) because the Huber loss will be the
  # smoothness constraint, not the spline basis.
  if (is.null(knots)) {
    knots <- max(20L, floor(n / 2L))
  }
  # Clamp knots to avoid singular spline bases with very short series
  max_knots <- max(1L, n - 4L)
  knots <- min(as.integer(knots), max_knots)
  mstop <- as.integer(mstop)

  # 4. Timer -----------------------------------------------------------------
  t0 <- proc.time()

  # 5. Base learners ---------------------------------------------------------
  time_idx <- seq_len(n)
  df_boost <- data.frame(y = y, time_idx = time_idx)
  bl_linear <- mboost::bols(time_idx, intercept = TRUE)
  bl_smooth <- mboost::bbs(time_idx, knots = knots, degree = 3, differences = 2)

  # 6. Fit -------------------------------------------------------------------
  fam  <- mboost::Huber(d = d)
  ctrl <- mboost::boost_control(mstop = mstop, nu = nu)

  utils::capture.output(suppressWarnings(
    mod <- mboost::mboost(
      y ~ bl_linear + bl_smooth,
      data = df_boost, family = fam, control = ctrl
    )
  ))

  # 7. Extract ---------------------------------------------------------------
  trend_num <- as.numeric(fitted(mod))
  cycle_num <- y - trend_num

  elapsed <- (proc.time() - t0)[["elapsed"]]

  # 8. Package into macrofilter ----------------------------------------------
  result <- new_macrofilter(
    cycle = cycle_num,
    trend = trend_num,
    data  = y,
    meta  = list(
      method       = "MBH",
      knots        = knots,
      d            = d,
      mstop        = mstop,
      nu           = nu,
      compute_time = elapsed
    )
  )

  # 9. Validate --------------------------------------------------------------
  validate_macrofilter(result)

  result
}
