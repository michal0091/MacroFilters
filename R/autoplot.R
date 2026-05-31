# ── autoplot method for macrofilter objects ───────────────────────────────────

# Re-export the ggplot2 generic so users can call autoplot() on a macrofilter
# without attaching ggplot2 first.
#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' Plot a macrofilter decomposition with ggplot2
#'
#' Produces a `ggplot` of the original series (attenuated) overlaid with the
#' estimated trend (emphasised). When the object carries bootstrap bands
#' (`$trend_lower` / `$trend_upper`, present when `mbh_filter()` was called
#' with `boot_iter > 0`), a confidence ribbon is drawn automatically beneath
#' the trend line.
#'
#' @param object A `macrofilter` object.
#' @param ... Currently ignored; present for S3 generic compatibility.
#'
#' @return A `ggplot` object.
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs theme_minimal
#'   scale_colour_manual scale_linewidth_manual scale_alpha_manual
#' @importFrom data.table data.table melt
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' y <- ts(cumsum(rnorm(120)), start = c(2000, 1), frequency = 4)
#' fit <- mbh_filter(y, mstop = 100L, boot_iter = 50L)
#' autoplot(fit)
#' }
autoplot.macrofilter <- function(object, ...) {
  has_ci <- !is.null(object$trend_lower) && !is.null(object$trend_upper)
  tax      <- .mf_time_axis(object)   # real dates / decimal time / index
  time_idx <- tax$time

  # Long-format data for the line layers (observed + trend) via melt
  wide <- data.table::data.table(
    time     = time_idx,
    Observed = as.numeric(object$data),
    Trend    = as.numeric(object$trend)
  )
  long <- data.table::melt(
    wide,
    id.vars       = "time",
    measure.vars  = c("Observed", "Trend"),
    variable.name = "series",
    value.name    = "value"
  )

  p <- ggplot2::ggplot()

  # Confidence ribbon (drawn first so it sits behind the lines)
  if (has_ci) {
    ribbon <- data.table::data.table(
      time  = time_idx,
      lower = as.numeric(object$trend_lower),
      upper = as.numeric(object$trend_upper)
    )
    p <- p + ggplot2::geom_ribbon(
      data    = ribbon,
      mapping = ggplot2::aes(x = time, ymin = lower, ymax = upper),
      fill    = "#2C7FB8",
      alpha   = 0.20
    )
  }

  p <- p +
    ggplot2::geom_line(
      data    = long,
      mapping = ggplot2::aes(
        x = time, y = value,
        colour    = series,
        linewidth = series,
        alpha     = series
      )
    ) +
    ggplot2::scale_colour_manual(
      values = c(Observed = "grey55", Trend = "#2C7FB8")
    ) +
    ggplot2::scale_linewidth_manual(
      values = c(Observed = 0.4, Trend = 1.1)
    ) +
    ggplot2::scale_alpha_manual(
      values = c(Observed = 0.6, Trend = 1.0)
    ) +
    ggplot2::labs(
      title    = sprintf("MacroFilter decomposition [%s]", object$meta$method),
      x        = tax$xlab,
      y        = "Value",
      colour   = NULL, linewidth = NULL, alpha = NULL
    ) +
    ggplot2::theme_minimal()

  p
}

# ── Internal: reconstruct the x-axis from the object's temporal identity ──────

#' Build the plotting time axis for a macrofilter
#'
#' Uses the temporal metadata stored by the filters (`meta$tsp` for `ts`
#' input, `meta$idx` for `xts`/`zoo`) to map observations back to real time.
#' Falls back to a positional index for plain numeric input.
#'
#' @param object A `macrofilter` object.
#' @return A list with `time` (axis values) and `xlab` (axis label).
#' @keywords internal
#' @noRd
.mf_time_axis <- function(object) {
  n <- length(object$trend)
  m <- object$meta

  if (!is.null(m$tsp)) {
    # ts input: evenly spaced decimal time from start to end
    list(time = seq(m$tsp[1L], m$tsp[2L], length.out = n), xlab = "Time")
  } else if (!is.null(m$idx) && length(m$idx) == n) {
    # xts / zoo input: real date index
    list(time = m$idx, xlab = "Date")
  } else {
    # plain numeric: positional index
    list(time = seq_len(n), xlab = "Time index")
  }
}
