# ── autoplot method for macrofilter objects ───────────────────────────────────

# Re-export the ggplot2 generic so users can call autoplot() on a macrofilter
# without attaching ggplot2 first.
#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' Plot a macrofilter decomposition with ggplot2
#'
#' Visualises any `macrofilter` object: the observed series (attenuated) over
#' the estimated trend (emphasised). When the object carries bootstrap bands
#' (`$trend_lower` / `$trend_upper`), a 95% confidence ribbon is drawn behind
#' the lines automatically.
#'
#' The x axis uses the time stamps of the original input: decimal time for a
#' `ts`, `yearqtr` or `yearmon` index, dates for a `Date` or `POSIXct` index
#' (`xts`/`zoo`), and the positions `1..n` for a plain numeric vector.
#'
#' @param object A `macrofilter` object.
#' @param ... Currently ignored; present for S3 generic compatibility.
#'
#' @return A `ggplot` object.
#'
#' @export
#' @method autoplot macrofilter
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon scale_color_manual
#' @importFrom ggplot2 labs theme_minimal theme element_blank
#' @importFrom data.table data.table
#'
#' @examples
#' \donttest{
#' y <- ts(cumsum(rnorm(120)), start = c(2000, 1), frequency = 4)
#' autoplot(mbh_filter(y, mstop = 100L, boot_iter = 50L))
#' }
autoplot.macrofilter <- function(object, ...) {

  # 4.1 Time axis + wide data ------------------------------------------------
  if (!inherits(object, "macrofilter")) {
    stop("`object` must be a 'macrofilter' object.", call. = FALSE)
  }

  time_vec <- .time_axis(object)

  DT <- data.table::data.table(
    time  = time_vec,
    Data  = as.numeric(object$data),
    Trend = as.numeric(object$trend)
  )

  has_ci <- !is.null(object$trend_lower) && !is.null(object$trend_upper)
  if (has_ci) {
    DT[, lower := as.numeric(object$trend_lower)]
    DT[, upper := as.numeric(object$trend_upper)]
  }

  # 4.2 Layered ggplot (ribbon at the back) ----------------------------------
  p <- ggplot2::ggplot(DT, ggplot2::aes(x = time))

  if (has_ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "#B0C4DE", alpha = 0.4, na.rm = TRUE
    )
  }

  p <- p +
    ggplot2::geom_line(
      ggplot2::aes(y = Data, color = "Observed Data"),
      alpha = 0.5, linewidth = 0.6, na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = Trend, color = "Structural Trend"),
      linewidth = 0.9, na.rm = TRUE
    )

  # 4.3 Aesthetics + metadata ------------------------------------------------
  m_name <- if (!is.null(object$meta$method)) object$meta$method else "Filter"

  p +
    ggplot2::scale_color_manual(
      values = c("Observed Data" = "gray40", "Structural Trend" = "#0047AB")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title    = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title    = paste(m_name, "Trend Extraction"),
      subtitle = if (has_ci) "With 95% Bootstrap Confidence Intervals" else NULL,
      x = "Time",
      y = "Value"
    )
}


#' Build the x axis of a macrofilter plot from its stored time identity
#'
#' Uses `meta$tsp` for `ts` input, building the sequence with `length.out` so
#' floating-point steps can never add or drop a point. Otherwise uses
#' `meta$idx` (`xts`/`zoo` input): `Date` and `POSIXct` stay as they are so
#' ggplot2 draws a date axis, and `yearqtr`/`yearmon` become decimal time,
#' which ggplot2 can scale without zoo's scale functions. Anything else, or an
#' index whose length does not match the series, falls back to `1..n`.
#'
#' @param object A `macrofilter` object.
#' @return A vector of length `length(object$trend)`.
#' @keywords internal
#' @noRd
.time_axis <- function(object) {
  n   <- length(object$trend)
  m   <- object$meta
  idx <- m$idx

  if (!is.null(m$tsp)) {
    return(seq(from = m$tsp[1L], by = 1 / m$tsp[3L], length.out = n))
  }
  if (!is.null(idx) && length(idx) == n) {
    # Rebuild from the numeric values: an xts index carries `tclass`/`tzone`
    # attributes that should not leak into the plot data.
    if (inherits(idx, "Date")) return(structure(as.numeric(idx), class = "Date"))
    if (inherits(idx, "POSIXct")) {
      tz <- attr(idx, "tzone")
      return(as.POSIXct(as.numeric(idx), origin = "1970-01-01",
                        tz = if (is.null(tz)) "" else tz[1L]))
    }
    if (inherits(idx, c("yearqtr", "yearmon")) || is.numeric(idx)) {
      return(as.numeric(idx))
    }
  }
  seq_len(n)
}
