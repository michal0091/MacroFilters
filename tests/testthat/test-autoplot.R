# ── Tests: autoplot.macrofilter() ─────────────────────────────────────────────

geoms_of <- function(p) vapply(p$layers, function(l) class(l$geom)[1L],
                               character(1L))

# ── 1. Returns a ggplot ──────────────────────────────────────────────────────

test_that("autoplot returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  y <- ts(cumsum(rnorm(80)) + 10, start = c(2000, 1), frequency = 4)
  p <- ggplot2::autoplot(hp_filter(y))
  expect_s3_class(p, "ggplot")
})

# ── 2. Input validation ──────────────────────────────────────────────────────

test_that("autoplot rejects non-macrofilter input", {
  skip_if_not_installed("ggplot2")
  # call the method directly: dispatch would never route a plain list here
  expect_error(autoplot.macrofilter(list(a = 1)), "macrofilter")
})

# ── 3. Ribbon present iff bands exist ─────────────────────────────────────────

test_that("ribbon layer appears only when bands are present", {
  skip_if_not_installed("ggplot2")
  y <- ts(cumsum(rnorm(80)) + 10, start = c(2000, 1), frequency = 4)

  p_plain <- ggplot2::autoplot(hp_filter(y))
  expect_false("GeomRibbon" %in% geoms_of(p_plain))

  p_band <- ggplot2::autoplot(hp_filter(y, boot_iter = 20L))
  expect_true("GeomRibbon" %in% geoms_of(p_band))
})

# ── 4. X-axis reconstruction ─────────────────────────────────────────────────

test_that("time axis is reconstructed from tsp when available", {
  skip_if_not_installed("ggplot2")

  # ts input: decimal time rebuilt from meta$tsp, starting at 2000
  y_ts <- ts(cumsum(rnorm(80)), start = c(2000, 1), frequency = 4)
  p_ts <- ggplot2::autoplot(hp_filter(y_ts))
  expect_equal(p_ts$data$time[1L], 2000)
  expect_equal(p_ts$labels$x, "Time")          # label hardcoded by design

  # plain numeric: positional index 1..n
  v <- cumsum(rnorm(80))
  p_v <- ggplot2::autoplot(suppressWarnings(hp_filter(v)))
  expect_equal(p_v$data$time, seq_len(80))
})

test_that("time axis uses the xts / zoo index", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("xts")
  dates <- seq(as.Date("1995-01-01"), by = "quarter", length.out = 40)
  y     <- cumsum(rnorm(40)) + 10

  p_x <- ggplot2::autoplot(hp_filter(xts::xts(y, order.by = dates)))
  expect_s3_class(p_x$data$time, "Date")
  expect_equal(p_x$data$time[1L], as.Date("1995-01-01"))
  expect_equal(p_x$data$time, dates)

  # POSIXct index keeps its time zone
  stamps <- seq(as.POSIXct("1995-01-01", tz = "UTC"), by = "month", length.out = 40)
  p_p <- ggplot2::autoplot(hp_filter(xts::xts(y, order.by = stamps)))
  expect_equal(p_p$data$time, stamps)

  # zoo::as.zoo(<ts>) carries a yearqtr index: drawn as decimal time
  y_zr <- zoo::as.zoo(ts(y, start = c(1995, 1), frequency = 4))
  p_z  <- ggplot2::autoplot(hp_filter(y_zr))
  expect_equal(p_z$data$time, 1995 + (0:39) / 4)
})

test_that("tsp axis has exactly n points for awkward frequencies", {
  skip_if_not_installed("ggplot2")
  for (f in c(4, 7, 12, 52)) {
    y <- ts(cumsum(rnorm(137)), start = c(2001, 3), frequency = f)
    p <- ggplot2::autoplot(hp_filter(y, lambda = 1600))
    expect_length(p$data$time, 137)
    expect_equal(p$data$time, as.numeric(stats::time(y)), tolerance = 1e-10,
                 info = paste("freq", f))
  }
})

test_that("an index of the wrong length falls back to 1..n", {
  obj <- hp_filter(cumsum(rnorm(20)), lambda = 100)
  obj$meta$idx <- as.Date("2000-01-01") + 0:4
  expect_equal(.time_axis(obj), seq_len(20))
})

# ── 5. Subtitle only with confidence bands ───────────────────────────────────

test_that("CI subtitle appears only when bands are present", {
  skip_if_not_installed("ggplot2")
  y <- ts(cumsum(rnorm(80)) + 10, start = c(2000, 1), frequency = 4)

  expect_null(ggplot2::autoplot(hp_filter(y))$labels$subtitle)
  expect_match(
    ggplot2::autoplot(hp_filter(y, boot_iter = 20L))$labels$subtitle,
    "Bootstrap Confidence Intervals"
  )
})
