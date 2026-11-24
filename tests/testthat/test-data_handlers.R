# ── Tests: ensure_computable() dispatch and validation ───────────────────────

# ── 1. Subclasses dispatch through inherits() ────────────────────────────────

test_that("zoo::as.zoo(<ts>) (class zooreg) is accepted and matches the ts", {
  skip_if_not_installed("zoo")
  set.seed(1)
  y_ts <- ts(cumsum(rnorm(80)) + 10, start = c(1995, 1), frequency = 4)
  y_zr <- zoo::as.zoo(y_ts)
  expect_s3_class(y_zr, "zooreg")

  expect_no_warning(r_zr <- hp_filter(y_zr))
  r_ts <- hp_filter(y_ts)
  expect_equal(r_zr$trend, r_ts$trend, tolerance = 1e-12)
  expect_equal(r_zr$meta$freq, 4)
  expect_equal(hamilton_filter(y_zr)$meta$h, 8L)
})

test_that("xts built from a ts is handled as xts, not zoo", {
  skip_if_not_installed("xts")
  y_ts <- ts(cumsum(rnorm(40)), start = c(1995, 1), frequency = 4)
  inputs <- ensure_computable(xts::as.xts(y_ts))
  expect_equal(inputs$class, "xts")
  expect_equal(.detect_frequency(inputs), 4)
})

test_that("integer vectors and single-column matrices are accepted as vectors", {
  y <- c(3L, 5L, 4L, 6L, 8L, 7L, 9L, 11L, 10L, 12L)
  r_int <- hp_filter(y, lambda = 100)
  expect_equal(r_int$trend, hp_filter(as.double(y), lambda = 100)$trend)
  r_mat <- hp_filter(matrix(as.double(y), ncol = 1), lambda = 100)
  expect_equal(r_mat$trend, r_int$trend)
})

# ── 2. Multi-column input is rejected, never flattened ───────────────────────

test_that("multi-column inputs error instead of being flattened", {
  m <- matrix(cumsum(rnorm(80)), ncol = 2)
  expect_error(hp_filter(m, lambda = 1600), "`x` has 2 columns; filters are univariate")
  expect_error(hp_filter(ts(m, frequency = 4)), "has 2 columns")

  skip_if_not_installed("xts")
  dates <- seq(as.Date("2000-01-01"), by = "quarter", length.out = 40)
  expect_error(mbh_filter(xts::xts(m, order.by = dates)), "has 2 columns")
  expect_error(hamilton_filter(zoo::zoo(m, order.by = dates)), "has 2 columns")
})

# ── 3. Unsupported classes ───────────────────────────────────────────────────

test_that("unsupported classes give an informative error", {
  df <- data.frame(y = cumsum(rnorm(20)))
  expect_error(hp_filter(df), "Unsupported input class 'data.frame'")
  expect_error(hp_filter(df), "df\\$column")
  expect_error(hp_filter(as.character(1:20)), "Unsupported input class 'character'")
  expect_error(hp_filter(Sys.Date() + 0:19), "Unsupported input class 'Date'")
})
