# ── Tests: hp_filter() ────────────────────────────────────────────────────────

# ── 1. Analytical properties ─────────────────────────────────────────────────

test_that("HP: linear input is a fixed point", {
  y <- 3 + 0.5 * (1:100)
  r <- hp_filter(y, lambda = 1600)
  expect_equal(r$trend, y, tolerance = 1e-10)
  expect_equal(r$cycle, rep(0, 100), tolerance = 1e-10)
})

test_that("HP: lambda -> 0 returns the data; lambda -> Inf returns the OLS line", {
  set.seed(1)
  y <- cumsum(rnorm(100))
  t <- 1:100
  expect_equal(hp_filter(y, lambda = 1e-12)$trend, y, tolerance = 1e-6)
  ols <- stats::lm.fit(cbind(1, t), y)$fitted.values
  # By linearity and the linear fixed point, HP(y) - ols = HP(y - ols). The
  # system has condition number ~1.6e13 at lambda = 1e12, so comparing HP(y)
  # with ols directly measures rounding error, which depends on the platform
  # (1.1e-4 on x86_64, 1.3e-3 on macOS arm64). Smoothing the residual keeps
  # the solution small, and with it the rounding error; the exact gap is
  # ~8e-8.
  expect_lt(max(abs(hp_filter(y - ols, lambda = 1e12)$trend)), 1e-6)
})

# ── 2. Invariance: same result regardless of input class ─────────────────────

test_that("HP filter produces identical results across input classes", {
  set.seed(42)
  y_num <- cumsum(rnorm(100))
  y_ts  <- ts(y_num, start = c(2000, 1), frequency = 4)

  res_num <- hp_filter(y_num, lambda = 1600)
  res_ts  <- hp_filter(y_ts, lambda = 1600)

  expect_equal(res_num$trend, res_ts$trend, tolerance = 1e-12)
  expect_equal(res_num$cycle, res_ts$cycle, tolerance = 1e-12)
})

test_that("HP filter produces identical results with xts input", {
  skip_if_not_installed("xts")
  set.seed(42)
  y_num <- cumsum(rnorm(100))
  dates <- seq(as.Date("2000-01-01"), by = "quarter", length.out = 100)
  y_xts <- xts::xts(y_num, order.by = dates)

  res_num <- hp_filter(y_num, lambda = 1600)
  res_xts <- hp_filter(y_xts, lambda = 1600)

  expect_equal(res_num$trend, res_xts$trend, tolerance = 1e-12)
  expect_equal(res_num$cycle, res_xts$cycle, tolerance = 1e-12)
})

test_that("HP filter produces identical results with zoo input", {
  skip_if_not_installed("zoo")
  set.seed(42)
  y_num <- cumsum(rnorm(100))
  dates <- seq(as.Date("2000-01-01"), by = "quarter", length.out = 100)
  y_zoo <- zoo::zoo(y_num, order.by = dates)

  res_num <- hp_filter(y_num, lambda = 1600)
  res_zoo <- hp_filter(y_zoo, lambda = 1600)

  expect_equal(res_num$trend, res_zoo$trend, tolerance = 1e-12)
  expect_equal(res_num$cycle, res_zoo$cycle, tolerance = 1e-12)
})

# ── 3. Stress tests ──────────────────────────────────────────────────────────

test_that("HP filter works with minimum-length series (n = 3, 4, 5)", {
  for (n in 3:5) {
    y <- seq_len(n) * 1.0
    result <- hp_filter(y, lambda = 1600)
    expect_equal(result$trend + result$cycle, result$data, tolerance = 1e-9)
  }
})

test_that("HP filter errors when n < 3", {
  expect_error(hp_filter(c(1, 2), lambda = 1600), "at least 3")
  expect_error(hp_filter(1, lambda = 1600), "at least 3")
})

test_that("HP filter handles long series (n = 10000)", {
  set.seed(99)
  y <- cumsum(rnorm(10000))
  result <- hp_filter(y, lambda = 1600)
  expect_equal(result$trend + result$cycle, result$data, tolerance = 1e-9)
})

test_that("Auto-lambda from ts frequency = 4 gives lambda = 1600", {
  y_ts <- ts(1:50 * 1.0, frequency = 4)
  result <- hp_filter(y_ts)
  expect_equal(result$meta$lambda, 6.25 * 4^4)
})

test_that("Auto-lambda from ts frequency = 12 gives lambda = 129600", {
  y_ts <- ts(1:50 * 1.0, frequency = 12)
  result <- hp_filter(y_ts)
  expect_equal(result$meta$lambda, 6.25 * 12^4)
})

test_that("freq parameter overrides default when lambda is NULL", {
  y <- 1:50 * 1.0
  result <- suppressWarnings(hp_filter(y, freq = 12))
  expect_equal(result$meta$lambda, 6.25 * 12^4)
})

test_that("Warning emitted when frequency cannot be determined", {
  y <- cumsum(rnorm(50))
  expect_warning(hp_filter(y), "Cannot determine series frequency")
})

# ── 4. Output format ─────────────────────────────────────────────────────────

test_that("HP filter returns macrofilter object with correct meta", {
  set.seed(7)
  y <- cumsum(rnorm(100))
  result <- hp_filter(y, lambda = 1600)

  expect_s3_class(result, "macrofilter")
  expect_equal(result$meta$method, "HP")
  expect_equal(result$meta$lambda, 1600)
  expect_true(!is.null(result$meta$compute_time))
  expect_true(is.numeric(result$meta$compute_time))
})

test_that("Explicit lambda is stored exactly in meta", {
  y <- 1:20 * 1.0
  result <- hp_filter(y, lambda = 6.25)
  expect_equal(result$meta$lambda, 6.25)
})
