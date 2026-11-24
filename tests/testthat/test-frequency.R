# ── Tests: frequency resolution and bootstrap block size ─────────────────────

q_dates <- function(n) seq(as.Date("1995-01-01"), by = "quarter", length.out = n)

# ── 1. Index inference ───────────────────────────────────────────────────────

test_that(".infer_freq_from_index reads calendar indexes", {
  skip_if_not_installed("zoo")
  expect_equal(.infer_freq_from_index(q_dates(12)), 4)
  expect_equal(.infer_freq_from_index(
    seq(as.Date("1995-01-01"), by = "month", length.out = 12)), 12)
  expect_equal(.infer_freq_from_index(
    seq(as.Date("1995-01-01"), by = "year", length.out = 12)), 1)
  expect_equal(.infer_freq_from_index(
    seq(as.Date("1995-01-02"), by = "week", length.out = 12)), 52)
  expect_equal(.infer_freq_from_index(
    seq(as.POSIXct("1995-01-01", tz = "UTC"), by = "month", length.out = 12)), 12)
  expect_equal(.infer_freq_from_index(zoo::as.yearqtr(1995 + (0:11) / 4)), 4)
  expect_equal(.infer_freq_from_index(zoo::as.yearmon(1995 + (0:11) / 12)), 12)
  expect_equal(.infer_freq_from_index(1995 + (0:11) / 4), 4)
})

test_that(".infer_freq_from_index returns NULL without calendar information", {
  expect_null(.infer_freq_from_index(NULL))
  expect_null(.infer_freq_from_index(q_dates(2)))
  expect_null(.infer_freq_from_index(1:20))                 # zoo() default index
  expect_null(.infer_freq_from_index(
    seq(as.Date("1995-01-01"), by = "day", length.out = 20)))
})

# ── 2. Precedence and validation ─────────────────────────────────────────────

test_that("a ts attribute takes precedence over a conflicting freq, with a warning", {
  y <- ts(cumsum(rnorm(40)), start = c(2000, 1), frequency = 4)
  expect_warning(r <- hp_filter(y, freq = 12),
                 "`freq = 12` ignored: `x` has frequency 4 from its ts attribute")
  expect_equal(r$meta$lambda, 1600)
  expect_warning(hamilton_filter(y, freq = 12), "ignored")
})

test_that("a calendar index takes precedence over a conflicting freq, with a warning", {
  skip_if_not_installed("xts")
  yx <- xts::xts(cumsum(rnorm(40)), order.by = q_dates(40))
  expect_warning(r <- hp_filter(yx, freq = 12), "from its time index")
  expect_equal(r$meta$freq, 4)
})

test_that("a freq that agrees with x does not warn", {
  y <- ts(cumsum(rnorm(40)), start = c(2000, 1), frequency = 4)
  expect_no_warning(hp_filter(y, freq = 4))
  expect_no_warning(suppressMessages(mbh_filter(y, freq = 4, mstop = 20)))
})

test_that("freq must be a single positive number", {
  y <- cumsum(rnorm(40))
  expect_error(hp_filter(y, freq = 0), "single positive number")
  expect_error(hp_filter(y, freq = c(4, 12)), "single positive number")
})

# ── 3. meta$freq ─────────────────────────────────────────────────────────────

test_that("meta$freq records the frequency read from x or freq", {
  set.seed(1)
  y <- cumsum(rnorm(60))
  expect_equal(hp_filter(ts(y, frequency = 4))$meta$freq, 4)
  expect_equal(hp_filter(y, lambda = 1600)$meta$freq, NA_real_)
  expect_equal(hp_filter(y, freq = 12)$meta$freq, 12)
  expect_equal(hamilton_filter(y, freq = 12)$meta$freq, 12)
  expect_equal(bhp_filter(y, freq = 4, stopping = "fixed", iter_max = 2)$meta$freq, 4)
  expect_equal(mbh_filter(y, d = 0.5, mstop = 20)$meta$freq, NA_real_)
})

# ── 4. No new warnings when the frequency is not needed ──────────────────────

test_that("explicit lambda / h on a plain vector does not warn", {
  set.seed(2)
  y <- cumsum(rnorm(60))
  expect_no_warning(hp_filter(y, lambda = 1600, boot_iter = 5))
  expect_no_warning(bhp_filter(y, lambda = 1600, stopping = "fixed",
                               iter_max = 2, boot_iter = 5))
  expect_no_warning(hamilton_filter(y, h = 8, boot_iter = 5))
  expect_no_warning(suppressMessages(mbh_filter(y, mstop = 20, boot_iter = 2)))
})

test_that("hamilton_filter(freq = ) sets the default horizon silently", {
  set.seed(3)
  y <- cumsum(rnorm(120))
  expect_no_warning(r <- hamilton_filter(y, freq = 12))
  expect_equal(r$meta$h, 24L)
})

test_that("mbh_filter(freq = ) sets the calibration lambda", {
  set.seed(4)
  y <- cumsum(rnorm(80))
  d_ts  <- suppressMessages(mbh_filter(ts(y, frequency = 4), mstop = 20))$meta$d
  d_num <- suppressMessages(mbh_filter(y, freq = 4, mstop = 20))$meta$d
  expect_equal(d_num, d_ts)
})

# ── 5. xts / zoo now match ts ────────────────────────────────────────────────

test_that("quarterly xts gets the same defaults as ts(freq = 4)", {
  skip_if_not_installed("xts")
  set.seed(5)
  y   <- cumsum(rnorm(80)) + 10
  yts <- ts(y, start = c(1995, 1), frequency = 4)
  yx  <- xts::xts(y, order.by = q_dates(80))

  expect_no_warning(hx <- hp_filter(yx))
  expect_equal(hx$meta$lambda, 1600)
  expect_no_warning(gx <- hamilton_filter(yx))
  expect_equal(gx$meta$h, 8L)
  expect_equal(bhp_filter(yx)$meta$lambda, 1600)
  expect_equal(suppressMessages(mbh_filter(yx, mstop = 20))$meta$d,
               suppressMessages(mbh_filter(yts, mstop = 20))$meta$d)
})

test_that("zoo with a yearqtr index gets the same defaults as ts(freq = 4)", {
  skip_if_not_installed("zoo")
  set.seed(6)
  y  <- cumsum(rnorm(80)) + 10
  yz <- zoo::zoo(y, order.by = zoo::as.yearqtr(1995 + (0:79) / 4))
  expect_no_warning(r <- hp_filter(yz))
  expect_equal(r$meta$lambda, 1600)
  expect_equal(r$meta$freq, 4)
})

# ── 6. Block size (task 0.3) ─────────────────────────────────────────────────

test_that(".resolve_block_size uses 2 * freq, capped at n_eff / 3", {
  expect_equal(.resolve_block_size(4, "auto", 300), 8L)
  expect_equal(.resolve_block_size(1, "auto", 300), 2L)
  expect_equal(.resolve_block_size(12, "auto", 30), 10L)
  expect_equal(.resolve_block_size(4, 5, 300), 5L)
  expect_equal(.resolve_block_size(4, "auto", 2), 1L)
})

test_that("quarterly xts and ts(freq = 4) give identical bootstrap bands", {
  skip_if_not_installed("xts")
  set.seed(7)
  y   <- cumsum(rnorm(80)) + 10
  yts <- ts(y, start = c(1995, 1), frequency = 4)
  yx  <- xts::xts(y, order.by = q_dates(80))

  fits <- list(
    hp  = function(v) hp_filter(v, boot_iter = 20),
    ham = function(v) hamilton_filter(v, boot_iter = 20),
    mbh = function(v) suppressMessages(mbh_filter(v, mstop = 20, boot_iter = 3))
  )
  for (nm in names(fits)) {
    set.seed(99); a <- fits[[nm]](yts)
    set.seed(99); b <- fits[[nm]](yx)
    expect_equal(b$trend_upper, a$trend_upper, tolerance = 1e-12, info = nm)
    expect_equal(b$trend_lower, a$trend_lower, tolerance = 1e-12, info = nm)
  }
})
