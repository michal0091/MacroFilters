# ── Tests: new_macrofilter() and validate_macrofilter() ──────────────────────

base_meta <- function(na_lead = 0L) list(method = "Test", na_lead = na_lead)

# ── 1. Constructor ───────────────────────────────────────────────────────────

test_that("new_macrofilter derives the cycle and sets the class", {
  y  <- c(1, 3, 2, 5, 4, 6)
  tr <- c(1.5, 2, 2.5, 3.5, 4.5, 5.5)
  obj <- new_macrofilter(tr, y, base_meta())
  expect_s3_class(obj, "macrofilter")
  expect_identical(class(obj), c("macrofilter", "list"))
  expect_identical(obj$cycle, y - tr)
  expect_named(obj, c("trend", "cycle", "data", "meta"))
  expect_null(obj$trend_lower)
})

test_that("new_macrofilter attaches bands after the core fields", {
  y <- 1:6 + 0.5; tr <- y - 0.1
  obj <- new_macrofilter(tr, y, base_meta(),
                         bands = list(lower = tr - 1, upper = tr + 1))
  expect_named(obj, c("trend", "cycle", "data", "meta",
                      "trend_lower", "trend_upper"))
})

test_that("new_macrofilter rejects malformed arguments", {
  expect_error(new_macrofilter(1:3, 1:4, base_meta()))
  expect_error(new_macrofilter(1:3, 1:3, list(na_lead = 0L)))   # no method
})

# ── 2. Validator: each check can fail ────────────────────────────────────────

test_that("validator catches a length mismatch", {
  obj <- new_macrofilter(1:5 + 0, 1:5 + 0, base_meta())
  obj$cycle <- obj$cycle[-1L]
  expect_error(validate_macrofilter(obj), "Length mismatch")
})

test_that("validator checks na_lead", {
  expect_error(new_macrofilter(1:5 + 0, 1:5 + 0, list(method = "Test")),
               "na_lead")
  expect_error(new_macrofilter(1:5 + 0, 1:5 + 0, base_meta(5L)), "na_lead")
})

test_that("validator requires NA exactly over the lead-in", {
  y <- as.double(1:8)
  expect_error(new_macrofilter(y, y, base_meta(2L)), "must be NA over the first 2")
  tr <- c(NA, NA, NA, y[4:8])
  expect_error(new_macrofilter(tr, y, base_meta(2L)), "1 non-finite values")
  expect_silent(new_macrofilter(c(NA, NA, y[3:8]), y, base_meta(2L)))
})

test_that("validator rejects non-finite trends", {
  y <- as.double(1:6)
  expect_error(new_macrofilter(c(1, 2, Inf, 4, 5, 6), y, base_meta()), "non-finite")
  expect_error(new_macrofilter(c(1, 2, NaN, 4, 5, 6), y, base_meta()), "non-finite")
})

test_that("validator checks the bands", {
  y <- as.double(1:6); tr <- y
  expect_error(new_macrofilter(tr, y, base_meta(),
                               bands = list(lower = tr - 1, upper = tr - 2)),
               "inverted")
  expect_error(new_macrofilter(tr, y, base_meta(),
                               bands = list(lower = tr[-1], upper = tr[-1])),
               "same length")
  obj <- new_macrofilter(tr, y, base_meta(),
                         bands = list(lower = tr - 1, upper = tr + 1))
  obj$trend_upper <- NULL
  expect_error(validate_macrofilter(obj), "both present or both absent")
  expect_error(new_macrofilter(tr, y, base_meta(),
                               bands = list(lower = c(NA, tr[-1]), upper = tr + 1)),
               "trend_lower")
})

# ── 3. Every filter declares its lead-in ─────────────────────────────────────

test_that("filters set meta$na_lead and keep it out of print()", {
  set.seed(1)
  y <- ts(cumsum(rnorm(60)) + 10, frequency = 4)
  expect_equal(hp_filter(y)$meta$na_lead, 0L)
  expect_equal(bhp_filter(y)$meta$na_lead, 0L)
  expect_equal(suppressMessages(mbh_filter(y, mstop = 20))$meta$na_lead, 0L)
  r <- hamilton_filter(y, h = 8, p = 4)
  expect_equal(r$meta$na_lead, 11L)
  expect_equal(sum(is.na(r$trend)), 11L)
  expect_false(any(grepl("na_lead", capture.output(print(r)))))
})
