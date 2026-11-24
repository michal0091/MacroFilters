# Golden fixtures (roadmap task 4.2)
#
# Freezes the numerical output of MacroFilters 0.2.1 so that refactors in
# milestone 0A (0.2.2) can be checked to leave `ts` and `numeric` inputs
# untouched. Run from a checkout of tag v0.2.1 (R/ in release-0.2.2 was
# verified identical to v0.2.1 when these files were generated):
#
#   Rscript data-raw/golden_fixtures.R
#
# Output: tests/testthat/fixtures/{hp,hamilton,bhp,mbh}_{us,fr,es}.rds
#   Each file is a list with
#     $ts   trend for a ts(freq = 4) input, default arguments
#     $num  trend for a plain numeric input, default arguments (this is the
#           frequency-1 path of mbh_filter() that 0.2.b changes in 0.3.0)
#   plus tests/testthat/fixtures/boot_{hp,hamilton,bhp,mbh}_us.rds with the
#   bootstrap bands for a fixed seed.
# All values are rounded to 10 decimals.

devtools::load_all(quiet = TRUE)

out_dir <- file.path("tests", "testthat", "fixtures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

datasets <- list(us = us_gdp_vintage, fr = fr_gdp, es = es_gdp)
filters  <- list(hp = hp_filter, hamilton = hamilton_filter,
                 bhp = bhp_filter, mbh = mbh_filter)

as_quarterly_ts <- function(dt) {
  d0 <- as.POSIXlt(dt$date[1L])
  stats::ts(dt$gdp_log, start = c(d0$year + 1900, d0$mon %/% 3 + 1),
            frequency = 4)
}

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
r10   <- function(v) round(as.numeric(v), 10)

for (f in names(filters)) {
  for (k in names(datasets)) {
    y_ts  <- as_quarterly_ts(datasets[[k]])
    y_num <- as.numeric(datasets[[k]]$gdp_log)
    fx <- list(
      ts  = r10(quiet(filters[[f]](y_ts))$trend),
      num = r10(quiet(filters[[f]](y_num))$trend)
    )
    saveRDS(fx, file.path(out_dir, sprintf("%s_%s.rds", f, k)), version = 2)
  }
}

boot_iter <- c(hp = 50L, hamilton = 50L, bhp = 50L, mbh = 20L)
y_us <- as_quarterly_ts(us_gdp_vintage)
for (f in names(filters)) {
  set.seed(20260521)
  fit <- quiet(filters[[f]](y_us, boot_iter = boot_iter[[f]]))
  fx <- list(boot_iter = boot_iter[[f]], seed = 20260521L,
             trend = r10(fit$trend),
             lower = r10(fit$trend_lower), upper = r10(fit$trend_upper))
  saveRDS(fx, file.path(out_dir, sprintf("boot_%s_us.rds", f)), version = 2)
}
