# Changelog

## MacroFilters 0.2.2

### Bug fixes

- `xts` and `zoo` inputs now get the same defaults as the equivalent
  `ts`. The frequency is read from a regular `Date`, `POSIXct`,
  `yearqtr` or `yearmon` index. Before, a `Date`-indexed series had no
  frequency:
  [`hp_filter()`](https://michal0091.github.io/MacroFilters/reference/hp_filter.md),
  [`bhp_filter()`](https://michal0091.github.io/MacroFilters/reference/bhp_filter.md)
  and
  [`hamilton_filter()`](https://michal0091.github.io/MacroFilters/reference/hamilton_filter.md)
  warned and assumed quarterly, while
  [`mbh_filter()`](https://michal0091.github.io/MacroFilters/reference/mbh_filter.md)
  (for `d = "auto"`) and the `"auto"` bootstrap block size treated it as
  annual.
- The `"auto"` bootstrap block size now follows the `freq` argument when
  `x` is a plain numeric vector. Before, it was always 2 for such
  vectors. A vector without `freq` still gets a block of 2.
- `zoo` objects of class `zooreg`, which is what
  [`zoo::as.zoo()`](https://rdrr.io/pkg/zoo/man/as.zoo.html) returns for
  a `ts`, are now accepted. Before, they were rejected as an unsupported
  class.
- Multi-column input (`mts`, a multi-column `xts` or `zoo`, or a matrix)
  is now an error. Before, a multi-column `xts` or `zoo` was silently
  flattened into one long series. A single-column matrix is accepted as
  a vector.
- [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  now draws `xts` and `zoo` input against its own time index: a date
  axis for a `Date` or `POSIXct` index, decimal time for `yearqtr` or
  `yearmon`. Before, it plotted the positions `1..n`.

### New features

- [`hamilton_filter()`](https://michal0091.github.io/MacroFilters/reference/hamilton_filter.md)
  and
  [`mbh_filter()`](https://michal0091.github.io/MacroFilters/reference/mbh_filter.md)
  gain a `freq` argument, as
  [`hp_filter()`](https://michal0091.github.io/MacroFilters/reference/hp_filter.md)
  and
  [`bhp_filter()`](https://michal0091.github.io/MacroFilters/reference/bhp_filter.md)
  already had. It is used only when `x` carries no frequency of its own.
  A `freq` that contradicts the frequency of `x` is now ignored with a
  warning instead of silently.
- All filters store in `meta$freq` the frequency read from `x` or
  `freq`, or `NA` when it is unknown.
- All filters store in `meta$na_lead` the number of leading positions
  with no estimate: `h + p - 1` for
  [`hamilton_filter()`](https://michal0091.github.io/MacroFilters/reference/hamilton_filter.md),
  `0` for the others.

### Other changes

- [`mbh_filter()`](https://michal0091.github.io/MacroFilters/reference/mbh_filter.md)
  no longer prints the calibrated Huber threshold when `d = "auto"`. The
  value is in `meta$d`; set `options(MacroFilters.verbose = TRUE)` to
  get the message back.

### Documentation

- The package description no longer claims an automated diagnostic stage
  for unit roots and structural breaks, which the package does not have.
- The README example used a `us_gdp_xts` object that the package does
  not ship. It now uses `us_gdp_vintage` from 2000 onwards.
- The BIC formula in the
  [`bhp_filter()`](https://michal0091.github.io/MacroFilters/reference/bhp_filter.md)
  help now uses the trace of the boosted smoother, `tr(B_m)`, which is
  what the code computes.
- The
  [`mbh_filter()`](https://michal0091.github.io/MacroFilters/reference/mbh_filter.md)
  help announces that `select_mstop` will be deprecated in 0.3.0.

### Internal changes

- All four filters now build their result through one constructor, which
  derives the cycle as `data - trend` and validates the object. The
  checks are now ones that can fail: equal lengths, a finite trend
  outside the `NA` lead-in, and well-formed bootstrap bands. The old
  check of `trend + cycle == data` held by construction. The cycle of
  [`bhp_filter()`](https://michal0091.github.io/MacroFilters/reference/bhp_filter.md)
  and
  [`hamilton_filter()`](https://michal0091.github.io/MacroFilters/reference/hamilton_filter.md)
  can change in the last floating-point digit; the trend and the bands
  do not change.
- Removed the unused `restore_class()` helpers.
- [`hp_filter()`](https://michal0091.github.io/MacroFilters/reference/hp_filter.md)
  and
  [`bhp_filter()`](https://michal0091.github.io/MacroFilters/reference/bhp_filter.md)
  factor the HP system once and reuse it for the fit and the bootstrap.
- Dropped the unused `strucchange` and `usethis` from Suggests, added
  `mFilter` to Suggests for a reference test, and set minimum versions
  for `data.table`, `ggplot2` and `Matrix`.

## MacroFilters 0.2.1

CRAN release: 2026-06-11

### Other changes

- [`mbh_filter()`](https://michal0091.github.io/MacroFilters/reference/mbh_filter.md)’s
  automatic knot count is now capped at 250
  (`min(max(20, floor(n / 2)), 250)`). Series of 500 observations or
  fewer are unaffected; the cap only bounds the B-spline basis for long
  or high-frequency inputs, where extra knots inflate memory and runtime
  without adding flexibility (in a P-spline the difference penalty, not
  the knot count, controls smoothness).

### Documentation

- Corrected the MBH parameter tables in the *Introduction* vignette: the
  `d = "auto"` default is calibrated from the MAD of the HP cyclical
  residual (not first differences), and the default learning rate is
  `nu = 0.1`.
- Fixed the COVID-19 highlight in the *Introduction* cycle plot, which
  was anchored to stale fixed indices instead of the 2020 date window.

## MacroFilters 0.2.0

### New features

- Confidence bands via block bootstrap for **all four filters**
  ([`hp_filter()`](https://michal0091.github.io/MacroFilters/reference/hp_filter.md),
  [`hamilton_filter()`](https://michal0091.github.io/MacroFilters/reference/hamilton_filter.md),
  [`bhp_filter()`](https://michal0091.github.io/MacroFilters/reference/bhp_filter.md),
  [`mbh_filter()`](https://michal0091.github.io/MacroFilters/reference/mbh_filter.md)).
  The new `boot_iter` and `block_size` arguments add `$trend_lower` /
  `$trend_upper` to the result: a 95% normal-approximation band
  (`trend ± 1.96 * sd`) built from a Circular Block Bootstrap of the
  cycle, with each replicate refit by the same estimator as the base
  fit.
- New
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  method for `macrofilter` objects (ggplot2): draws the observed series,
  the estimated trend, and the confidence ribbon when present, with the
  time axis reconstructed from the stored temporal identity.
- [`mbh_filter()`](https://michal0091.github.io/MacroFilters/reference/mbh_filter.md)
  gains `hp_lambda` to control the HP-based auto-calibration of the
  Huber threshold `d` when the input is a plain numeric vector whose
  true frequency is not annual.

### Performance

- The HP system matrix is now Cholesky-factorized **once** and reused
  across every bootstrap replicate (and every bHP inner iteration),
  instead of being re-factorized on each solve. This markedly speeds up
  [`hp_filter()`](https://michal0091.github.io/MacroFilters/reference/hp_filter.md)
  and
  [`bhp_filter()`](https://michal0091.github.io/MacroFilters/reference/bhp_filter.md)
  with `boot_iter > 0` (and the base bHP fit), with bit-identical
  results.

### Other changes

- The `d = "auto"` calibration in
  [`mbh_filter()`](https://michal0091.github.io/MacroFilters/reference/mbh_filter.md)
  now uses the MAD of the HP cyclical residual (output-gap scale)
  instead of `mad(diff(y))`, and reports the chosen value via a
  [`message()`](https://rdrr.io/r/base/message.html).
- Filters now return a list of class `c("macrofilter", "list")` and
  store the temporal identity (`meta$ts_class`, `meta$tsp`, `meta$idx`)
  so trend, cycle and bands can all be mapped back to dates for
  plotting.

### Documentation

- New vignette *Uncertainty Bands via Block Bootstrap* covering
  `boot_iter`, `block_size`, the end-point fan and the Hamilton
  conditional band.
- [`mbh_filter()`](https://michal0091.github.io/MacroFilters/reference/mbh_filter.md)
  documents the `mstop`–`d` interaction (reducing `mstop` on long
  log-level series under-smooths the trend);
  [`hamilton_filter()`](https://michal0091.github.io/MacroFilters/reference/hamilton_filter.md)
  documents the conditional bootstrap band behaviour.
