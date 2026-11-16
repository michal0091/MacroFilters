# MacroFilters (development version)

## Bug fixes

* `xts` and `zoo` inputs now get the same defaults as the equivalent `ts`.
  The frequency is read from a regular `Date`, `POSIXct`, `yearqtr` or
  `yearmon` index. Before, a `Date`-indexed series had no frequency:
  `hp_filter()`, `bhp_filter()` and `hamilton_filter()` warned and assumed
  quarterly, while `mbh_filter()` (for `d = "auto"`) and the `"auto"`
  bootstrap block size treated it as annual.
* The `"auto"` bootstrap block size now follows the `freq` argument when `x`
  is a plain numeric vector. Before, it was always 2 for such vectors. A
  vector without `freq` still gets a block of 2.
* `zoo` objects of class `zooreg`, which is what `zoo::as.zoo()` returns for
  a `ts`, are now accepted. Before, they were rejected as an unsupported
  class.
* Multi-column input (`mts`, a multi-column `xts` or `zoo`, or a matrix) is
  now an error. Before, a multi-column `xts` or `zoo` was silently flattened
  into one long series. A single-column matrix is accepted as a vector.
* `autoplot()` now draws `xts` and `zoo` input against its own time index:
  a date axis for a `Date` or `POSIXct` index, decimal time for `yearqtr` or
  `yearmon`. Before, it plotted the positions `1..n`.

## New features

* `hamilton_filter()` and `mbh_filter()` gain a `freq` argument, as
  `hp_filter()` and `bhp_filter()` already had. It is used only when `x`
  carries no frequency of its own. A `freq` that contradicts the frequency of
  `x` is now ignored with a warning instead of silently.
* All filters store in `meta$freq` the frequency read from `x` or `freq`,
  or `NA` when it is unknown.

# MacroFilters 0.2.1

## Other changes

* `mbh_filter()`'s automatic knot count is now capped at 250
  (`min(max(20, floor(n / 2)), 250)`). Series of 500 observations or fewer are
  unaffected; the cap only bounds the B-spline basis for long or high-frequency
  inputs, where extra knots inflate memory and runtime without adding
  flexibility (in a P-spline the difference penalty, not the knot count,
  controls smoothness).

## Documentation

* Corrected the MBH parameter tables in the *Introduction* vignette: the
  `d = "auto"` default is calibrated from the MAD of the HP cyclical residual
  (not first differences), and the default learning rate is `nu = 0.1`.
* Fixed the COVID-19 highlight in the *Introduction* cycle plot, which was
  anchored to stale fixed indices instead of the 2020 date window.

# MacroFilters 0.2.0

## New features

* Confidence bands via block bootstrap for **all four filters**
  (`hp_filter()`, `hamilton_filter()`, `bhp_filter()`, `mbh_filter()`). The new
  `boot_iter` and `block_size` arguments add `$trend_lower` / `$trend_upper`
  to the result: a 95% normal-approximation band (`trend ± 1.96 * sd`) built
  from a Circular Block Bootstrap of the cycle, with each replicate refit by the
  same estimator as the base fit.
* New `autoplot()` method for `macrofilter` objects (ggplot2): draws the
  observed series, the estimated trend, and the confidence ribbon when present,
  with the time axis reconstructed from the stored temporal identity.
* `mbh_filter()` gains `hp_lambda` to control the HP-based auto-calibration of
  the Huber threshold `d` when the input is a plain numeric vector whose true
  frequency is not annual.

## Performance

* The HP system matrix is now Cholesky-factorized **once** and reused across
  every bootstrap replicate (and every bHP inner iteration), instead of being
  re-factorized on each solve. This markedly speeds up `hp_filter()` and
  `bhp_filter()` with `boot_iter > 0` (and the base bHP fit), with bit-identical
  results.

## Other changes

* The `d = "auto"` calibration in `mbh_filter()` now uses the MAD of the HP
  cyclical residual (output-gap scale) instead of `mad(diff(y))`, and reports
  the chosen value via a `message()`.
* Filters now return a list of class `c("macrofilter", "list")` and store the
  temporal identity (`meta$ts_class`, `meta$tsp`, `meta$idx`) so trend, cycle
  and bands can all be mapped back to dates for plotting.

## Documentation

* New vignette *Uncertainty Bands via Block Bootstrap* covering `boot_iter`,
  `block_size`, the end-point fan and the Hamilton conditional band.
* `mbh_filter()` documents the `mstop`–`d` interaction (reducing `mstop` on long
  log-level series under-smooths the trend); `hamilton_filter()` documents the
  conditional bootstrap band behaviour.
