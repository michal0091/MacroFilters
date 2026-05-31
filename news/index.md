# Changelog

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
