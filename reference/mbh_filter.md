# MacroBoost Hybrid (MBH) Filter

Decomposes a time series into trend and cycle using a robust boosting
algorithm. Unlike the HP filter, MBH uses the Huber loss function to
automatically downweight outliers (like the COVID-19 shock), preventing
them from distorting the trend.

## Usage

``` r
mbh_filter(
  x,
  d = "auto",
  boot_iter = 0,
  block_size = "auto",
  knots = NULL,
  mstop = 500L,
  nu = 0.1,
  df = 4L,
  select_mstop = FALSE,
  boundary.knots = NULL,
  hp_lambda = NULL
)
```

## Arguments

- x:

  Numeric vector, `ts`, `xts`, or `zoo` object.

- d:

  Numeric or `"auto"`. The delta parameter for Huber loss. If `"auto"`
  (default), it is calibrated as `stats::mad(.hp_fast(x))`, i.e. the MAD
  of the HP cyclical residual. This anchors the threshold to the
  output-gap scale rather than the growth-rate scale, avoiding the
  under-truncation failure mode of the legacy `mad(diff(y))` heuristic.
  A [`message()`](https://rdrr.io/r/base/message.html) is emitted
  reporting the exact value chosen. Supply an explicit positive numeric
  to override.

- boot_iter:

  Non-negative integer. Number of block-bootstrap iterations for
  uncertainty quantification (default `0`, bootstrap disabled). When
  `> 0`, the function adds `$trend_lower` and `$trend_upper`: a 95%
  normal-approximation band, `trend +/- 1.96 * sd(bootstrap trends)`,
  centred on the estimated trend. The bootstrap sd is used instead of
  empirical percentiles because it is smooth and stable at practical
  `boot_iter`. Each bootstrap refit uses the same `mstop` as the base
  fit, so larger `boot_iter` raises cost linearly. See also
  `block_size`.

- block_size:

  Positive integer or `"auto"`. Block length for the moving-block
  bootstrap (used only when `boot_iter > 0`). If `"auto"` (default), it
  is set to `2 * stats::frequency(x)` (two full cycles), bounded above
  by `floor(length(x) / 3)` to keep at least three blocks.

- knots:

  Integer. Number of interior knots for the P-Spline. If `NULL`
  (default), it is calculated as `max(20, floor(n / 2))`. High knot
  density is required for the trend to be flexible enough.

- mstop:

  Integer. Maximum number of boosting iterations (default 500). If
  `select_mstop = TRUE` this is the upper bound; the actual stopping
  point is chosen by AICc.

  **Under-smoothing warning (`mstop` vs `d`):** when `d` is small
  relative to the trend's range – the typical case for long log-level
  series, where the cycle (hence the auto-calibrated `d`) is tiny but
  the trend spans a large range – the Huber loss caps the gradient from
  the first iteration, so each boosting step advances the trend only
  slightly. Reducing `mstop` then leaves the trend unable to climb its
  full range: it collapses to a nearly flat curve while the cycle
  absorbs the long-run variation. Keep the default `mstop = 500` (or
  higher) for such series; lower it only for short or
  high-cycle-variance inputs.

- nu:

  Numeric. The learning rate (shrinkage) for boosting (default 0.1).

- df:

  Integer. Effective degrees of freedom per boosting step for the
  P-Spline base learner (default 4). This enforces the *weak-learner*
  constraint of Bühlmann & Hothorn (2007): each boosting step
  contributes only a small, smooth update so that the trend is built up
  gradually over many iterations rather than fitted in one pass.

  **End-point instability warning:** Higher `df` values cause the
  B-spline basis matrix to shift drastically when the sample size
  changes by even one observation (the "rubber-band effect"). The last
  few data points pull the estimated trend non-smoothly, producing
  unreliable end-of-sample estimates. Keep `df = 4` (the default) unless
  you have a specific reason to deviate.

- select_mstop:

  Logical. If `TRUE`, the optimal number of boosting iterations is
  selected automatically via AICc (corrected AIC), following Bühlmann &
  Hothorn (2007). The `mstop` argument acts as the search upper bound.
  Default `FALSE`.

  **AICc underfitting warning:** In the combination of Huber
  quasi-likelihood + P-splines, AICc penalises model complexity
  hyper-aggressively. In practice the algorithm stops at iteration ~5–15
  instead of the intended ~500. The resulting trend is nearly a straight
  line; all long-run variance is pushed into the cycle component,
  defeating the purpose of the filter. Treat `select_mstop = TRUE` as an
  experimental option and validate visually before relying on it.

- boundary.knots:

  A numeric vector of length 2 specifying the global domain for the
  B-spline basis (e.g., `c(1, T_max)`). If `NULL` (default), the range
  of `time_idx` is used. For real-time stability, fix this to the
  full-sample domain so the basis does not shift as the sample grows.

- hp_lambda:

  Numeric or `NULL`. Smoothing parameter for the internal HP filter used
  to auto-calibrate `d` (only relevant when `d = "auto"`). If `NULL`
  (default), it is derived from `stats::frequency(x)` via the Ravn-Uhlig
  rule. **Supply this when `x` is a plain numeric vector whose true
  frequency is not annual**, since
  [`frequency()`](https://rdrr.io/r/stats/time.html) returns `1` for
  unclassed vectors and would otherwise under-smooth the calibration
  cycle (e.g. monthly data: `hp_lambda = 129600`).

## Value

A list of class `c("macrofilter", "list")` with:

- `$trend`:

  Numeric trend vector.

- `$cycle`:

  Numeric cycle vector.

- `$data`:

  Original input as numeric.

- `$meta`:

  Named list: `method`, `knots`, `d`, `mstop`, `nu`, `df`,
  `select_mstop`, `compute_time`.

- `$trend_lower`, `$trend_upper`:

  95% normal-approximation bootstrap band (`trend +/- 1.96 * sd`).
  Present only when `boot_iter > 0`.

## Details

The model estimated is an additive model: \$\$y_t = \text{Linear}(t) +
\text{Smooth}(t) + \epsilon_t\$\$

It is fitted using
[`mboost::mboost()`](https://rdrr.io/pkg/mboost/man/gamboost.html) with:

- **Base Learners:** A linear time trend
  ([`mboost::bols()`](https://rdrr.io/pkg/mboost/man/baselearners.html))
  to capture the global path, plus a B-spline
  ([`mboost::bbs()`](https://rdrr.io/pkg/mboost/man/baselearners.html))
  to capture local curvature.

- **Loss Function:** Huber loss
  ([`mboost::Huber()`](https://rdrr.io/pkg/mboost/man/Family.html)) with
  parameter `d`. This is the key to robustness.

The default parameters (`knots = n/2`, `mstop = 500`) are calibrated to
mimic the flexibility of a standard HP filter while retaining the
robustness of the Huber loss.

## Calibration Guidance

Three failure modes were discovered through empirical stress-testing.
The defaults guard against all three:

- 1\. Huber delta scale mismatch (`d`):

  The automatic fallback `mad(diff(y))` operates on the scale of growth
  rates, not the output gap. For log-level input this sets `d` one to
  two orders of magnitude too small, causing ordinary business-cycle
  swings to be treated as outliers. If the estimated cycle looks
  implausibly large or the trend is nearly linear, override with
  `d = mad(hp_filter(x)$cycle)` as a starting point.

- 2\. AICc underfitting (`select_mstop`):

  AICc + Huber quasi-likelihood + P-splines stops boosting at iteration
  ~5–15. The trend degenerates to a near-straight line and the cycle
  absorbs all long-run variance. Leave `select_mstop = FALSE` (the
  default) and set `mstop` explicitly instead.

- 3\. End-point instability (`df`):

  Values above 4 shift the B-spline basis matrix non-smoothly as the
  sample grows, producing a "rubber-band" distortion in the final
  observations. Keep `df = 4` (the default) for real-time applications.

## Examples

``` r
# Fast example with reduced series and iterations
set.seed(42)
y <- ts(cumsum(rnorm(80)), start = c(2000, 1), frequency = 4)
result <- mbh_filter(y, mstop = 100L)
#> Info: Huber threshold automatically calibrated to d = 1.627402 via HP cyclical MAD.
print(result)
#> -- MacroFilter [MBH] --
#>    Observations : 80
#>    Parameters   : knots = 40, d = 1.627, mstop = 100, mstop_initial = 100, nu = 0.1, df = 4, select_mstop = FALSE
#>    Cycle range  : [-2.66, 4.533]  sd = 1.518
#>    Compute time : 0.026 s

# \donttest{
# Full example with default parameters
y2 <- ts(cumsum(rnorm(200)), start = c(2000, 1), frequency = 4)
result2 <- mbh_filter(y2)
#> Info: Huber threshold automatically calibrated to d = 1.161531 via HP cyclical MAD.
print(result2)
#> -- MacroFilter [MBH] --
#>    Observations : 200
#>    Parameters   : knots = 100, d = 1.162, mstop = 500, mstop_initial = 500, nu = 0.1, df = 4, select_mstop = FALSE
#>    Cycle range  : [-4.077, 3.726]  sd = 1.419
#>    Compute time : 0.143 s
# }
```
