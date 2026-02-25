# MacroBoost Hybrid (MBH) Filter

Decomposes a time series into trend and cycle using a robust boosting
algorithm. Unlike the HP filter, MBH uses the Huber loss function to
automatically downweight outliers (like the COVID-19 shock), preventing
them from distorting the trend.

## Usage

``` r
mbh_filter(
  x,
  knots = NULL,
  mstop = 500L,
  d = NULL,
  nu = 0.1,
  df = 4L,
  select_mstop = FALSE
)
```

## Arguments

- x:

  Numeric vector, `ts`, `xts`, or `zoo` object.

- knots:

  Integer. Number of interior knots for the P-Spline. If `NULL`
  (default), it is calculated as `max(20, floor(n / 2))`. High knot
  density is required for the trend to be flexible enough.

- mstop:

  Integer. Maximum number of boosting iterations (default 500). If
  `select_mstop = TRUE` this is the upper bound; the actual stopping
  point is chosen by AICc.

- d:

  Numeric or `NULL`. The delta parameter for Huber loss. If `NULL`
  (default), it is auto-calibrated as the Median Absolute Deviation
  (MAD) of the first differences of the series (`mad(diff(y))`).

  **Scale-mismatch warning:** When `x` is a log-level series, `diff(y)`
  returns inter-period growth rates (typical scale 0.001–0.02), whereas
  the output gap (the residual the filter must explain) has a much
  larger scale (typical scale 0.01–0.05). Using `mad(diff(y))` as `d`
  therefore sets the Huber threshold far too low: the filter treats
  normal business-cycle oscillations as outliers, truncates their
  gradients, and blocks learning. The trend becomes over-smooth and the
  cycle absorbs too much long-run variance.

  **Recommendation:** For a quick preliminary calibration use
  `d = mad(hp_filter(x)$cycle)`, which sets the threshold on the
  residual scale. Supply an explicit numeric value to override the
  automatic fallback entirely.

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

## Value

A `macrofilter` object with `trend`, `cycle`, `data`, and `meta`
components. The `meta` list contains `method = "MBH"`, `knots`, `d`,
`mstop`, `nu`, `df`, `select_mstop`, and `compute_time`.

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
# Quarterly GDP-like series
set.seed(42)
y <- ts(cumsum(rnorm(200)), start = c(2000, 1), frequency = 4)
result <- mbh_filter(y)
print(result)
#> -- MacroFilter [MBH] --
#>    Observations : 200
#>    Parameters   : knots = 100, d = 0.9445, mstop = 500, mstop_initial = 500, nu = 0.1, df = 4, select_mstop = FALSE
#>    Cycle range  : [-4.17, 5.05]  sd = 1.721
#>    Compute time : 0.204 s
```
