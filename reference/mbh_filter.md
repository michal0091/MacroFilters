# MacroBoost Hybrid (MBH) Filter

Decomposes a time series into trend and cycle using a robust boosting
algorithm. Unlike the HP filter, MBH uses the Huber loss function to
automatically downweight outliers (like the COVID-19 shock), preventing
them from distorting the trend.

## Usage

``` r
mbh_filter(x, knots = NULL, mstop = 500L, d = NULL, nu = 0.2)
```

## Arguments

- x:

  Numeric vector, `ts`, `xts`, or `zoo` object.

- knots:

  Integer. Number of interior knots for the P-Spline. If `NULL`
  (default), it is calculated as `max(20, floor(n / 2))`. High knot
  density is required for the trend to be flexible enough.

- mstop:

  Integer. Number of boosting iterations (default 500). Robust loss
  functions require more iterations to converge than L2 loss.

- d:

  Numeric or `NULL`. The delta parameter for Huber loss. If `NULL`
  (default), it is auto-calibrated as the Median Absolute Deviation
  (MAD) of the first differences of the series. This makes the threshold
  robust and scale-invariant, working correctly for both log-differenced
  and level series. Supply an explicit numeric value to override this
  behaviour.

- nu:

  Numeric. The learning rate (shrinkage) for boosting (default 0.2).

## Value

A `macrofilter` object with `trend`, `cycle`, `data`, and `meta`
components. The `meta` list contains `method = "MBH"`, `knots`, `d`,
`mstop`, `nu`, and `compute_time`.

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
#>    Compute time : 0.202 s
```
