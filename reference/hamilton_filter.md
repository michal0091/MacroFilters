# Hamilton Filter

Decomposes a time series into trend and cycle components using the
regression-based filter proposed by Hamilton (2018). The trend is the
fitted value from an OLS regression of \\y\_{t+h}\\ on \\(1, y_t,
y\_{t-1}, \ldots, y\_{t-p+1})\\, and the cycle is the residual.

## Usage

``` r
hamilton_filter(x, h = NULL, p = 4L, boot_iter = 0, block_size = "auto")
```

## Arguments

- x:

  Numeric vector, `ts`, `xts`, or `zoo` object.

- h:

  Integer horizon (number of periods ahead). If `NULL` (default),
  auto-detected from the series frequency using Hamilton's rule: annual
  = 2, quarterly = 8, monthly = 24.

- p:

  Integer number of lags in the regression (default 4).

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

## Value

A list of class `c("macrofilter", "list")` with `trend`, `cycle`,
`data`, and `meta` (`h`, `p`, `coefficients`, `compute_time`). When
`boot_iter > 0` it also carries `trend_lower` and `trend_upper` (95%
normal-approximation band). The bootstrap is a residual bootstrap
conditional on the initial `h + p - 1` observations (the regression
lead-in is held fixed); band entries for those lead-in positions are
`NA`.

## Details

Hamilton (2018) proposes replacing the HP filter with a simple
regression: \$\$y\_{t+h} = \beta_0 + \beta_1 y_t + \beta_2 y\_{t-1} +
\cdots + \beta_p y\_{t-p+1} + v\_{t+h}\$\$ The fitted values
\\\hat{y}\_{t+h}\\ define the trend and the residuals \\\hat{v}\_{t+h}\\
define the cycle.

The first \\h + p - 1\\ observations have no computable trend or cycle
and are filled with `NA`.

The lag matrix is constructed vectorized via
[`embed()`](https://rdrr.io/r/stats/embed.html) and the regression is
solved with [`stats::lm.fit()`](https://rdrr.io/r/stats/lmfit.html) for
speed.

When `boot_iter > 0`, the confidence band comes from a residual
bootstrap that holds the observed lead-in fixed (conditional on initial
values) and resamples residuals only over the valid window. A direct
consequence is that the band is narrow at the start of the valid window
– where the predictors are entirely the (frozen) lead-in, so only the
regression coefficients vary across replicates – and widens forward as
the predictors themselves become resampled quantities. This is the
correct behaviour of a conditional bootstrap, not an artefact.

## References

Hamilton, J.D. (2018). Why You Should Never Use the Hodrick-Prescott
Filter. *Review of Economics and Statistics*, 100(5), 831–843.

## Examples

``` r
# Quarterly GDP-like series
y <- ts(cumsum(rnorm(200)), start = c(2000, 1), frequency = 4)
result <- hamilton_filter(y)
print(result)
#> -- MacroFilter [Hamilton] --
#>    Observations : 200
#>    Parameters   : h = 8, p = 4
#>    Cycle range  : [-6.344, 4.538]  sd = 2.119
#>    Compute time : 0.001 s
```
