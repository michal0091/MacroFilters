# Boosted HP Filter

Iteratively applies the Hodrick-Prescott filter on residuals to better
capture stochastic trends. At each iteration the HP smoother is applied
to the current residual and the resulting trend increment is added to
the cumulative trend estimate. Iteration stops according to one of three
rules: BIC minimisation (default), ADF stationarity test on residuals,
or a fixed number of iterations.

## Usage

``` r
bhp_filter(
  x,
  lambda = NULL,
  iter_max = 100L,
  stopping = c("bic", "adf", "fixed"),
  sig_level = 0.05,
  freq = NULL,
  boot_iter = 0,
  block_size = "auto"
)
```

## Arguments

- x:

  Numeric vector, `ts`, `xts`, or `zoo` object.

- lambda:

  Smoothing parameter. If `NULL` (default), it is auto-detected using
  the Ravn-Uhlig rule (`6.25 * freq^4`).

- iter_max:

  Integer. Maximum number of boosting iterations (default 100).

- stopping:

  Character. Stopping rule: `"bic"` (default), `"adf"`, or `"fixed"`.

- sig_level:

  Numeric. Significance level for the ADF test when `stopping = "adf"`
  (default 0.05).

- freq:

  Numeric frequency override (1 = annual, 4 = quarterly, 12 = monthly).
  Used only when `lambda` is `NULL` and the frequency cannot be inferred
  from `x`.

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
`data`, and `meta` (`method = "bHP"`, `lambda`, `iterations`,
`stopping_rule`, `compute_time`). When `boot_iter > 0` it also carries
`trend_lower` and `trend_upper` (95% normal-approximation bootstrap
band); each bootstrap refit runs a *fixed* `iterations` passes,
conditioning on the complexity selected by the base fit.

## Details

The boosted HP filter starts from the standard HP solution and then
re-applies the same HP smoother to the residual (cycle) component. The
trend increment from each pass is accumulated, and the procedure stops
when one of the following criteria is met:

- `"bic"`:

  Schwarz information criterion computed as \\n \log(\hat\sigma^2) +
  \log(n)\\\mathrm{tr}(S^m)\\, where \\S^m\\ is the iterated smoother.
  Iteration stops when the BIC increases relative to the previous best.

- `"adf"`:

  Augmented Dickey-Fuller test on the residual. Iteration stops when the
  residual is stationary at level `sig_level`.

- `"fixed"`:

  Runs exactly `iter_max` iterations.

## References

Phillips, P.C.B. and Shi, Z. (2021). Boosting: Why You Can Use the HP
Filter. *International Economic Review*, 62(2), 521–570.

## Examples

``` r
# Quarterly GDP-like series
y <- ts(cumsum(rnorm(200)), start = c(2000, 1), frequency = 4)
result <- bhp_filter(y)
print(result)
#> -- MacroFilter [bHP] --
#>    Observations : 200
#>    Parameters   : lambda = 1600, iterations = 100, stopping_rule = bic
#>    Cycle range  : [-1.657, 1.673]  sd = 0.6962
#>    Compute time : 0.006 s
```
