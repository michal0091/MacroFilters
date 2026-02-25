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
  freq = NULL
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

## Value

A `macrofilter` object with `trend`, `cycle`, `data`, and `meta`
components. The `meta` list contains `method = "bHP"`, `lambda`,
`iterations`, `stopping_rule`, and `compute_time`.

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
#>    Parameters   : lambda = 1600, iterations = 93, stopping_rule = bic
#>    Cycle range  : [-1.728, 1.688]  sd = 0.6928
#>    Compute time : 0.018 s
```
