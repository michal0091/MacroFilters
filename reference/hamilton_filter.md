# Hamilton Filter

Decomposes a time series into trend and cycle components using the
regression-based filter proposed by Hamilton (2018). The trend is the
fitted value from an OLS regression of \\y\_{t+h}\\ on \\(1, y_t,
y\_{t-1}, \ldots, y\_{t-p+1})\\, and the cycle is the residual.

## Usage

``` r
hamilton_filter(x, h = NULL, p = 4L)
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

## Value

A `macrofilter` object with `trend`, `cycle`, `data`, and `meta`
components. `meta` includes `h`, `p`, `coefficients`, and
`compute_time`.

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
#>    Cycle range  : [-5.501, 6.003]  sd = 2.437
#>    Compute time : 0.000 s
```
