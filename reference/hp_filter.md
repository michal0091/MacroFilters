# Hodrick-Prescott Filter (Sparse Matrix Implementation)

Decomposes a time series into trend and cycle components by solving the
HP penalized least-squares problem using a sparse Cholesky
factorization. This avoids the dense O(n^3) inversion used by other
implementations and scales linearly in the number of observations.

## Usage

``` r
hp_filter(x, lambda = NULL, freq = NULL, boot_iter = 0, block_size = "auto")
```

## Arguments

- x:

  Numeric vector, `ts`, `xts`, or `zoo` object.

- lambda:

  Smoothing parameter. If `NULL` (default), it is auto-detected using
  the Ravn-Uhlig rule (`6.25 * freq^4`).

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
`data`, and `meta`. When `boot_iter > 0` it also carries `trend_lower`
and `trend_upper` (95% normal-approximation bootstrap band).

## Details

The HP filter minimises \$\$\sum (y_t - \tau_t)^2 + \lambda \sum
(\Delta^2 \tau_t)^2\$\$ which admits the closed-form solution \$\$(I +
\lambda D'D)\\\tau = y\$\$ where \\D\\ is the second-difference
operator.

The implementation builds \\D\\ as a banded sparse matrix
([`Matrix::bandSparse()`](https://rdrr.io/pkg/Matrix/man/bandSparse.html))
and solves the symmetric positive-definite system with a sparse Cholesky
decomposition
([`Matrix::solve()`](https://rdrr.io/pkg/Matrix/man/solve-methods.html)).

When `lambda` is not supplied the Ravn-Uhlig (2002) rule is applied:
`lambda = 6.25 * freq^4`, yielding 6.25 (annual), 1600 (quarterly), and
129 600 (monthly).

## References

Hodrick, R.J. and Prescott, E.C. (1997). Postwar U.S. Business Cycles:
An Empirical Investigation. *Journal of Money, Credit and Banking*,
29(1), 1–16.

Ravn, M.O. and Uhlig, H. (2002). On Adjusting the Hodrick-Prescott
Filter for the Frequency of Observations. *Review of Economics and
Statistics*, 84(2), 371–376.

## Examples

``` r
# Quarterly GDP-like series
y <- ts(cumsum(rnorm(200)), start = c(2000, 1), frequency = 4)
result <- hp_filter(y)
print(result)
#> -- MacroFilter [HP] --
#>    Observations : 200
#>    Parameters   : lambda = 1600
#>    Cycle range  : [-2.803, 4.151]  sd = 1.189
#>    Compute time : 0.002 s
```
