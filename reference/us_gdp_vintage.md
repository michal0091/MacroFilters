# US Real GDP — FRED Vintage

Quarterly US Real Gross Domestic Product from the Federal Reserve Bank
of St. Louis (FRED) public data API (series **GDPC1**), expressed in
billions of chained 2017 US dollars, seasonally adjusted annual rate.

## Usage

``` r
us_gdp_vintage
```

## Format

A `data.table` with one row per quarter and three columns:

- `date`:

  `Date`. Quarter start date (e.g. `1947-01-01` = 1947 Q1).

- `gdp_real`:

  `numeric`. Real GDP level, billions of chained 2017 USD.

- `gdp_log`:

  `numeric`. Natural logarithm of `gdp_real`, pre-computed for
  convenience.

## Source

Federal Reserve Bank of St. Louis — FRED Economic Data, series GDPC1.
Downloaded via the public CSV endpoint
`https://fred.stlouisfed.org/graph/fredgraph.csv?id=GDPC1`. See
`data-raw/us_gdp_vintage.R` for the reproducible download script.

## Details

The dataset covers 1947 Q1 through the latest vintage available at
download time (approximately 316 rows as of 2025). Rows with `NA` in
`gdp_real` are excluded.

The log-level column (`gdp_log`) is particularly useful for trend-cycle
decomposition because log-differences approximate quarter-on-quarter
percentage growth rates: \$\$\Delta \log(\text{GDP}\_t) \approx g_t\$\$

## Examples

``` r
data("us_gdp_vintage", package = "MacroFilters")
head(us_gdp_vintage)
#>         date gdp_real  gdp_log
#> 1 1947-01-01 2182.681 7.688309
#> 2 1947-04-01 2176.892 7.685653
#> 3 1947-07-01 2172.432 7.683603
#> 4 1947-10-01 2206.452 7.699141
#> 5 1948-01-01 2239.682 7.714089
#> 6 1948-04-01 2276.690 7.730478
plot(us_gdp_vintage$date, us_gdp_vintage$gdp_log,
     type = "l", xlab = "Date", ylab = "Log Real GDP",
     main = "US Real GDP (log level)")
```
