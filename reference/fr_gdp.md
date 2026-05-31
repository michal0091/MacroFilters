# France Real GDP — FRED Vintage

Quarterly France Real Gross Domestic Product from the Federal Reserve
Bank of St. Louis (FRED) public data API (series **CLVMNACSCAB1GQFR**),
expressed in millions of chained 2015 EUR, seasonally adjusted. Original
source: Eurostat/OECD National Accounts via FRED.

## Usage

``` r
fr_gdp
```

## Format

A `data.table` with one row per quarter and three columns:

- `date`:

  `Date`. Quarter start date (e.g. `2000-01-01` = 2000 Q1).

- `gdp_real`:

  `numeric`. Real GDP level, millions of chained 2015 EUR.

- `gdp_log`:

  `numeric`. Natural logarithm of `gdp_real`, pre-computed for
  convenience.

## Source

Federal Reserve Bank of St. Louis — FRED Economic Data, series
CLVMNACSCAB1GQFR. Downloaded via the public CSV endpoint
`https://fred.stlouisfed.org/graph/fredgraph.csv?id=CLVMNACSCAB1GQFR`.
See `data-raw/intl_gdp.R` for the reproducible download script.

## Details

France experienced a sharp COVID-19 contraction of approximately
\\-14\\\\ quarter-on-quarter in 2020 Q2, followed by a rapid V-shaped
recovery. Together with Spain (`es_gdp`), the two series serve as a
demanding stress test for trend filters in the introduction vignette.

## Examples

``` r
data("fr_gdp", package = "MacroFilters")
head(fr_gdp)
#>          date gdp_real  gdp_log
#>        <Date>    <num>    <num>
#> 1: 1980-01-01 276309.0 12.52928
#> 2: 1980-04-01 274083.4 12.52119
#> 3: 1980-07-01 274497.8 12.52270
#> 4: 1980-10-01 274069.2 12.52114
#> 5: 1981-01-01 275088.0 12.52485
#> 6: 1981-04-01 277168.0 12.53238
```
