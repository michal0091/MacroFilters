# Spain Real GDP — FRED Vintage

Quarterly Spain Real Gross Domestic Product from the Federal Reserve
Bank of St. Louis (FRED) public data API (series **CLVMNACSCAB1GQES**),
expressed in millions of chained 2015 EUR, seasonally adjusted. Original
source: Eurostat/OECD National Accounts via FRED.

## Usage

``` r
es_gdp
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
CLVMNACSCAB1GQES. Downloaded via the public CSV endpoint
`https://fred.stlouisfed.org/graph/fredgraph.csv?id=CLVMNACSCAB1GQES`.
See `data-raw/intl_gdp.R` for the reproducible download script.

## Details

Spain experienced one of the sharpest COVID-19 contractions in the EU:
approximately \\-18\\\\ quarter-on-quarter in 2020 Q2, followed by a
strong V-shaped recovery in 2020 Q3. The series starts in 1995 Q1 under
ESA 2010 methodology.

## Examples

``` r
data("es_gdp", package = "MacroFilters")
head(es_gdp)
#>          date gdp_real  gdp_log
#>        <Date>    <num>    <num>
#> 1: 1995-01-01 178588.8 12.09284
#> 2: 1995-04-01 179766.0 12.09941
#> 3: 1995-07-01 180561.7 12.10383
#> 4: 1995-10-01 181928.6 12.11137
#> 5: 1996-01-01 182974.5 12.11710
#> 6: 1996-04-01 184239.5 12.12399
```
