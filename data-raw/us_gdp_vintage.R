## ── us_gdp_vintage dataset ────────────────────────────────────────────────────
## Reproducible download script.  Run once to (re)generate data/us_gdp_vintage.rda
##
## Source : FRED series GDPC1 — Real Gross Domestic Product
##          Billions of Chained 2017 Dollars, Seasonally Adjusted Annual Rate
## URL    : https://fred.stlouisfed.org/series/GDPC1
## Period : 1947 Q1 — latest available (~308+ rows)
##
## No API key required; uses the FRED public CSV endpoint.

library(data.table)

url <- "https://fred.stlouisfed.org/graph/fredgraph.csv?id=GDPC1"
dt  <- data.table::fread(url, col.names = c("date", "gdp_real"))
dt[, date     := as.Date(date)]
dt[, gdp_real := as.numeric(gdp_real)]
dt[, gdp_log  := log(gdp_real)]
us_gdp_vintage <- dt[!is.na(gdp_real)]

usethis::use_data(us_gdp_vintage, overwrite = TRUE)

message(sprintf(
  "us_gdp_vintage: %d rows, %s to %s",
  nrow(us_gdp_vintage),
  format(min(us_gdp_vintage$date)),
  format(max(us_gdp_vintage$date))
))
