## ── fr_gdp / es_gdp datasets ─────────────────────────────────────────────────
## Reproducible download script.  Run once to (re)generate the bundled .rda files.
##
## France: FRED series CLVMNACSCAB1GQFR
##         Real GDP, millions of chained 2015 EUR, SA quarterly
##         Source: Eurostat/OECD via FRED — https://fred.stlouisfed.org/series/CLVMNACSCAB1GQFR
##
## Spain:  FRED series CLVMNACSCAB1GQES
##         Real GDP, millions of chained 2015 EUR, SA quarterly
##         Source: Eurostat/OECD via FRED — https://fred.stlouisfed.org/series/CLVMNACSCAB1GQES
##
## No API key required; uses the FRED public CSV endpoint.

library(data.table)

fetch_fred <- function(id) {
  url <- sprintf("https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s", id)
  dt  <- data.table::fread(url, col.names = c("date", "gdp_real"), na.strings = ".")
  dt[, date     := as.Date(date)]
  dt[, gdp_real := as.numeric(gdp_real)]
  dt[, gdp_log  := log(gdp_real)]
  dt[!is.na(gdp_real)]
}

fr_gdp <- fetch_fred("CLVMNACSCAB1GQFR")
es_gdp <- fetch_fred("CLVMNACSCAB1GQES")

usethis::use_data(fr_gdp, overwrite = TRUE)
usethis::use_data(es_gdp, overwrite = TRUE)

message(sprintf("fr_gdp: %d rows, %s to %s", nrow(fr_gdp),
                format(min(fr_gdp$date)), format(max(fr_gdp$date))))
message(sprintf("es_gdp: %d rows, %s to %s", nrow(es_gdp),
                format(min(es_gdp$date)), format(max(es_gdp$date))))
