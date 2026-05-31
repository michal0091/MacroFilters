# ── Package-level declarations ────────────────────────────────────────────────

#' @importFrom data.table :=
NULL

# Silence R CMD check NOTEs for data.table symbols used via non-standard
# evaluation inside .boot_engine() (column names referenced bare in `[.data.table`).
utils::globalVariables(c(
  # .boot_engine() data.table symbols
  "iter", "time_idx", "y_boot", "boot_trend",
  # autoplot.macrofilter() melt + aes symbols
  "time", "value", "series", "lower", "upper"
))
