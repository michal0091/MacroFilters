## Submission

This is a maintenance update of MacroFilters, from 0.1.0 (the current CRAN
version) to 0.2.1. The main user-visible change bounds the automatic knot
count of `mbh_filter()` to keep the B-spline basis from growing without limit
on long / high-frequency series; quarterly macro series are unaffected. See
NEWS.md for the full list. (Version 0.2.0 was a GitHub-only development release
and was not submitted to CRAN.)

## Test environments

* local: Ubuntu 24.04 (WSL2), R 4.6.0
* win-builder: R-devel (2026-06-09 r90126 ucrt) and R-release (4.6.0) — Status: OK
* R-hub v2 (GitHub Actions):
  - linux, R-devel (2026-06-09 r90126), Ubuntu 24.04 — Status: OK
  - windows, R-devel (2026-06-09 r90126 ucrt), Windows Server 2022 — Status: OK
  - macos-arm64, R-devel (2026-05-03 r89994), macOS Sequoia 15.7 (aarch64) — Status: OK

## R CMD check results

0 errors | 0 warnings | 0 notes

on every environment listed above.

## Notes

CRAN's incoming check may flag the following domain-specific terms in
DESCRIPTION as possibly misspelled — all are correct and are listed in
`inst/WORDLIST`: AICc, Hodrick, Kinel, MBH, MacroBoost.

## Reverse dependencies

There are no reverse dependencies for this package on CRAN.
