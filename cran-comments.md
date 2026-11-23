## Submission

This is a maintenance update of MacroFilters, from 0.2.1 (the current CRAN
version) to 0.2.2. It fixes how `xts` and `zoo` input is handled, so that it
gives the same results as the equivalent `ts`, rejects multi-column input,
and corrects the package description, which claimed a diagnostic stage the
package does not have. Results for `ts` and numeric input are unchanged.
See NEWS.md for the full list.

## Test environments

* local: Ubuntu 24.04 (WSL2), R 4.6.2 — Status: TODO
* win-builder: R-devel and R-release — Status: TODO
* R-hub v2 (GitHub Actions): linux, windows, macos-arm64 (R-devel) — Status: TODO

## R CMD check results

TODO

## Notes

CRAN's incoming check may flag the following domain-specific terms in
DESCRIPTION as possibly misspelled — all are correct and are listed in
`inst/WORDLIST`: AICc, Hodrick, Kinel, MBH, MacroBoost.

## Reverse dependencies

There are no reverse dependencies for this package on CRAN.
