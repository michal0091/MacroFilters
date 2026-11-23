## Submission

This is a maintenance update of MacroFilters, from 0.2.1 (the current CRAN
version) to 0.2.2. It fixes how `xts` and `zoo` input is handled, so that it
gives the same results as the equivalent `ts`, rejects multi-column input,
and corrects the package description, which claimed a diagnostic stage the
package does not have. Results for `ts` and numeric input are unchanged.
See NEWS.md for the full list.

## Test environments

* local: Ubuntu 24.04 (WSL2), R 4.6.2
* win-builder:
  - R-devel (2026-11-17 r90926 ucrt) — Status: TODO
  - R-release (4.6.2) — Status: TODO
* R-hub v2 (GitHub Actions):
  - linux, R-devel (2026-11-21 r90959), Ubuntu 26.04.1 LTS — Status: OK
  - windows, R-devel (2026-11-17 r90926 ucrt), Windows Server 2022 — Status: OK
  - macos, R-devel (2026-11-23 r90962), macOS Sequoia 15.7.9 (x86_64) — Status: OK

## R CMD check results

0 errors | 0 warnings | 0 notes

on every environment listed above.

## Notes

CRAN's incoming check may flag the following domain-specific terms in
DESCRIPTION as possibly misspelled — all are correct: AICc, Hodrick, Kinel,
MBH, MacroBoost.

## Reverse dependencies

There are no reverse dependencies for this package on CRAN.
