## Test environments

* local Ubuntu 24.04 (WSL2), R 4.5.1
* win-builder R-devel (2026-05-20): 0 errors | 0 warnings | 1 note

## R CMD check results

0 errors | 0 warnings | 1 note

### Note — CRAN incoming feasibility (expected for new submission)

```
New submission

Possibly misspelled words in DESCRIPTION:
  AICc (15:67), Hodrick (11:5), Kinel (16:55), MBH (9:61), MacroBoost (9:42)
```

This is a first CRAN submission.

The flagged words are all legitimate domain-specific terms:
- **AICc**: corrected Akaike Information Criterion (standard statistical abbreviation)
- **Hodrick**: proper name, co-inventor of the Hodrick-Prescott filter
- **Kinel**: author's surname in a self-citation
- **MBH**: acronym for MacroBoost Hybrid (the package's flagship filter)
- **MacroBoost**: the package's gradient-boosting methodology

These terms are included in `inst/WORDLIST`.

## Downstream dependencies

None — this is the first CRAN submission of this package.
