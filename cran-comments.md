## Test environments

* local Ubuntu 24.04 (WSL2), R 4.5.1
* win-builder (devel and release) — results pending
* R-hub (ubuntu-latest, macos-latest) — results pending

## R CMD check results

0 errors | 0 warnings | 3 notes

### Note 1 — CRAN incoming feasibility (expected for new submission)

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

### Note 2 — Future file timestamps

```
unable to verify current time
```

This is a transient network issue on the local check machine (NTP server unreachable).
It does not reflect a problem with the package.

### Note 3 — HTML manual validation

```
Skipping checking HTML validation: no command 'tidy' found.
Skipping checking math rendering: package 'V8' unavailable.
```

This is a local environment limitation (`tidy` and `V8` not installed).
CRAN check servers have these tools available.

## Downstream dependencies

None — this is the first CRAN submission of this package.
