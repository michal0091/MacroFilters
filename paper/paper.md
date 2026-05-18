---
title: 'MacroFilters: High-Performance Robust Trend Extraction in R'
tags:
  - R
  - time series
  - macroeconomics
  - econometrics
  - gradient boosting
authors:
  - name: Michal Kinel
    orcid: 0009-0007-3295-7199
    affiliation: 1
affiliations:
 - name: Independent Researcher
   index: 1
date: 18 May 2026
bibliography: paper.bib
---

# Summary

Extracting the underlying structural trend from macroeconomic time series is a foundational task in econometrics. However, real-time estimation is chronically compromised by the "end-point problem": traditional $L_2$-loss filters mechanically deform the long-run trend to absorb recent, transitory outliers (e.g., the COVID-19 shock). `MacroFilters` is an R package designed to solve this instability by introducing the MacroBoost Hybrid (MBH) filter. It integrates penalized P-splines within a gradient boosting architecture, optimized via an adaptive Huber quasi-likelihood. The result is a high-performance, robust filtering tool that preserves trend integrity during black swan events without sacrificing baseline accuracy.

# Statement of need

In the R ecosystem, practitioners typically rely on legacy packages (e.g., `mFilter`) to apply standard tools like the Hodrick-Prescott (HP) filter [@hodrick1997postwar]. However, these methods are notoriously unreliable at the end of the sample [@hamilton2018never], leading to severe backward revisions as new data arrives. This forces analysts and policymakers to continuously rewrite recent economic history.

`MacroFilters` addresses this gap by providing a modern, machine-learning-driven alternative explicitly designed for real-time robustness. The core algorithms are optimized using `data.table` [@dowle2021datatable] to handle large-scale datasets and high-frequency simulations with minimal memory overhead.

The package provides a streamlined API built around the `mbh_filter()` function, seamlessly accepting standard R time-series objects (`ts`, `zoo`, `xts`). To ensure immediate usability, `MacroFilters` includes an extensive documentation website built with `pkgdown`, featuring comprehensive vignettes on basic usage, hyperparameter tuning, and visualizing real-time revision spreads. The underlying econometric methodology and proofs justifying the MBH approach are detailed in @kinel2026robust.

# Acknowledgements

The author acknowledges the R community for the continuous development of the foundational data manipulation and visualization tools that made this package possible.

# References
