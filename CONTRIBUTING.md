# Contributing to MacroFilters

Thank you for your interest in contributing to `MacroFilters`! This
package aims to provide high-performance, robust macroeconomic trend
extraction in R. We welcome contributions from the community, whether
they are bug reports, feature requests, or code contributions.

## 1. Reporting Bugs

If you find a bug, please open an issue on GitHub. To help us understand
and fix the problem quickly, please provide: \* A brief description of
the issue. \* A minimal reproducible example (`reprex`) using the
`reprex` package if possible. \* Your R session information (run
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)).

## 2. Suggesting Enhancements

We are always looking for ways to improve the MacroBoost Hybrid (MBH)
filter. If you have an idea for a new feature or a performance
enhancement, please open an issue to discuss it before writing any code.
This ensures your effort aligns with the package’s roadmap.

## 3. Contributing Code (Pull Requests)

If you’d like to contribute directly to the codebase, please follow
these steps:

1.  **Fork** the repository and clone it locally.
2.  **Branch:** Create a new branch for your feature or bugfix
    (`git checkout -b feature/your-feature-name`).
3.  **Develop:** Make your changes.
    - **Performance:** `MacroFilters` relies heavily on `data.table` for
      speed and memory efficiency. Please ensure new data manipulations
      follow `data.table` best practices (e.g., update by reference when
      appropriate).
    - **Style:** Try to follow the standard R tidyverse style guide for
      code readability, even though we use `data.table` under the hood.
4.  **Document:** If you add or modify a function, update the `roxygen2`
    comments above the function. Run `devtools::document()` to update
    the `NAMESPACE` and `.Rd` files.
5.  **Test:** Add unit tests for your new feature in the
    `tests/testthat/` directory. Run `devtools::test()` to ensure all
    tests (old and new) pass.
6.  **Check:** Run `devtools::check()` locally to ensure there are no
    warnings, errors, or notes.
7.  **Commit & Push:** Commit your changes with a clear, descriptive
    message and push them to your fork.
8.  **Pull Request:** Open a Pull Request against the `main` branch of
    this repository. Link any relevant issues in the PR description.

## 4. Code of Conduct

Please note that this project is released with a Contributor Code of
Conduct. By participating in this project you agree to abide by its
terms. We expect all contributors to maintain a respectful and
constructive environment.
