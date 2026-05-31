# ── Internal utilities (not exported) ─────────────────────────────────────────

#' Fast HP cycle for internal calibration
#'
#' Solves (I + lambda * D'D) t = y via sparse Cholesky and returns only the
#' cycle (y - trend). No input validation, no class wrapping, no timing.
#' `stats::frequency(x)` drives the Ravn-Uhlig lambda so the threshold is
#' always on the correct residual scale.
#'
#' @param x Numeric vector or `ts` object.
#' @return Numeric vector of cyclical residuals, same length as `x`.
#' @keywords internal
#' @noRd
.hp_fast <- function(x) {
  n <- length(x)
  if (n < 3L) return(rep(0, n))          # bandSparse needs n-2 >= 1 rows

  freq   <- stats::frequency(x)          # 1 for numeric/annual, 4/12 for ts
  y      <- as.double(x)
  lambda <- 6.25 * freq^4                # Ravn-Uhlig rule

  D <- Matrix::bandSparse(
    n - 2L, n,
    k         = c(0L, 1L, 2L),
    diagonals = list(rep(1, n - 2L), rep(-2, n - 2L), rep(1, n - 2L))
  )
  Q <- Matrix::Diagonal(n) + lambda * Matrix::crossprod(D)
  y - as.numeric(Matrix::solve(Q, y))
}
