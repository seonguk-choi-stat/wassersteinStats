# Internal helpers -------------------------------------------------------

#' Check a one-dimensional sample
#' @noRd
check_sample_1d <- function(x, name = "x") {
  if (!is.numeric(x) || length(x) < 1L || anyNA(x) || any(!is.finite(x))) {
    stop("`", name, "` must be a non-empty numeric vector of finite values.",
         call. = FALSE)
  }
  invisible(x)
}

#' Check a two-dimensional point cloud
#' @noRd
check_cloud_2d <- function(x, name = "x") {
  if (!is.matrix(x) || !is.numeric(x) || ncol(x) != 2L || nrow(x) < 1L ||
      anyNA(x) || any(!is.finite(x))) {
    stop("`", name, "` must be a numeric matrix with two columns and finite values.",
         call. = FALSE)
  }
  invisible(x)
}

#' Check the order of a Wasserstein distance
#' @noRd
check_p <- function(p) {
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p < 1) {
    stop("`p` must be a single number >= 1.", call. = FALSE)
  }
  invisible(p)
}

#' Quantile function of an empirical measure
#'
#' Left-continuous inverse of the empirical CDF of `x`, evaluated at levels
#' `u` in (0, 1]: Q(u) = x_(ceiling(n u)). Equivalent to
#' `stats::quantile(x, u, type = 1)` without its rounding conventions.
#' @noRd
empirical_quantile <- function(x, u) {
  xs <- sort(x)
  n <- length(xs)
  xs[pmin(n, pmax(1L, ceiling(n * u - 1e-12)))]
}

#' Midpoint levels of an equal-weight grid
#' @noRd
midpoint_levels <- function(n) {
  (seq_len(n) - 0.5) / n
}
