#' Wasserstein distance between two one-dimensional samples
#'
#' Computes the p-Wasserstein distance between the empirical measures of two
#' samples on the real line. In one dimension the optimal coupling is the
#' monotone (quantile) coupling, so
#' \deqn{W_p(\mu, \nu)^p = \int_0^1 |F_\mu^{-1}(u) - F_\nu^{-1}(u)|^p \, du.}
#' Both quantile functions are step functions, so the integral is evaluated
#' exactly on the merged breakpoints \eqn{i/n} and \eqn{j/m}. The samples may
#' have different sizes.
#'
#' @param x,y Numeric vectors: the two samples, each given equal weight per
#'   point.
#' @param p Order of the distance, a number \eqn{\ge 1}. Defaults to 2.
#'
#' @return The distance \eqn{W_p}, a non-negative number.
#' @export
#'
#' @examples
#' set.seed(1)
#' wasserstein_1d(rnorm(200), rnorm(150, mean = 2))
#'
#' # For two normal laws, W_2 equals sqrt((m1 - m2)^2 + (s1 - s2)^2)
#' u <- (1:1000 - 0.5) / 1000
#' wasserstein_1d(qnorm(u, 0, 1), qnorm(u, 2, 3))
wasserstein_1d <- function(x, y, p = 2) {
  check_sample_1d(x, "x")
  check_sample_1d(y, "y")
  check_p(p)

  n <- length(x)
  m <- length(y)
  breaks <- sort(unique(c(seq_len(n) / n, seq_len(m) / m)))
  widths <- diff(c(0, breaks))
  mids <- breaks - widths / 2

  cost <- sum(widths * abs(empirical_quantile(x, mids) -
                             empirical_quantile(y, mids))^p)
  cost^(1 / p)
}

#' Optimal transport map between two one-dimensional samples
#'
#' Evaluates the monotone transport map \eqn{T = F_\nu^{-1} \circ F_\mu} that
#' pushes the empirical measure of `x` onto that of `y`. Each point of `x` is
#' assigned the level of its rank cell, \eqn{(r - 1/2)/n}, and sent to the
#' target quantile at that level. Ties in `x` are broken by order of
#' appearance. When `x` and `y` have the same size, the result is the sorted
#' matching: the \eqn{i}-th smallest point of `x` goes to the \eqn{i}-th
#' smallest point of `y`.
#'
#' @param x Numeric vector: the source sample.
#' @param y Numeric vector: the target sample.
#'
#' @return A numeric vector of the same length and order as `x`: the image
#'   \eqn{T(x_i)} of each source point.
#' @export
#'
#' @examples
#' transport_map_1d(c(3, 1, 2), c(10, 30, 20))
transport_map_1d <- function(x, y) {
  check_sample_1d(x, "x")
  check_sample_1d(y, "y")

  u <- (rank(x, ties.method = "first") - 0.5) / length(x)
  empirical_quantile(y, u)
}

#' Logarithmic map in one-dimensional Wasserstein space
#'
#' Lifts the target measure `mu` to the tangent space at the base measure
#' `gamma`. The tangent vector is the displacement of the optimal map from
#' the identity, \eqn{\log_\gamma(\mu) = T - \mathrm{id}}, evaluated at the
#' points of `gamma`.
#'
#' @param gamma Numeric vector: sample of the base measure.
#' @param mu Numeric vector: sample of the target measure.
#'
#' @return A numeric vector of the same length and order as `gamma`: the
#'   displacement \eqn{T(\gamma_i) - \gamma_i} of each base point.
#' @export
#'
#' @examples
#' log_map_1d(c(0, 1, 2), c(1, 3, 5))
log_map_1d <- function(gamma, mu) {
  transport_map_1d(gamma, mu) - gamma
}

#' Exponential map in one-dimensional Wasserstein space
#'
#' Pushes the base measure forward along a tangent vector,
#' \eqn{\exp_\gamma(t v) = (\mathrm{id} + t v)_\# \gamma}. With
#' `v = log_map_1d(gamma, mu)`, varying `t` from 0 to 1 traces the
#' Wasserstein geodesic (displacement interpolation) from `gamma` to `mu`.
#'
#' @param gamma Numeric vector: sample of the base measure.
#' @param v Numeric vector of the same length as `gamma`: the tangent vector.
#' @param t Step along the tangent vector. Defaults to 1. For a tangent vector
#'   produced by [log_map_1d()], values in \eqn{[0, 1]} stay on the geodesic.
#'
#' @return A numeric vector: the points of the pushed-forward sample, in the
#'   order of `gamma`.
#' @export
#'
#' @examples
#' gamma <- c(0, 1, 2)
#' v <- log_map_1d(gamma, c(1, 3, 5))
#' exp_map_1d(gamma, v, t = 0.5)
exp_map_1d <- function(gamma, v, t = 1) {
  check_sample_1d(gamma, "gamma")
  check_sample_1d(v, "v")
  if (length(v) != length(gamma)) {
    stop("`v` must have the same length as `gamma`.", call. = FALSE)
  }
  if (!is.numeric(t) || length(t) != 1L || !is.finite(t)) {
    stop("`t` must be a single finite number.", call. = FALSE)
  }
  gamma + t * v
}
