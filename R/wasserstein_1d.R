#' Compute 1D Wasserstein Distance
#'
#' Calculates the p-Wasserstein distance between two 1D empirical distributions.
#' The distance is defined as the minimal total transportation cost to reconfigure 
#' the mass of the first distribution into the second.
#'
#' @param x A numeric vector representing the source sample.
#' @param y A numeric vector representing the target sample.
#' @param p A numeric value >= 1 specifying the order of the metric (default is 2).
#'
#' @return A numeric value representing the W_p distance.
#' @export
#'
#' @examples
#' wasserstein1d(rnorm(100, 0, 1), rnorm(100, 5, 1))
wasserstein1d <- function(x, y, p = 2) {
  if (length(x) == length(y)) {
    # If samples are equal length, the optimal coupling is exactly matched quantiles
    cost <- mean(abs(sort(x) - sort(y))^p)
  } else {
    # For unequal sizes, evaluate inverse CDFs on a common probability grid
    probs <- seq(0, 1, length.out = max(length(x), length(y)))
    qx <- stats::quantile(x, probs, names = FALSE)
    qy <- stats::quantile(y, probs, names = FALSE)
    cost <- mean(abs(qx - qy)^p)
  }
  return(cost^(1/p))
}

#' 1D Optimal Transport Map
#'
#' Constructs the optimal transport map T(x) pushing distribution x to y.
#' In 1D, this corresponds to applying the target quantile function to the 
#' source cumulative distribution.
#'
#' @param x A numeric vector representing the source sample.
#' @param y A numeric vector representing the target sample.
#'
#' @return A numeric vector of the same length as x, representing mapped values.
#' @export
transport_map_1d <- function(x, y) {
  # F(x): Empirical CDF of the source
  F_x <- stats::ecdf(x)
  u <- F_x(x)
  
  # G^{-1}(u): Quantile function of the target evaluated at u
  mapped_x <- stats::quantile(y, probs = u, names = FALSE, type = 1)
  
  return(mapped_x)
}

#' Logarithmic Map in 1D Wasserstein Space
#'
#' Projects a target measure into the tangent space at a base measure.
#' The mapping is defined by the displacement of the optimal transport map 
#' from the identity.
#'
#' @param gamma Numeric vector for the base measure.
#' @param mu Numeric vector for the target measure.
#'
#' @return A numeric vector representing the tangent vector.
#' @export
log_map_1d <- function(gamma, mu) {
  t_map <- transport_map_1d(gamma, mu)
  
  # The tangent vector is the optimal map minus the identity (the base measure)
  tangent_vector <- t_map - gamma
  return(tangent_vector)
}

#' Exponential Map in 1D Wasserstein Space
#'
#' Pushes a tangent vector forward from a base measure to recover a distribution.
#' This is defined as the push-forward of the base measure under the identity 
#' plus the tangent vector.
#'
#' @param gamma Numeric vector for the base measure.
#' @param v Numeric vector representing the tangent vector.
#'
#' @return A numeric vector representing the reconstructed measure.
#' @export
exp_map_1d <- function(gamma, v) {
  return(gamma + v)
}

