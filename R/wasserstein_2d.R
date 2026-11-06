#' Optimal assignment between two planar point clouds
#'
#' Solves the optimal transport problem between two point clouds of the same
#' size with equal weight on each point and cost \eqn{|x_i - y_j|^p}. For
#' equal-size uniform clouds an optimal plan is a permutation, so the result
#' is a one-to-one assignment. The problem is solved exactly with the network
#' flow solver of the \pkg{transport} package.
#'
#' @param x,y Numeric matrices with two columns and the same number of rows.
#' @param p Order of the cost, a number \eqn{\ge 1}. Defaults to 2.
#'
#' @return An integer vector `sigma` with `y[sigma[i], ]` the point assigned
#'   to `x[i, ]`.
#' @noRd
optimal_assignment_2d <- function(x, y, p = 2) {
  check_cloud_2d(x, "x")
  check_cloud_2d(y, "y")
  check_p(p)
  if (nrow(x) != nrow(y)) {
    stop("`x` and `y` must have the same number of points.", call. = FALSE)
  }
  if (nrow(x) == 1L) {
    return(1L)
  }

  plan <- transport::transport(transport::pp(x), transport::pp(y),
                               p = p, method = "networkflow")
  if (anyDuplicated(plan$from) || nrow(plan) != nrow(x)) {
    stop("The transport solver did not return a one-to-one assignment.",
         call. = FALSE)
  }
  sigma <- integer(nrow(x))
  sigma[plan$from] <- plan$to
  sigma
}

#' Wasserstein distance between two planar point clouds
#'
#' Computes the p-Wasserstein distance between the empirical measures of two
#' point clouds in \eqn{R^2} of the same size, each point carrying equal
#' weight:
#' \deqn{W_p = \Big(\min_\sigma \frac{1}{n}\sum_{i=1}^n
#'   |x_i - y_{\sigma(i)}|^p\Big)^{1/p}.}
#'
#' @param x,y Numeric matrices with two columns and the same number of rows.
#' @param p Order of the distance, a number \eqn{\ge 1}. Defaults to 2.
#'
#' @return The distance \eqn{W_p}, a non-negative number.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(100), ncol = 2)
#' wasserstein_2d(x, sweep(x, 2, c(3, 4), "+"))  # translation by (3, 4): 5
wasserstein_2d <- function(x, y, p = 2) {
  sigma <- optimal_assignment_2d(x, y, p)
  dist <- sqrt(rowSums((x - y[sigma, , drop = FALSE])^2))
  mean(dist^p)^(1 / p)
}

#' Optimal transport map between two planar point clouds
#'
#' Returns the optimal (\eqn{W_2}) assignment of the points of `y` to the
#' points of `x`.
#'
#' @param x,y Numeric matrices with two columns and the same number of rows.
#'
#' @return A numeric matrix with the same shape as `x`: row `i` is the point
#'   of `y` matched to `x[i, ]`.
#' @export
#'
#' @examples
#' x <- cbind(c(0, 1), c(0, 0))
#' y <- cbind(c(1, 0), c(1, 1))
#' transport_map_2d(x, y)
transport_map_2d <- function(x, y) {
  sigma <- optimal_assignment_2d(x, y, p = 2)
  y[sigma, , drop = FALSE]
}
