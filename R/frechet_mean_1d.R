#' Fréchet mean of one-dimensional samples
#'
#' Computes the Wasserstein-2 Fréchet mean (barycenter) of several
#' one-dimensional empirical measures. In one dimension it has a closed form:
#' its quantile function is the average of the sample quantile functions,
#' \deqn{F_{\bar\mu}^{-1}(u) = \frac{1}{K} \sum_{k=1}^K F_{\mu_k}^{-1}(u).}
#'
#' @param samples A list of numeric vectors, one sample per measure. The
#'   samples may have different sizes.
#' @param n_grid Number of quantile levels at which to evaluate the mean.
#'   The levels are the cell midpoints \eqn{(i - 1/2)/n_{grid}}.
#'
#' @return A sorted numeric vector of length `n_grid`: an equal-weight sample
#'   representing the Fréchet mean.
#' @export
#'
#' @seealso [frechet_mean_1d_descent()] for the iterative algorithm.
#'
#' @examples
#' set.seed(1)
#' samples <- list(rnorm(100, -2), rnorm(100, 2), rnorm(100, 0, 3))
#' summary(frechet_mean_1d(samples))
frechet_mean_1d <- function(samples, n_grid = 100) {
  check_samples_1d(samples)
  n_grid <- check_count(n_grid, "n_grid")

  u <- midpoint_levels(n_grid)
  quantiles <- lapply(samples, empirical_quantile, u = u)
  Reduce(`+`, quantiles) / length(samples)
}

#' Fréchet mean of one-dimensional samples by gradient descent
#'
#' Approximates the Wasserstein-2 Fréchet mean by gradient descent on the
#' Fréchet functional \eqn{F(\gamma) = \frac{1}{2K}\sum_k W_2^2(\gamma, \mu_k)}.
#' Its Wasserstein gradient at \eqn{\gamma} is
#' \eqn{\mathrm{id} - \bar T}, where \eqn{\bar T} averages the optimal maps
#' from \eqn{\gamma} to each \eqn{\mu_k}, so each step moves the points of the
#' current estimate by
#' \deqn{\gamma \leftarrow \gamma + \eta (\bar T(\gamma) - \gamma).}
#' With step size \eqn{\eta = 1} the iteration reaches the closed-form mean of
#' [frechet_mean_1d()] after one step; smaller steps show the path of the
#' descent.
#'
#' @param samples A list of numeric vectors, one sample per measure.
#' @param init Numeric vector: the starting sample. Its size fixes the size of
#'   the estimate. Defaults to the first sample.
#' @param step Step size \eqn{\eta} in \eqn{(0, 1]}. Defaults to 1.
#' @param max_iter Maximum number of iterations.
#' @param tol Stop when the root-mean-square norm of the gradient falls below
#'   this value.
#'
#' @return A list with components
#'   \describe{
#'     \item{`mean`}{the final estimate, a sorted numeric vector;}
#'     \item{`path`}{a list of the iterates, starting with the sorted `init`;}
#'     \item{`grad_norm`}{the gradient norm before each step;}
#'     \item{`converged`}{whether the tolerance was reached.}
#'   }
#' @export
#'
#' @examples
#' set.seed(1)
#' samples <- list(rnorm(100, -2), rnorm(100, 2), rnorm(100, 0, 3))
#' fit <- frechet_mean_1d_descent(samples, step = 0.5)
#' fit$grad_norm
frechet_mean_1d_descent <- function(samples, init = samples[[1]], step = 1,
                                    max_iter = 50, tol = 1e-8) {
  check_samples_1d(samples)
  check_sample_1d(init, "init")
  check_step(step)
  max_iter <- check_count(max_iter, "max_iter")

  gamma <- sort(init)
  path <- list(gamma)
  grad_norm <- numeric(0)
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    t_avg <- Reduce(`+`, lapply(samples, transport_map_1d, x = gamma)) /
      length(samples)
    direction <- t_avg - gamma
    grad_norm <- c(grad_norm, sqrt(mean(direction^2)))
    if (grad_norm[iter] < tol) {
      converged <- TRUE
      break
    }
    gamma <- gamma + step * direction
    path[[length(path) + 1L]] <- gamma
  }

  list(mean = gamma, path = path, grad_norm = grad_norm, converged = converged)
}

# Internal checks --------------------------------------------------------

#' @noRd
check_samples_1d <- function(samples) {
  if (!is.list(samples) || length(samples) < 1L) {
    stop("`samples` must be a non-empty list of numeric vectors.", call. = FALSE)
  }
  for (k in seq_along(samples)) {
    check_sample_1d(samples[[k]], paste0("samples[[", k, "]]"))
  }
  invisible(samples)
}

#' @noRd
check_count <- function(n, name) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 ||
      n != round(n)) {
    stop("`", name, "` must be a positive whole number.", call. = FALSE)
  }
  as.integer(n)
}

#' @noRd
check_step <- function(step) {
  if (!is.numeric(step) || length(step) != 1L || !is.finite(step) ||
      step <= 0 || step > 1) {
    stop("`step` must be a single number in (0, 1].", call. = FALSE)
  }
  invisible(step)
}
