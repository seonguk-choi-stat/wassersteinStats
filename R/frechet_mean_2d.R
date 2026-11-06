#' Fréchet mean of planar point clouds
#'
#' Approximates the Wasserstein-2 Fréchet mean (barycenter) of several point
#' clouds in \eqn{R^2} by the fixed-point iteration of Álvarez-Esteban et al.
#' (2016) with free support: at each step, match the current estimate to every
#' cloud by an optimal assignment and move each point toward the average of
#' its matched points,
#' \deqn{\gamma \leftarrow \gamma + \eta (\bar T(\gamma) - \gamma).}
#' With step size \eqn{\eta = 1} this is the plain fixed-point iteration.
#' Unlike the one-dimensional case there is no closed form, and the iteration
#' may stop at a stationary point that is not the global minimizer; trying
#' several starting clouds helps.
#'
#' @param samples A list of numeric matrices with two columns, all with the
#'   same number of rows.
#' @param init Numeric matrix with two columns and that same number of rows:
#'   the starting cloud. Defaults to the first sample.
#' @param step Step size \eqn{\eta} in \eqn{(0, 1]}. Defaults to 1.
#' @param max_iter Maximum number of iterations.
#' @param tol Stop when the root-mean-square length of the step direction
#'   falls below this value.
#'
#' @return A list with components
#'   \describe{
#'     \item{`mean`}{the final estimate, a numeric matrix;}
#'     \item{`path`}{a list of the iterates, starting with `init`;}
#'     \item{`grad_norm`}{the gradient norm before each step;}
#'     \item{`converged`}{whether the tolerance was reached.}
#'   }
#' @references Álvarez-Esteban, P. C., del Barrio, E., Cuesta-Albertos, J. A.,
#'   and Matrán, C. (2016). A fixed-point approach to barycenters in
#'   Wasserstein space. \emph{Journal of Mathematical Analysis and
#'   Applications}, 441(2), 744--762.
#' @export
#'
#' @examples
#' set.seed(1)
#' base <- matrix(rnorm(60), ncol = 2)
#' samples <- list(base, sweep(base, 2, c(4, 0), "+"), sweep(base, 2, c(0, 4), "+"))
#' fit <- frechet_mean_2d(samples)
#' colMeans(fit$mean)
frechet_mean_2d <- function(samples, init = samples[[1]], step = 1,
                            max_iter = 100, tol = 1e-8) {
  if (!is.list(samples) || length(samples) < 1L) {
    stop("`samples` must be a non-empty list of two-column matrices.",
         call. = FALSE)
  }
  for (k in seq_along(samples)) {
    check_cloud_2d(samples[[k]], paste0("samples[[", k, "]]"))
  }
  check_cloud_2d(init, "init")
  if (any(vapply(samples, nrow, integer(1)) != nrow(init))) {
    stop("All samples and `init` must have the same number of points.",
         call. = FALSE)
  }
  check_step(step)
  max_iter <- check_count(max_iter, "max_iter")

  gamma <- init
  path <- list(gamma)
  grad_norm <- numeric(0)
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    t_avg <- Reduce(`+`, lapply(samples, transport_map_2d, x = gamma)) /
      length(samples)
    direction <- t_avg - gamma
    grad_norm <- c(grad_norm, sqrt(mean(rowSums(direction^2))))
    if (grad_norm[iter] < tol) {
      converged <- TRUE
      break
    }
    gamma <- gamma + step * direction
    path[[length(path) + 1L]] <- gamma
  }

  list(mean = gamma, path = path, grad_norm = grad_norm, converged = converged)
}
