#' Exact 1D Fréchet Mean
#'
#' Computes the closed-form Fréchet mean for a collection of 1D distributions.
#' The result is the measure whose quantile function is the arithmetic average 
#' of the sample quantile functions.
#'
#' @param samples A list of numeric vectors, each representing an empirical distribution.
#' @param n_grid Integer. The number of probability grid points to evaluate.
#'
#' @return A numeric vector representing the quantiles of the Fréchet mean.
#' @export
frechet_mean_1d_exact <- function(samples, n_grid = 100) {
  # Define a common probability grid
  probs <- seq(0, 1, length.out = n_grid)
  
  # Extract quantile functions (inverse CDFs) evaluated on the common grid
  quantiles_list <- lapply(samples, function(x) {
    stats::quantile(x, probs = probs, names = FALSE, type = 1)
  })
  
  # Average the quantile functions in the tangent/Hilbert space
  mean_quantiles <- Reduce("+", quantiles_list) / length(samples)
  
  return(mean_quantiles)
}

#' Steepest Descent for Fréchet Mean
#'
#' Iteratively approximates the Fréchet mean by averaging optimal transport maps.
#'
#' @param gamma_init Numeric vector. The initial starting measure.
#' @param target_samples A list of numeric vectors representing the target measures.
#' @param max_iter Integer. Maximum number of iterations.
#' @param tol Numeric. Tolerance for the gradient norm to determine convergence.
#'
#' @return A numeric vector representing the converged Fréchet mean.
#' @export
frechet_mean_1d_steepest_descent <- function(gamma_init, target_samples, max_iter = 20, tol = 1e-4) {
  gamma_current <- sort(gamma_init)
  N <- length(target_samples)
  
  for (iter in seq_len(max_iter)) {
    # 1. Compute optimal maps from the current estimate to all targets
    t_maps <- lapply(target_samples, function(mu) {
      transport_map_1d(gamma_current, mu)
    })
    
    # 2. Average the transport maps
    t_avg <- Reduce("+", t_maps) / N
    
    # 3. Check for convergence using the norm of the gradient
    # The gradient is the displacement: t_avg - identity (gamma_current)
    gradient <- t_avg - gamma_current
    grad_norm <- sqrt(mean(gradient^2))
    
    if (grad_norm < tol) {
      break
    }
    
    # 4. Retract back to the manifold (Exponential map step)
    gamma_current <- sort(t_avg)
  }
  
  return(gamma_current)
}