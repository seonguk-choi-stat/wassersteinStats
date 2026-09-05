#' Steepest Descent for 2D Fréchet Mean
#'
#' Iteratively approximates the Fréchet mean of multiple 2D point clouds.
#'
#' @param gamma_init A numeric matrix (n x 2). The initial starting point cloud.
#' @param target_samples A list of numeric matrices (n x 2). 
#' @param max_iter Integer. Maximum number of iterations.
#' @param tol Numeric. Tolerance for the gradient norm.
#'
#' @return A numeric matrix representing the converged Fréchet mean point cloud.
#' @export
frechet_mean_2d_steepest <- function(gamma_init, target_samples, max_iter = 20, tol = 1e-4) {
  gamma_current <- gamma_init
  N <- length(target_samples)
  
  for (iter in seq_len(max_iter)) {
    
    # 1. Compute 2D optimal maps from current estimate to all targets
    mapped_targets <- lapply(target_samples, function(mu) {
      transport_map_2d(gamma_current, mu)
    })
    
    # 2. Average the coordinates of the mapped point clouds
    t_avg <- Reduce("+", mapped_targets) / N
    
    # 3. Compute gradient norm (mean squared displacement)
    gradient <- t_avg - gamma_current
    grad_norm <- sqrt(mean(rowSums(gradient^2)))
    
    if (grad_norm < tol) {
      break
    }
    
    # 4. Update the barycenter estimate
    gamma_current <- t_avg
  }
  
  return(gamma_current)
}