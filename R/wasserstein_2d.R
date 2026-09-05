#' Compute 2D Wasserstein Distance
#'
#' Calculates the W_p distance between two 2D empirical point clouds of equal size.
#'
#' @param x A numeric matrix with 2 columns representing the source point cloud.
#' @param y A numeric matrix with 2 columns representing the target point cloud.
#' @param p A numeric value >= 1 specifying the order of the metric.
#'
#' @return A numeric value representing the W_p distance.
#' @export
wasserstein_2d <- function(x, y, p = 2) {
  # Convert matrices to 'pgrid' objects for the transport package
  px <- transport::pgrid(x)
  py <- transport::pgrid(y)
  
  # Compute the Wasserstein distance directly
  dist <- transport::wasserstein(transport::pp(x), transport::pp(y), p = p)
  return(dist)
}

#' 2D Optimal Transport Map
#'
#' Finds the 1-to-1 assignment mapping points in x to points in y.
#' Assumes x and y have the same number of points.
#'
#' @param x A numeric matrix with 2 columns.
#' @param y A numeric matrix with 2 columns.
#'
#' @return A numeric matrix representing y reordered to match x.
#' @export
transport_map_2d <- function(x, y) {
  # Compute the optimal coupling matrix (assignment)
  coupling <- transport::transport(transport::pp(x), transport::pp(y), p = 2)
  
  # Reorder the target points based on the optimal assignment indices
  # The 'transport' output contains 'from' and 'to' indices
  mapped_y <- y[coupling$to, , drop = FALSE]
  
  return(mapped_y)
}