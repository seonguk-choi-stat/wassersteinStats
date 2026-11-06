# Helpers shared by the Shiny modules ------------------------------------

#' Distribution families offered in the app
#'
#' Each family is standardized to mean 0 and variance 1, so that the `mean`
#' and `sd` inputs of the app mean what they say.
#' @noRd
distribution_families <- c("Normal", "Uniform", "Exponential", "Bimodal")

#' Draw a sample from a standardized family
#'
#' @param family One of `distribution_families`.
#' @param n Sample size.
#' @param mean,sd Mean and standard deviation of the law.
#' @return A numeric vector of length `n`.
#' @noRd
sample_distribution <- function(family, n, mean = 0, sd = 1) {
  n <- check_count(n, "n")
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0) {
    stop("`sd` must be a single positive number.", call. = FALSE)
  }
  z <- switch(
    family,
    Normal = stats::rnorm(n),
    Uniform = stats::runif(n, -sqrt(3), sqrt(3)),
    Exponential = stats::rexp(n) - 1,
    # equal mixture of N(-0.9, 0.19) and N(0.9, 0.19): mean 0, variance 1
    Bimodal = sample(c(-0.9, 0.9), n, replace = TRUE) +
      stats::rnorm(n, sd = sqrt(0.19)),
    stop("Unknown distribution family: ", family, call. = FALSE)
  )
  mean + sd * z
}

#' Inputs that describe one distribution
#' @noRd
distribution_inputs <- function(ns, prefix, label, family = "Normal",
                                mean = 0, sd = 1) {
  tagList(
    tags$h5(label),
    selectInput(ns(paste0(prefix, "_family")), "Family",
                choices = distribution_families, selected = family),
    sliderInput(ns(paste0(prefix, "_mean")), "Mean",
                min = -5, max = 5, value = mean, step = 0.1),
    sliderInput(ns(paste0(prefix, "_sd")), "Standard deviation",
                min = 0.2, max = 3, value = sd, step = 0.1)
  )
}

#' Common look for the app's plots
#' @noRd
theme_ws <- function() {
  ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank())
}

#' Gradient norms of a descent, with the shown iteration marked
#' @noRd
plot_grad_norm <- function(grad_norm, current) {
  df <- data.frame(iteration = seq_along(grad_norm) - 1L, norm = grad_norm)
  df <- df[df$norm > 0, , drop = FALSE]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$norm)) +
    ggplot2::geom_line(colour = "grey40") +
    ggplot2::geom_point(colour = "grey40") +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Iteration", y = "Gradient norm (log scale)") +
    theme_ws()
  marked <- df[df$iteration == current, , drop = FALSE]
  if (nrow(marked) > 0) {
    p <- p + ggplot2::geom_point(data = marked, colour = "darkorange", size = 4)
  }
  p
}

#' Keep an iteration slider within the length of a descent path
#' @noRd
sync_iteration_slider <- function(session, input_id, path_length, current) {
  last <- path_length - 1L
  updateSliderInput(session, input_id, max = last,
                    value = min(if (is.null(current)) last else current, last))
}
