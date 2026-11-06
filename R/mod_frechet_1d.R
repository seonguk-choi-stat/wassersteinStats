#' Module: Fréchet mean of distributions on the line
#'
#' Several distributions, their Wasserstein barycenter, the ordinary
#' (Euclidean) average of their densities, and the gradient descent that
#' reaches the barycenter.
#'
#' @param id Module id.
#' @noRd
mod_frechet_1d_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput(ns("scenario"), "Distributions",
                  choices = c("Shifted normals", "Shifted bimodals",
                              "Mixed families")),
      sliderInput(ns("k"), "Number of distributions", min = 2, max = 8,
                  value = 4, step = 1),
      sliderInput(ns("n"), "Points per sample", min = 50, max = 500,
                  value = 200, step = 50),
      tags$hr(),
      sliderInput(ns("step"), "Step size of gradient descent", min = 0.1,
                  max = 1, value = 0.3, step = 0.1),
      numericInput(ns("seed"), "Random seed", value = 2026, min = 1, step = 1)
    ),
    mainPanel(
      width = 9,
      tags$p(
        "The Wasserstein barycenter averages quantile functions, so it keeps",
        "the shape the distributions share. Averaging the densities instead",
        "gives a mixture with one bump per distribution. Gradient descent",
        "starts from the first distribution and moves its points toward the",
        "average of their transport maps."
      ),
      sliderInput(ns("iter"), "Iteration of gradient descent", min = 0,
                  max = 1, value = 0, step = 1, width = "100%",
                  animate = animationOptions(interval = 500)),
      textOutput(ns("status")),
      plotOutput(ns("density_plot"), height = "420px"),
      plotOutput(ns("grad_plot"), height = "260px")
    )
  )
}

#' @rdname mod_frechet_1d_ui
#' @noRd
mod_frechet_1d_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    samples <- reactive({
      req(input$scenario, input$k, input$n, input$seed)
      set.seed(input$seed)
      lapply(seq_len(input$k), function(i) {
        family <- switch(input$scenario,
                         "Shifted normals" = "Normal",
                         "Shifted bimodals" = "Bimodal",
                         "Mixed families" = sample(distribution_families, 1))
        sample_distribution(family, input$n,
                            mean = stats::runif(1, -4, 4),
                            sd = stats::runif(1, 0.5, 1.5))
      })
    })

    barycenter <- reactive(frechet_mean_1d(samples(), n_grid = input$n))

    fit <- reactive({
      req(input$step)
      frechet_mean_1d_descent(samples(), init = samples()[[1]],
                              step = input$step, max_iter = 60, tol = 1e-6)
    })

    observeEvent(fit(), {
      sync_iteration_slider(session, "iter", length(fit()$path), 0L)
    })

    shown_iter <- reactive({
      it <- if (is.null(input$iter)) 0L else as.integer(input$iter)
      min(it, length(fit()$path) - 1L)
    })
    current <- reactive(fit()$path[[shown_iter() + 1L]])

    output$status <- renderText({
      sprintf(paste("Showing iteration %d of %d.",
                    "W2(iterate, barycenter) = %.4f."),
              shown_iter(), length(fit()$path) - 1L,
              wasserstein_1d(current(), barycenter()))
    })

    output$density_plot <- renderPlot({
      inputs <- do.call(rbind, lapply(seq_along(samples()), function(i) {
        data.frame(value = samples()[[i]], group = paste("Input", i))
      }))
      labels <- c("Wasserstein barycenter", "Average of densities",
                  sprintf("Gradient descent, iteration %d", shown_iter()))
      summaries <- rbind(
        data.frame(value = barycenter(), group = labels[1]),
        data.frame(value = unlist(samples()), group = labels[2]),
        data.frame(value = current(), group = labels[3])
      )
      summaries$group <- factor(summaries$group, levels = labels)

      ggplot2::ggplot() +
        ggplot2::geom_density(
          data = inputs,
          ggplot2::aes(x = .data$value, group = .data$group),
          colour = "grey75", linewidth = 0.5
        ) +
        ggplot2::geom_density(
          data = summaries,
          ggplot2::aes(x = .data$value, colour = .data$group,
                       linetype = .data$group),
          linewidth = 1.1
        ) +
        ggplot2::scale_colour_manual(
          values = c("steelblue", "darkorange", "black"), name = NULL
        ) +
        ggplot2::scale_linetype_manual(
          values = c("solid", "solid", "dashed"), name = NULL
        ) +
        ggplot2::labs(title = "Inputs (grey) and their averages",
                      x = "Value", y = "Density") +
        theme_ws()
    })

    output$grad_plot <- renderPlot({
      plot_grad_norm(fit()$grad_norm, shown_iter())
    })

    invisible(list(samples = samples, barycenter = barycenter, fit = fit))
  })
}
