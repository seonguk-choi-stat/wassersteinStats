#' Module: Fréchet mean of point clouds in the plane
#'
#' Several Gaussian point clouds, their Wasserstein barycenter found by the
#' fixed-point iteration, and the optimal matching from the barycenter to
#' each cloud.
#'
#' @param id Module id.
#' @noRd
mod_cloud_2d_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 3,
      sliderInput(ns("k"), "Number of clouds", min = 2, max = 5, value = 3,
                  step = 1),
      sliderInput(ns("n"), "Points per cloud", min = 10, max = 100,
                  value = 40, step = 10),
      sliderInput(ns("spread"), "Distance between clouds", min = 1, max = 6,
                  value = 3, step = 0.5),
      sliderInput(ns("step"), "Step size", min = 0.1, max = 1, value = 0.5,
                  step = 0.1),
      selectInput(ns("matched"), "Show matching to", choices = c("None")),
      numericInput(ns("seed"), "Random seed", value = 2026, min = 1, step = 1)
    ),
    mainPanel(
      width = 9,
      tags$p(
        "In the plane the barycenter has no closed form. Each iteration",
        "matches the current estimate to every cloud by an optimal",
        "assignment and moves each point toward the average of its matched",
        "points. For Gaussian laws the barycenter is again Gaussian, with a",
        "covariance between theirs."
      ),
      sliderInput(ns("iter"), "Iteration", min = 0, max = 1, value = 0,
                  step = 1, width = "100%",
                  animate = animationOptions(interval = 500)),
      fluidRow(
        column(8, plotOutput(ns("cloud_plot"), height = "520px")),
        column(4,
               tags$h5("W2 from the shown iterate to each cloud"),
               tableOutput(ns("distance_table")),
               plotOutput(ns("grad_plot"), height = "260px"))
      )
    )
  )
}

#' @rdname mod_cloud_2d_ui
#' @noRd
mod_cloud_2d_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    clouds <- reactive({
      req(input$k, input$n, input$spread, input$seed)
      set.seed(input$seed)
      lapply(seq_len(input$k), function(i) {
        angle <- stats::runif(1, 0, pi)
        rotation <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)),
                           2, 2)
        axes <- diag(stats::runif(2, 0.3, 1.2))
        centre <- stats::runif(2, -input$spread, input$spread)
        z <- matrix(stats::rnorm(2 * input$n), ncol = 2)
        sweep(z %*% t(rotation %*% axes), 2, centre, "+")
      })
    })

    observeEvent(clouds(), {
      choices <- c("None", paste("Cloud", seq_along(clouds())))
      updateSelectInput(session, "matched", choices = choices,
                        selected = "Cloud 1")
    })

    fit <- reactive({
      req(input$step)
      frechet_mean_2d(clouds(), init = clouds()[[1]], step = input$step,
                      max_iter = 40, tol = 1e-6)
    })

    observeEvent(fit(), {
      sync_iteration_slider(session, "iter", length(fit()$path), 0L)
    })

    shown_iter <- reactive({
      it <- if (is.null(input$iter)) 0L else as.integer(input$iter)
      min(it, length(fit()$path) - 1L)
    })
    current <- reactive(fit()$path[[shown_iter() + 1L]])

    distances <- reactive(
      vapply(clouds(), wasserstein_2d, numeric(1), x = current())
    )

    output$distance_table <- renderTable({
      d <- distances()
      data.frame(Cloud = c(paste("Cloud", seq_along(d)), "Mean of W2\u00b2"),
                 Value = c(d, mean(d^2)))
    }, digits = 4)

    output$cloud_plot <- renderPlot({
      pts <- do.call(rbind, lapply(seq_along(clouds()), function(i) {
        data.frame(x = clouds()[[i]][, 1], y = clouds()[[i]][, 2],
                   cloud = paste("Cloud", i))
      }))
      est <- data.frame(x = current()[, 1], y = current()[, 2])

      p <- ggplot2::ggplot()
      selected <- if (is.null(input$matched)) "None" else input$matched
      chosen <- match(selected, paste("Cloud", seq_along(clouds())))
      if (!is.na(chosen)) {
        target <- transport_map_2d(current(), clouds()[[chosen]])
        seg <- data.frame(x = est$x, y = est$y,
                          xend = target[, 1], yend = target[, 2])
        p <- p + ggplot2::geom_segment(
          data = seg,
          ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                       yend = .data$yend),
          colour = "grey60", alpha = 0.7
        )
      }
      p +
        ggplot2::geom_point(
          data = pts,
          ggplot2::aes(x = .data$x, y = .data$y, colour = .data$cloud),
          alpha = 0.7, size = 2
        ) +
        ggplot2::geom_point(
          data = est, ggplot2::aes(x = .data$x, y = .data$y),
          shape = 17, size = 3, colour = "black"
        ) +
        ggplot2::coord_equal() +
        ggplot2::labs(
          title = sprintf("Barycenter estimate (black), iteration %d",
                          shown_iter()),
          x = NULL, y = NULL, colour = NULL
        ) +
        theme_ws()
    })

    output$grad_plot <- renderPlot({
      plot_grad_norm(fit()$grad_norm, shown_iter())
    })

    invisible(list(clouds = clouds, fit = fit, distances = distances))
  })
}
