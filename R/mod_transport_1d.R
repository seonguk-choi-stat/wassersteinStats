#' Module: optimal transport and geodesics on the line
#'
#' Two distributions, the optimal (monotone) transport map between them, and
#' the Wasserstein geodesic that morphs one into the other as `t` runs from 0
#' to 1.
#'
#' @param id Module id.
#' @noRd
mod_transport_1d_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 3,
      distribution_inputs(ns, "src", "Source", family = "Normal",
                          mean = -2, sd = 1),
      distribution_inputs(ns, "tgt", "Target", family = "Bimodal",
                          mean = 2, sd = 1.5),
      tags$hr(),
      sliderInput(ns("n"), "Points per sample", min = 50, max = 1000,
                  value = 300, step = 50),
      radioButtons(ns("p"), "Order p of the distance", choices = c("1", "2"),
                   selected = "2", inline = TRUE),
      numericInput(ns("seed"), "Random seed", value = 2026, min = 1, step = 1)
    ),
    mainPanel(
      width = 9,
      tags$p(
        "In one dimension the optimal transport map is monotone: the i-th",
        "smallest source point goes to the i-th smallest target point.",
        "Moving every point a fraction t of the way along its arrow traces",
        "the Wasserstein geodesic between the two laws."
      ),
      sliderInput(ns("t"), "Time t along the geodesic", min = 0, max = 1,
                  value = 0.5, step = 0.02, width = "100%",
                  animate = animationOptions(interval = 150, loop = TRUE)),
      tags$h4(textOutput(ns("distance"))),
      textOutput(ns("geodesic_distance")),
      fluidRow(
        column(6, plotOutput(ns("density_plot"), height = "340px")),
        column(6, plotOutput(ns("matching_plot"), height = "340px"))
      ),
      plotOutput(ns("map_plot"), height = "320px")
    )
  )
}

#' @rdname mod_transport_1d_ui
#' @noRd
mod_transport_1d_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    samples <- reactive({
      req(input$src_family, input$tgt_family, input$n, input$seed)
      set.seed(input$seed)
      list(
        x = sort(sample_distribution(input$src_family, input$n,
                                     input$src_mean, input$src_sd)),
        y = sort(sample_distribution(input$tgt_family, input$n,
                                     input$tgt_mean, input$tgt_sd))
      )
    })

    p_order <- reactive(as.numeric(req(input$p)))
    tangent <- reactive(log_map_1d(samples()$x, samples()$y))
    interpolated <- reactive({
      req(!is.null(input$t))
      exp_map_1d(samples()$x, tangent(), input$t)
    })
    distance <- reactive(wasserstein_1d(samples()$x, samples()$y,
                                        p = p_order()))

    output$distance <- renderText(
      sprintf("W%s(source, target) = %.4f", input$p, distance())
    )
    output$geodesic_distance <- renderText(
      sprintf(paste("Distance travelled at t = %.2f: W%s(source, geodesic)",
                    "= %.4f = t \u00d7 W%s(source, target)"),
              input$t, input$p,
              wasserstein_1d(samples()$x, interpolated(), p = p_order()),
              input$p)
    )

    output$density_plot <- renderPlot({
      n <- length(samples()$x)
      labels <- c("Source", "Target", sprintf("Geodesic at t = %.2f", input$t))
      df <- data.frame(
        value = c(samples()$x, samples()$y, interpolated()),
        law = factor(rep(labels, each = n), levels = labels)
      )
      ggplot2::ggplot(df, ggplot2::aes(x = .data$value, fill = .data$law,
                                       colour = .data$law)) +
        ggplot2::geom_density(alpha = 0.25, linewidth = 0.8) +
        ggplot2::scale_fill_manual(values = c("steelblue", "firebrick", "black"),
                                   name = NULL) +
        ggplot2::scale_colour_manual(values = c("steelblue", "firebrick", "black"),
                                     name = NULL) +
        ggplot2::labs(title = "Densities", x = "Value", y = "Density") +
        theme_ws()
    })

    output$matching_plot <- renderPlot({
      x <- samples()$x
      n <- length(x)
      idx <- unique(round(seq(1, n, length.out = min(n, 60))))
      seg <- data.frame(from = x[idx], to = samples()$y[idx],
                        now = interpolated()[idx])
      ggplot2::ggplot(seg) +
        ggplot2::geom_segment(
          ggplot2::aes(x = .data$from, xend = .data$to, y = 1, yend = 0),
          colour = "grey60", alpha = 0.6
        ) +
        ggplot2::geom_point(ggplot2::aes(x = .data$from, y = 1),
                            colour = "steelblue", size = 2) +
        ggplot2::geom_point(ggplot2::aes(x = .data$to, y = 0),
                            colour = "firebrick", size = 2) +
        ggplot2::geom_point(ggplot2::aes(x = .data$now, y = 1 - input$t),
                            colour = "black", size = 2.5) +
        ggplot2::scale_y_continuous(breaks = c(0, 1),
                                    labels = c("Target", "Source")) +
        ggplot2::labs(title = "Optimal matching (every few points shown)",
                      x = "Value", y = NULL) +
        theme_ws()
    })

    output$map_plot <- renderPlot({
      df <- data.frame(x = samples()$x, tx = samples()$y)
      ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$tx)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             colour = "grey50") +
        ggplot2::geom_step(colour = "black", linewidth = 0.8) +
        ggplot2::labs(title = "Transport map T(x) (dashed: identity)",
                      x = "Source point x", y = "T(x)") +
        theme_ws()
    })

    invisible(list(samples = samples, distance = distance,
                   interpolated = interpolated))
  })
}
