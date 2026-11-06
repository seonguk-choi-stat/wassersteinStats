#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    golem_add_external_resources(),
    navbarPage(
      title = "wassersteinStats",
      id = "tabs",
      tabPanel("1D Transport", mod_transport_1d_ui("transport_1d")),
      tabPanel("Fr\u00e9chet Mean (1D)", mod_frechet_1d_ui("frechet_1d")),
      tabPanel("2D Point Clouds", mod_cloud_2d_ui("cloud_2d")),
      tabPanel(
        "About",
        fluidRow(column(
          8, offset = 2,
          tags$h3("Statistics in Wasserstein space"),
          tags$p(
            "This app accompanies the R package wassersteinStats. Each tab",
            "illustrates one object of optimal transport: the transport map",
            "and geodesics between two distributions on the line, the",
            "Fr\u00e9chet mean (barycenter) of several distributions on the",
            "line, and the barycenter of point clouds in the plane."
          ),
          tags$p(
            "All computations call the package functions",
            tags$code("wasserstein_1d()"), ",", tags$code("log_map_1d()"), ",",
            tags$code("exp_map_1d()"), ",", tags$code("frechet_mean_1d()"), ",",
            tags$code("frechet_mean_1d_descent()"), ",",
            tags$code("wasserstein_2d()"), ", and",
            tags$code("frechet_mean_2d()"), "."
          ),
          tags$p(
            "Reference: V. M. Panaretos and Y. Zemel,",
            tags$em("An Invitation to Statistics in Wasserstein Space"),
            "(Springer, 2020)."
          ),
          tags$p(
            "Source code:",
            tags$a(href = "https://github.com/seonguk-choi-stat/wassersteinStats",
                   "github.com/seonguk-choi-stat/wassersteinStats")
          )
        ))
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "wassersteinStats"
    )
  )
}
