#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  mod_transport_1d_server("transport_1d")
  mod_frechet_1d_server("frechet_1d")
  mod_cloud_2d_server("cloud_2d")
}
