#' Plot the status of planetary boundaries
#'
#' Plot the status of planetary boundaries. The function takes the output from
#' calc_status and plots the status of the planetary boundaries depending on
#' the spatial scale applying different forms of plotting.
#'
#' @param x  output object from calc_* with the status of the
#' control variable for one point in time, incl. pb thresholds as attribute
#'
#' @param filename character string providing file name (including directory
#' and file extension). Defaults to NULL (plotting to screen)
#'
#' @param add_legend logical, specify whether a legend should be plotted
#'
#' @param stylized Logical. If `spatial_scale == "global"`, the function will
#'  plot the status of the planetary boundaries using a stylized plot.
#'
#' @param ... additional arguments passed to the plotting functions, see also
#' [`plot_status_global`], [`plot_status_maps`] and [`plot_status_stylized`]
#'
#' @examples
#' \dontrun{
#' pb_status <- calc_status(
#'   boundary = c("lsc", "biosphere", "bluewater", "greenwater", "nitrogen"),
#'   config_scenario = "./config_lu_1500_2016.json",
#'   config_reference = "./config_pnv_1500_2016.json",
#'   time_span_scenario = as.character(1986:2016),
#'   time_span_reference = as.character(1986:2016),
#'   spatial_scale = "global",
#'   approach = list(
#'     bluewater = "porkka2024",
#'     nitrogen = "schulte_uebbing2022"
#'   ),
#'   savanna_proxy = list(vegc = 7500),
#'   nyear_window = 1,
#'   path_baseline = "./pnv_1500_2016/",
#' )
#'
#' plot_status(
#'   x = pb_status,
#'   filename = "status.png",
#'   add_legend = TRUE,
#'   stylized = TRUE
#' )
#' }
#'
#' @md
#' @export
plot_status <- function(
    x,
    filename = NULL,
    add_legend = TRUE,
    stylized = FALSE,
    ...) {

  if (attributes(x[[1]])$spatial_scale == "global") {
    if (stylized) {
      plot <- plot_status_stylized(x, filename, add_legend, ...) # nolint:object_usage_linter
      if (is.null(filename)) {
        return(plot)
      }
    } else {
      plot <- plot_status_global(x, filename, ...) # nolint:object_usage_linter
      if (add_legend == TRUE) {
        print("Note: save legend seperately based on plot_legend function")
      }
      if (is.null(filename)) {
        return(plot)
      }
    }
  } else {
    plot <- plot_status_maps(x, filename, ...) # nolint:object_usage_linter
    if (add_legend == TRUE) {
      print("Note: save legend seperately based on plot_legend function")
    }
    if (is.null(filename)) {
      return(plot)
    }
  }
}
