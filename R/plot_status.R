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
#' @param stylized Logical. If `spatial_scale == "global"`, the function will
#'  plot the status of the planetary boundaries using a stylized plot.
#'
#' @examples
#' \dontrun{
#' plot_status(
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
  ...
) {

  if (!is.null(filename)) {
    file_extension <- file_ext(filename)
    if (!file_extension %in% c("pdf", "png")) {
      stop("File extension ", dQuote(file_extension), " not supported.")
    }
  }

  if (attributes(x[[1]])$spatial_scale == "global") {
    if (stylized) {
      plot_status_stylized(x, filename, add_legend, ...)
    } else {
      plot_status_global(x, filename, add_legend, ...)
    }
  } else {
    plot_status_maps(x, filename, add_legend, ...)
  }
}

file_ext <- function (x) {
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}