#' Plot the global status of planetary boundaries
#'
#' Plot global map(s) with the status of planetary boundaries for
#' a scenario LPJmL run and derived planetary boundary statuses
#'
#' @param x  output object from calc_* with the status of the
#' control variable for one point in time, incl. pb thresholds as attribute
#'
#' @param file_name character string providing file name (including directory
#' and file extension). Defaults to NULL (plotting to screen)
#'
#' @param legend logical, specify whether a legend should be plotted
#'
#' @param projection character string defining the projection, default set to
#' "+proj=robin"
#'
#' @param ncol integer, number of columns in the plot, default set to 2
#'
#' @param grid_path character string providing the path to a grid file
#'
#' @examples
#' \dontrun{
#'  plot_status(file_name = "./my_boundary_status.png",
#'              x = calc_output
#'              legend = FALSE
#'              grid_path = "/path/to/gridfile.bin.json")
#' }
#'
#' @md
#' @export

plot_status <- function(
  x,
  file_name = NULL,
  legend = TRUE,
  projection = "+proj=robin",
  ncol = 2,
  grid_path = NULL
) {

  if (class(x[[1]]) != "control_variable") {
    stop("x elements must be of class control variable")
  }
  # TODO check for spatial scale!

  # plot settings
  if (length(x) == 1) {
    ncol <- 1
  }
  nrow <- ceiling(length(x) / ncol)
  leg_adj <- ifelse(legend, 2, 0)

  plot_nat <- to_raster(lpjml_array = array(0, length(x[[which(!is.na(x))[1]]])), # nolint
                        projection = projection,
                        grid_path = grid_path)


  if (!is.null(file_name)) {

    file_extension <- strsplit(
      file_name, split = "\\."
    )[[1]][-1] %>%
      tail(1)
    switch(file_extension,
      `png` = {
        png(file_name,
            width = 8 * ncol,
            height = 4 * nrow + leg_adj,
            units = "cm",
            res = 600,
            pointsize = 7)
      },
      `pdf` = {
        pdf(file_name,
            width = 8 * ncol / 2.54,
            height = (4 * nrow + leg_adj) / 2.54,
            pointsize = 7)
      }, {
        stop("File extension ", dQuote(file_extension), " not supported.")
      }
    )
  }

  NA_col <- c("grey92")

  # get country outlines
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

  plot_list <- list()
  for (i in seq_len(length(x))) {

    # convert lpjml vector with continuous control variable status to risk level
    plot_data <- x[[i]] %>%
      as_risk_level(type = "continuous", normalize = "increasing risk")

    # convert lpjml vector to raster with defined projection
    plotvar <- to_raster(lpjml_array = plot_data,
                         projection = projection,
                         grid_path = grid_path)

    #combine NA cells with cells in safe zone
    plot_nat[plotvar < 1] <- 1

    # prepare plotting of values > pb threshold
    plotvar_risk <- plotvar
    plotvar_risk[plotvar_risk <= 1] <- NA
    plotvar_risk[plotvar_risk > 3.5] <- 3.5

    p <- ggplot2::ggplot() +
      ggspatial::layer_spatial(data = plot_nat) +
      ggplot2::scale_fill_continuous(na.value = NA,
                                     low = NA_col,
                                     high = "#7ac4a7") +
      ggnewscale::new_scale_fill() +
      ggspatial::layer_spatial(plotvar_risk) +
      ggplot2::scale_fill_viridis_c(na.value = NA, direction = -1,
                                    option = "A") +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.1,
                                                              color = "#8d8b8b"),
                     axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     legend.position = "none") +
      ggplot2::geom_sf(data = world, fill = NA, linewidth = 0.12,
                       color = "#7e7d7d") +
      xlim(terra::ext(plotvar)[1], terra::ext(plotvar)[2]) +
      ylim(terra::ext(plotvar)[3], terra::ext(plotvar)[4]) +
      xlab(names(x)[i])

    # combine all plots in one list
    plot_list[[i]] <- p
  }

  if (legend) {
    plot_list["legend"] <- plot_legend()
  }

  gridExtra::grid.arrange(grobs = plot_list, ncol = ncol)

  if (!is.null(file_name)) dev.off()
}


# convert lpjml vector to raster and change projection to robinson
to_raster <- function(lpjml_array, projection, grid_path) {
  grid <- read_io(grid_path)$data %>% drop
  lon <- grid[, 1]
  lat <- grid[, 2]
  ra <- terra::rast(ncols = 720, nrows = 360)
  ra[terra::cellFromXY(ra, cbind(lon, lat))] <- c(lpjml_array)
  extent <- terra::ext(c(-180, 180, -53, 85))
  ra <- crop(ra, extent)
  ra <- terra::project(ra, projection)
  return(ra)
}