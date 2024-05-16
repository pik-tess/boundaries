#' Plot the global status of planetary boundaries
#'
#' Plot global map(s) with the status of planetary boundaries for
#' a scenario LPJmL run and derived planetary boundary statuses. Legend can be
#' plotted seperately based on the plot_legend() function
#'
#' @param x  output object from calc_* with the status of the
#' control variable for one point in time, incl. pb thresholds as attribute
#'
#' @param filename character string providing file name (including directory
#' and file extension). Defaults to NULL (plotting to screen and returning
#' the ggplot object)
#'
#' @param risk_level logical, specify whether the status should be plotted as
#' risk level. Default set to TRUE.
#'
#' @param high_risk numeric, specify the normalized PB status value for the
#' upper end of the risk color scale (default 3.5). Only relevant if risk_level
#' is set to TRUE
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
#'  plot_status_maps(
#'   filename = "./my_boundary_status.png",
#'   x = calc_output
#'   grid_path = "/path/to/gridfile.bin.json"
#' )
#' }
#'
#' @md
#' @export

plot_status_maps <- function(
  x,
  filename = NULL,
  risk_level = TRUE,
  projection = "+proj=robin",
  high_risk = 3.5,
  ncol = 2,
  grid_path = NULL
) {

  if (class(x[[1]]) != "control_variable") {
    stop("x elements must be of class control variable")
  }

  # create empty raster for plotting
  plot_nat <- to_raster(lpjml_array = array(0, length(x[[which(!is.na(x))[1]]])), # nolint
                        projection = projection,
                        grid_path = grid_path)

  # define color for cells where PB definition does not apply
  na_col <- c("#ebebeb")

  # get country outlines
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

  plot_list <- list()

  # define order or PBs to be plotted
  pbs_ordered <- factor(names(x),
    levels = c("lsc", "biosphere", "bluewater", "greenwater", "nitrogen")
  ) %>%
    sort(decreasing = FALSE) %>%
    as.character()

  for (i in pbs_ordered) {
    if (risk_level) {
      # convert lpjml vector with continuous control variable status to risk
      # level
      plot_data <- x[[i]] %>%
        as_risk_level(type = "continuous", normalize = "increasing risk")
    } else {
      plot_data <- x[[i]]
      plot_data[plot_data > quantile(plot_data, 0.95, na.rm = TRUE)] <-
        quantile(plot_data, 0.95,  na.rm = TRUE)
      plot_data[plot_data < quantile(plot_data, 0.05, na.rm = TRUE)] <-
        quantile(plot_data, 0.05,  na.rm = TRUE)

      legend_title <- paste0(attr(x[[i]], "control_variable"),
                          " (", attr(x[[i]], "unit"), ")")
    }

    # convert lpjml vector to raster with defined projection
    plotvar <- to_raster(lpjml_array = plot_data,
                         projection = projection,
                         grid_path = grid_path)

    if (risk_level) {
      #combine NA cells with cells in safe zone
      plot_nat[plotvar < 1] <- 1
      levels(plot_nat) <- data.frame(id = c(0:1),
                                     name = c("no value", "safe"))
      terra::is.factor(plot_nat)

      # define color map
      coltb <- data.frame(value = c(0, 1), col = c(na_col, green))
      terra::coltab(plot_nat) <- coltb

      # prepare plotting of values > pb threshold
      plotvar_risk <- plotvar
      plotvar_risk[plotvar_risk <= 1] <- NA
      plotvar_risk[plotvar_risk > high_risk] <- high_risk

      # define viridis color scale end value, depending on max value
      max_value <- terra::minmax(plotvar_risk)[2]
      end_value <- max_value / high_risk
      if (end_value > 1) {
        end_value <- 1
      }
      end_value <- 0.9 * end_value #(0.9 = end of viridis scale)

      p <- ggplot2::ggplot() +
        tidyterra::geom_spatraster(data = plot_nat) +
        ggnewscale::new_scale_fill() +
        tidyterra::geom_spatraster(data = plotvar_risk) +
        ggplot2::scale_fill_viridis_c(na.value = NA, direction = -1,
                                      option = "inferno", begin = 1 - end_value,
                                      end = 0.95) +
        ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
                       panel.grid.major = ggplot2::element_line(linewidth = 0.1,
                                                                color = "#8d8b8b"), # nolint:line_length_linter
                       axis.text = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       legend.position = "none",
                       axis.title.x = ggplot2::element_text(size = 9)) +
        ggplot2::geom_sf(data = world, fill = NA, linewidth = 0.12,
                         color = "#7e7d7d") +
        ggplot2::xlim(terra::ext(plotvar)[1], terra::ext(plotvar)[2]) +
        ggplot2::ylim(terra::ext(plotvar)[3], terra::ext(plotvar)[4]) +
        ggplot2::xlab(attr(x[[i]], "long_name"))

    } else {
      # plotting of continous control variable without conversion to risk level
      levels(plot_nat) <- data.frame(id = c(0),
                                     name = c("no value"))
      terra::is.factor(plot_nat)

      # define color map
      coltb <- data.frame(value = c(0), col = c(na_col))
      terra::coltab(plot_nat) <- coltb

      p <- ggplot2::ggplot() +
        tidyterra::geom_spatraster(data = plot_nat, show.legend = FALSE) +
        ggnewscale::new_scale_fill() +
        tidyterra::geom_spatraster(data = plotvar) +
        ggplot2::scale_fill_viridis_c(na.value = NA, name = legend_title,
                                      option = "D", direction = -1, end = 0.9) +
        ggplot2::guides(fill = ggplot2::guide_colourbar(title.position = "top",
                                                        title.hjust = 0.5,
                                                        barwidth = 11,
                                                        barheight = 0.5)) +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill = "#ffffff"),
          panel.grid.major = ggplot2::element_line(linewidth = 0.1,
                                                   color = "#8d8b8b"),
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          legend.position = "bottom",
          legend.title = ggplot2::element_text(size = 6),
          legend.text = ggplot2::element_text(size = 6),
          legend.margin = ggplot2::margin(0, 0, 0, 0),
          legend.box.margin = ggplot2::margin(-10, -10, -10, -10),
          plot.title = ggplot2::element_text(hjust = 0.5, size = 8),
          plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
          panel.spacing = ggplot2::unit(0, "lines")
        ) +
        ggplot2::geom_sf(data = world, fill = NA, linewidth = 0.12,
                         color = "#7e7d7d") +
        ggplot2::xlim(terra::ext(plotvar)[1], terra::ext(plotvar)[2]) +
        ggplot2::ylim(terra::ext(plotvar)[3], terra::ext(plotvar)[4]) +
        ggplot2::ggtitle(attr(x[[i]], "long_name"))
    }
    # combine all plots in one list
    plot_list[[i]] <- p
  }

  # define nrows and ncols for arranging maps
  if (length(x) == 1) {
    ncol <- 1
  }
  nrow <- ceiling(length(x) / ncol)
  if (!risk_level) {
    nrow <- nrow + 1
  }

  # arrange maps
  plot <- ggpubr::as_ggplot(gridExtra::arrangeGrob(grobs = plot_list,
                                                   ncol = ncol)) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "#ffffff", color = NA)
    )

  if (!is.null(filename)) {
    ggplot2::ggsave(
      filename,
      plot,
      width = 8 * ncol,
      height = 4 * nrow,
      dpi = 600,
      units = "cm",
      pointsize = 7
    )
  } else {
    # plot maps to screen
    print(plot)
    return(plot)
  }
}


# convert lpjml vector to raster and change projection to robinson
to_raster <- function(lpjml_array, projection, grid_path) {
  grid <- lpjmlkit::read_io(grid_path)$data %>% drop
  lon <- grid[, 1]
  lat <- grid[, 2]
  ra <- terra::rast(ncols = 720, nrows = 360)
  ra[terra::cellFromXY(ra, cbind(lon, lat))] <- c(lpjml_array)
  extent <- terra::ext(c(-180, 180, -53, 85))
  ra <- terra::crop(ra, extent)
  ra <- terra::project(ra, projection, method = "near")
  return(ra)
}
