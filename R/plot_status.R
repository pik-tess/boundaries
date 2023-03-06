#' Plot the global status of planetary boundaries
#'
#' Plot global map(s) with the status of 1-4 planetary boundaries for
#' a scenario LPJmL run and derived planetary boundary statuses
#'
#' @param file_name character string providing file name (including directory
#' and file extension). Defaults to NULL (plotting to screen)
#'
#' @param status_data list with one to four LPJmL vectors (length: ncells) with
#' the outputs from the calc_* functions for one point in time.
#' Given element names within the list will be displayed underneath the
#' respective map. Order in the list determines the plotting order (by row).
#'
#' @param colors definition of colors for plotting of the safe, increasing
#' risk and high risk zone (3 elements, do not change the order). Colors used
#' for previous plots: "#008B0099" (green), "#FFEB00", "brown3"
#'
#' @param legend logical, specify whether a legend should be plotted
#'
#' @param to_robinson logical to define if robinson projection should be used
#' for plotting
#'
#' @param bg_col character, specify background possible (`NA` for transparent)
#'
#' @examples
#' \dontrun{
#'  plot_status(file_name = "./my_boundary_status.png",
#'                 status_data = list("land system change" = lsc_status,
#'                                    "nitrogen" = nitrogen_status,
#'                                    "bluewater" = water_status),
#'  legend = FALSE
#'  bg_col = NA)
#' }
#'
#' @md
#' @export

plot_status <- function(file_name = NULL,
                        status_data = list("biosphere" = NA,
                                           "lsc" = NA,
                                           "bluewater" = NA,
                                           "nitrogen" = NA),
                        colors = c("safe zone" = "#74bca9e1",
                                   "increasing risk" = "#f6ee0f",
                                   "high risk" = "#e23a50"),
                        bg_col = "white",
                        to_robinson = TRUE,
                        legend = TRUE) {
  # checking
  if (length(status_data) == 0 || length(status_data) > 4) {
    stop(paste0("Number of elements in status data (", length(status_data),
            ") is out of scope. 1 - 4 PB status maps can be",
            " plotted."))
  }
  if (length(colors) != 3) {
    stop(paste0("Length of color vector (", length(colors), ") incorrect. ",
                 "Three colors have to be provided"))
  }

  # load required data: bbox, countries
  lpjml_extent <- c(-180, 180, -60, 85)

  bounding_box <- system.file("extdata", "ne_110m_wgs84_bounding_box.shp",
                              package = "boundaries") %>%
      rgdal::readOGR(layer = "ne_110m_wgs84_bounding_box", verbose = FALSE) %>%
      { if(to_robinson) sp::spTransform(., sp::CRS("+proj=robin")) else . } # nolint

  countries <- system.file("extdata", "ne_110m_admin_0_countries.shp",
                              package = "boundaries") %>%
      rgdal::readOGR(layer = "ne_110m_admin_0_countries", verbose = FALSE) %>%
      crop(., lpjml_extent) %>%
      { if(to_robinson) sp::spTransform(., CRS("+proj=robin")) else . } # nolint

  pb_names <- names(status_data)

  not_na <- !is.na(status_data)

  plot_nat <- to_raster(lpjml_array = array(0, length(status_data[[which(not_na)[1]]])), # nolint
                        boundary_box = bounding_box,
                        ext = lpjml_extent,
                        to_robinson = to_robinson)

  # plot settings
  if (length(status_data) == 1) {
    n_row <- 1
    n_col <- 1
  } else if (length(status_data) == 2) {
    n_row <- 1
    n_col <- 2
  } else {
    n_row <- 2
    n_col <- 2
  }

  if (legend) {
    lfrac <- 0.1
    leg_adj <- 1
  } else {
    lfrac <- 0
    leg_adj <- 0
  }
  fig_params <- definefig(n_row, n_col, lfrac)
  if (!is.null(file_name)) {
    file_extension <- strsplit(file_name, split = "\\.")[[1]][-1]
    switch(file_extension,
      `png` = {
        png(file_name,
            width = 8 * n_col,
            height = 4 * n_row + leg_adj,
            units = "cm",
            res = 600,
            pointsize = 7)
      },
      `pdf` = {
        pdf(file_name,
            width = 8 * n_col / 2.54,
            height = (4 * n_row + leg_adj) / 2.54,
            pointsize = 7)
      }, {
        stop("File extension ", dQuote(file_extension), " not supported.")
      }
    )
  }
  textcex <- 0.5 + 0.125 * (n_row + n_col)
  par(mar = rep(0, 4), xpd = TRUE, bg = bg_col)
  brk <- c(-1:4)
  cols <- c("grey92", colors, "darkgrey")

  for (i in seq_len(length(status_data))) {
    if (i == 1) {
      par(fig = fig_params[i, ])
    } else {
      par(fig = fig_params[i, ], new = TRUE)
    }
    # convert lpjml vector to raster with robinson projection
    plotvar <- to_raster(lpjml_array = unlist(status_data[i]),
                         boundary_box = bounding_box,
                         ext = lpjml_extent,
                         to_robinson = to_robinson)

    image(plot_nat, zlim = c(-1, 0), asp = 1, xaxt = "n", yaxt = "n",
          xlab = "", ylab = "", col = "grey90", lwd = 0.1, bty = "n")
    image(plotvar, asp = 1, xaxt = "n",
          yaxt = "n", xlab = "", ylab = "", col = cols, breaks = brk,
          lwd = 0.1, bty = "n", add = TRUE)
    raster::plot(countries, add = TRUE, lwd = 0.3, border = "#33333366",
         usePolypath = FALSE)
    mtext(pb_names[i], 1, -3, cex = textcex, font = 1,
          adj = 0.57)
  }

  if (legend) {
    # Legend
    par(fig = fig_params[n_row * n_col + 1, ], new = TRUE, xpd = TRUE)
    legend(
        x = -1500000, y = 35000000 / (n_row), cex = textcex,
        legend = c("safe zone", "increasing risk", "high risk"), horiz = F,
        pch = 22, pt.bg = cols[2:4], pt.lwd = 0.0, pt.cex = 1.6,
        box.col = NA, xpd = NA, bg = NA
    )
  }
  if (!is.null(file_name)) dev.off()
}


# convert lpjml vector to raster and change projection to robinson
to_raster <- function(lpjml_array, boundary_box, ext, to_robinson) {

  crs_init <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  lpj_ras <- raster::raster(res = 0.5, crs = crs_init)
  lpj_ras[cellFromXY(lpj_ras, cbind(lpjmliotools::lon, lpjmliotools::lat))] <-
        lpjml_array
  if (to_robinson) {
    ras_to <- raster(xmn = -18000000,
                     xmx = 18000000,
                     ymn = -9000000,
                     ymx = 9000000,
                     crs = "+proj=robin",
                     nrows = 2 * 360,
                     ncols = 2 * 720)

    out_ras <- raster::crop(lpj_ras, ext) %>%
               raster::projectRaster(to = ras_to, crs = "+proj=robin",
                                     na.rm = TRUE, method = "ngb") %>%
               raster::mask(boundary_box) %>%
               suppressWarnings()
  } else {
    out_ras <- raster::crop(lpj_ras, ext)
  }
  return(out_ras)
}

# define plotting locations for par(fig = ...), depending on the number
#   of maps and whether or not a legend is plotted
definefig <- function(n_row, n_col, legfrac) {
  n <- (n_row * n_col) + 1 # number of plots
  rs <- (1 - legfrac) / n_row # rowsize
  top <- 1
  figs <- array(0, dim = c(n, 4)) # dim=c(number of plots, length(x1,x2,y1,y2))
  for (nr in 0:(n_row - 1)) {
    for (nc in 1:n_col) {
      figs[(nr * n_col) + nc, ] <- c(
        (nc - 1) / n_col, (nc) / n_col, top - (rs * (nr + 1)), top - (rs * nr)
      )
    } # of n_col loop
  } # of n_row loop
  figs[n, ] <- c(0, 1, 0, legfrac)
  return(figs)
}
