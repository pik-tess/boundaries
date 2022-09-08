#' Plot the global status of planetary boundaries
#'
#' Plot global map(s) with the status of 1-4 planetary boundaries for
#' a scenario LPJmL run and derived planetary boundary statuses
#'
#' @param status_data list with one to four LPJmL vectors (length: ncells) with
#' the outputs from the calc_[...] functions for one point in time.
#' Given element names within the list will be displayed underneath the
#' respective map. Order in the list determines the plotting order (by row).
#'
#' @param colors definition of colors for plotting of the safe, increasing
#' risk and high risk zone (3 elements, do not change the order). Colors used
#' for previous plots: "#008B0099" (green), "#FFEB00", "brown3"
#'
#' @param legend boolean to specify whether a legend should be plotted
#'
#' @param plot_to_screen if true, the plot is not saved but displayed on
#' the screen
#'
#' @param pdf If true, the plot will be saved in pdf format; if false in png
#' format
#'
#' @param plotpath directory for saving the plot (character string)
#'
#' @param plotname character string for the name of the plot for saving
#'
#' @examples
#' \dontrun{
#'  plot_pb_status(status_data = list("land system change" = lsc_status,
#'  "nitrogen" = nitrogen_status, "bluewater" = water_status),
#'  plotpath = "/home/plots", plotname = "test", legend = FALSE)
#' }
#'
#' @md
#' @export

plot_pb_status <- function(status_data = list("biosphere" = NA,
                                              "land-system change" = NA,
                                              "water" = NA,
                                              "nitrogen" = NA),
                           colors = c("safe zone" = "#74bca9e1",
                                      "increasing risk" = "#f6ee0f",
                                      "high risk" = "#e23a50"),
                           legend = TRUE,
                           plot_to_screen = FALSE,
                           pdf = FALSE,
                           plotpath = NA,
                           plotname = NA) {
  # checking
  if (length(status_data) > 4) {
    stop(paste0("Number of elements in status data (", length(status_data),
            ") larger than 4. Only up to four PB status maps can be",
            " plotted -- stopping"))
  }
  if (length(colors) != 3) {
    stop(paste0("Length of color vector (", length(colors), ") incorrect. ",
                 "Three colors have to be provided"))
  }

  # load required data: bbox, countries
    extent <- c(-180, 180, -60, 85)

    bounding_box <- system.file("extdata", "ne_110m_wgs84_bounding_box.shp",
                                package = "pbfunctions") %>%
        rgdal::readOGR(layer = "ne_110m_wgs84_bounding_box") %>%
        spTransform(., sp::CRS("+proj=robin"))

    countries <- system.file("extdata", "ne_110m_admin_0_countries.shp",
                                package = "pbfunctions") %>%
        rgdal::readOGR(layer = "ne_110m_admin_0_countries") %>%
        crop(., extent) %>%
        spTransform(., CRS("+proj=robin"))


  pb_names <- names(status_data)

  plot_nat <- toraster_rob(array(0, 67420), bounding_box, extent)

  # plot settings
  if (length(status_data) == 1) {
    nrow <- 1
    ncol <- 1
  } else if (length(status_data) == 2) {
    nrow <- 1
    ncol <- 2
  } else {
    nrow <- 2
    ncol <- 2
  }

if (legend == TRUE) {
  lfrac <- 0.1
  leg_adj <- 1
} else {
  lfrac <- 0
  leg_adj <- 0
}
  fig_params <- definefig(nrow, ncol, lfrac)
  if (plot_to_screen == FALSE) {
    if (pdf == FALSE) {
      png(paste0(plotpath, "/", plotname, ".png"), width = 8 * ncol,
      height = 4 * nrow + leg_adj, units = "cm", res = 600, pointsize = 7)
    } else {
      pdf(paste0(plotpath, "/", plotname, ".pdf"), width = 8 * ncol / 2.54,
      height = (4 * nrow + leg_adj) / 2.54, pointsize = 7)
    }
  }
    textcex <- 0.5 + 0.125 * (nrow + ncol)
    par(mar = rep(0, 4), xpd = TRUE) #, bg = NA)
    brk <- c(-1:4)
    cols <- c("grey92", colors, "darkgrey")

    for (i in seq_len(length(status_data))) {
      if (i == 1) {
        par(fig = fig_params[i, ])
      } else {
        par(fig = fig_params[i, ], new = TRUE)
      }
      # convert lpjml vector to raster with robinson projection
      plotvar <- toraster_rob(unlist(status_data[i]), bounding_box,
                              ext = extent)

      image(plot_nat, zlim = c(-1, 0), asp = 1, xaxt = "n", yaxt = "n",
            xlab = "", ylab = "", col = "grey90", lwd = 0.1, bty = "n")
      image(plotvar, asp = 1, xaxt = "n",
            yaxt = "n", xlab = "", ylab = "", col = cols, breaks = brk,
            lwd = 0.1, bty = "n", add = TRUE)
      plot(countries, add = TRUE, lwd = 0.3, border = "#33333366",
           usePolypath = FALSE)
      mtext(pb_names[i], 1, -3, cex = textcex, font = 1,
            adj = 0.57)
    }

    if (legend == TRUE) {
      # Legend
      par(fig = fig_params[nrow * ncol + 1, ], new = TRUE, xpd = TRUE)
      legend(
          x = -1500000, y = 35000000 / (nrow), cex = textcex,
          legend = c("safe zone", "increasing risk", "high risk"), horiz = F,
          pch = 22, pt.bg = cols[2:4], pt.lwd = 0.0, pt.cex = 1.6,
          box.col = NA, xpd = NA, bg = NA
      )
    }
    if (plot_to_screen == FALSE) dev.off()
}


# convert lpjml vector to raster and change projection to robinson
toraster_rob <- function(lpjml_raster, boundary_box, ext) {
  ras_to <- raster(xmn = -18000000,
                      xmx = 18000000,
                      ymn = -9000000,
                      ymx = 9000000,
                      crs = "+proj=robin",
                      nrows = 2 * 360,
                      ncols = 2 * 720)
  crs_init <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  lpj_ras <- raster::raster(res = 0.5, crs = crs_init)
  lpj_ras[cellFromXY(lpj_ras, cbind(lpjmliotools::lon, lpjmliotools::lat))] <-
        lpjml_raster
  rob_ras <- raster::crop(lpj_ras, ext) %>%
             raster::projectRaster(to = ras_to, crs = "+proj=robin",
                                   na.rm = TRUE, method = "ngb") %>%
             raster::mask(boundary_box) %>%
             suppressWarnings()
  return(rob_ras)
}

# define plotting locations for par(fig = ...), depending on the number
# of maps and whether or not a legend is plotted
definefig <- function(nrow, ncol, legfrac) {
  n <- (nrow * ncol) + 1 # number of plots
  rs <- (1 - legfrac) / nrow # rowsize
  top <- 1
  figs <- array(0, dim = c(n, 4)) # dim=c(number of plots, length(x1,x2,y1,y2))
  for (nr in 0:(nrow - 1)) {
    for (nc in 1:ncol) {
      figs[(nr * ncol) + nc, ] <- c(
        (nc - 1) / ncol, (nc) / ncol, top - (rs * (nr + 1)), top - (rs * nr)
      )
    } # of ncol loop
  } # of nrow loop
  figs[n, ] <- c(0, 1, 0, legfrac)
  return(figs)
}
