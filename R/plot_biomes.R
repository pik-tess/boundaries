#' Plot global distribution of lpjml simulated biomes
#'
#' Plots a map with the biome distribution as derived from a lpjml run based
#' on the "classify_biomes" function
#'
#' @param biome_data output (list) from classify_biomes()
#'
#' @param file_name directory for saving the plot (character string)
#'
#' @param to_robinson logical to define if robinson projection should be used
#' for plotting
#'
#' @param bg_col character, specify background possible (`NA` for transparent)#
#'
#' @param grid_path character string providing the path to a grid file
#'
#' @examples
#' \dontrun{
#' plot_biomes(
#'   biome_data = biomes,
#'   file_name = "/p/projects/open/Johanna/R/biomes.pfd"
#'   grid_path = ".grid.bin.json"
#' )
#' }
#'
#' @md
#' @export

plot_biomes <- function(biome_data,
                        file_name = NULL,
                        to_robinson = TRUE,
                        bg_col = "white",
                        grid_path = NULL) {
  # load required data: bbox, countries
  lpjml_extent <- c(-180, 180, -60, 85) # nolint:commented_code_linter

  countries <- system.file(
    "extdata",
    "ne_110m_admin_0_countries.shp",
    package = "boundaries"
  ) %>%
    rgdal::readOGR(layer = "ne_110m_admin_0_countries", verbose = FALSE) %>%
    raster::crop(., lpjml_extent) %>%
    {
      if (to_robinson) sp::spTransform(., sp::CRS("+proj=robin")) else .
    } # nolint

  biome_cols <- c(
    "#993404", "#D95F0E", "#004529", "#238443",
    "#D9F0A3", "#4EB3D3", "#2B8CBE", "#c4e2f4",
    "#FE9929", "#FEC44F", "#FEE391", "#A8DDB5",
    "#E0F3DB", "#F7FCF0", "#c79999", "#0868AC",
    "#FFFFD4", "white", "#dad4d4"
  )

  biome_mapping <- system.file(
    "extdata",
    "biomes.csv",
    package = "boundaries"
  ) %>%
    readr::read_delim(delim = ";", col_types = readr::cols())
  names(biome_cols) <- biome_mapping$short_name

  order_legend <- c(
    1, 2, 9, 10, 11, 3, 4, 5, 12, 13, 14, 6, 7, 8, 15, 16, 17, 18, 19
  )

  biome_cols_legend <- biome_cols[order_legend]

  biome_names_legend <- biome_mapping$short_name[order_legend]

  # TODO: check if this is correct (function arguments changed meanwhile)
  biomes_lpjml <- to_raster(
    lpjml_array = biome_data$biome_id,
    projection = ifelse(to_robinson, "+proj=robin", "+proj=longlat"),
    grid_path = grid_path
  )

  if (!is.null(file_name)) {
    file_extension <- strsplit(file_name, split = "\\.")[[1]][-1]
    switch(file_extension,
      `png` = {
        grDevices::png(file_name,
          width = 8 * 1.8,
          height = 4 * 2,
          units = "cm",
          res = 600,
          pointsize = 7
        )
      },
      `pdf` = {
        grDevices::pdf(file_name,
          width = 8 * 1.8 / 2.54,
          height = (4 * 2) / 2.54,
          pointsize = 7
        )
      },
      {
        stop("File extension ", dQuote(file_extension), " not supported.")
      }
    )
  }
  brk <- seq(
    min(biome_mapping$id) - 0.5,
    max(biome_mapping$id, na.rm = TRUE) + 0.5, 1
  )
  graphics::par(mar = c(4, 0, 0, 0), xpd = TRUE, bg = bg_col) # nolint:undesirable_function_linter
  raster::image(biomes_lpjml,
    asp = 1, xaxt = "n", yaxt = "n",
    xlab = "", ylab = "", col = biome_cols, breaks = brk, lwd = 0.1,
    bty = "n"
  )
  raster::plot(
    countries,
    add = TRUE,
    lwd = 0.3,
    border = "#5c565667",
    usePolypath = FALSE
  )
  if (to_robinson == TRUE) {
    ypoint <- (-6736039)
  } else {
    ypoint <- (-67)
  }
  graphics::legend(0,
    y = ypoint, xjust = 0.45, yjust = 1, cex = 0.8,
    biome_names_legend[1:19],
    fill = biome_cols_legend[1:19],
    horiz = FALSE, border = NULL, bty = "o", box.col = "white",
    bg = bg_col, ncol = 4
  )
  if (!is.null(file_name)) grDevices::dev.off()
}
