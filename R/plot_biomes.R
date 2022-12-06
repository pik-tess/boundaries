#' Plot global distribution of lpjml simulated biomes
#'
#' Plots a map with the biome distribution as derived from a lpjml run based
#' on the "classify_biomes" function
#'
#' @param biome_data output (list) from classify_biomes()
#'
#' @param plotpath directory for saving the plot (character string)
#'
#' @param to_robinson logical to define if robinson projection should be used
#' for plotting
#'
#' @param bg_col character, specify background possible (`NA` for transparent)#
#'
#' @examples
#' \dontrun{
#'  plot_biomes(biome_data = biomes,
#'              file_name ="/p/projects/open/Johanna/R/biomes.pfd")
#' }
#'
#' @md
#' @export

plot_biomes <- function(biome_data,
                        file_name = NULL,
                        to_robinson = TRUE,
                        bg_col = NA) {

  # load required data: bbox, countries
  lpjml_extent <- c(-180, 180, -60, 85)

  bounding_box <- system.file("extdata", "ne_110m_wgs84_bounding_box.shp",
                              package = "pbfunctions") %>%
      rgdal::readOGR(layer = "ne_110m_wgs84_bounding_box", verbose = FALSE) %>%
      { if(to_robinson) sp::spTransform(., sp::CRS("+proj=robin")) else . } # nolint

  countries <- system.file("extdata", "ne_110m_admin_0_countries.shp",
                              package = "pbfunctions") %>%
      rgdal::readOGR(layer = "ne_110m_admin_0_countries", verbose = FALSE) %>%
      crop(., lpjml_extent) %>%
      { if(to_robinson) sp::spTransform(., CRS("+proj=robin")) else . } # nolint

  biome_cols <-  c("#993404", "#D95F0E", "#004529", "#238443",
                   "#78C679", "#D9F0A3", "#4EB3D3", "#2B8CBE",
                   "#FE9929", "#FEC44F", "#FEE391", "#A8DDB5",
                   "#E0F3DB", "#F7FCF0", "#c79999", "#0868AC",
                   "#FFFFD4", "white", "#dad4d4")

  biome_mapping <- system.file("extdata", "biomes.csv",
                              package = "pbfunctions") %>%
                   readr::read_csv2()
  names(biome_cols) <- biome_mapping$short_name

  order_legend <- c(1, 2, 9, 10, 11, 3, 4, 5, 6, 12, 13, 14, 7, 8, 15, 16, 17,
                      18, 19)

  biome_cols_legend <- biome_cols[order_legend]

  biome_names_legend <- biome_mapping$short_name[order_legend]

  biomes_lpjml <- to_raster(lpjml_array = biome_data$biome_id,
                         boundary_box = bounding_box,
                         ext = lpjml_extent,
                         to_robinson = to_robinson)

  if (!is.null(file_name)) {
    file_extension <- strsplit(file_name, split = "\\.")[[1]][-1]
    switch(file_extension,
      `png` = {
        png(file_name,
            width = 8 * 1.8,
            height = 4 * 2,
            units = "cm",
            res = 600,
            pointsize = 7)
      },
      `pdf` = {
        pdf(file_name,
            width = 8 * 1.8 / 2.54,
            height = (4 * 2) / 2.54,
            pointsize = 7)
      }, {
        stop("File extension ", dQuote(file_extension), " not supported.")
      }
    )
  }
  brk <- seq(0.5, 19.5, 1)
  par(mar = c(4, 0, 0, 0), xpd = T, bg = bg_col)
  image(biomes_lpjml, asp = 1, xaxt = "n", yaxt = "n",
          xlab = "", ylab = "", col = biome_cols, breaks = brk, lwd = 0.1,
          bty = "n")
  plot(countries, add = TRUE, lwd = 0.3,
         border = "#5c565667", usePolypath = FALSE)
  legend(0, y = -6736039, xjust = 0.45, yjust = 1, cex = 0.8,
    biome_names_legend[1:19],
    fill = biome_cols_legend[1:19],
    horiz = F, border = NULL, bty = "o", box.col = "white",
    bg = bg_col, ncol = 4)
  if (!is.null(file_name)) dev.off()
}