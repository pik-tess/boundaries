#' Plot global distribution of lpjml simulated biomes
#'
#' Plots a map with the biome distribution as derived from a lpjml run based
#' on the "classify_biomes" function
#'
#' @param x output (list) from classify_biomes()
#'
#' @param filename directory for saving the plot (character string)
#'
#' @param projection character string defining the projection, default set to
#' "+proj=robin"
#'
#' @param grid_path character string providing the path to a grid file
#'
#' @examples
#' \dontrun{
#'
#' biomes <- classify_biomes(
#'  config_reference = path_reference,
#'  time_span_reference = as.character(2008:2017),
#'  savanna_proxy = list(vegc = 7500)
#' )
#'
#' plot_biomes(
#'   x = biomes,
#'   filename = "/p/projects/open/Johanna/R/biomes.pfd"
#'   grid_path = ".grid.bin.json"
#' )
#' }
#'
#' @md
#' @export
plot_biomes <- function(x,
                        filename = NULL,
                        projection = "+proj=robin",
                        grid_path = NULL) {

  biome_mapping <- system.file(
    "extdata",
    "biomes.csv",
    package = "boundaries"
  ) %>%
    readr::read_delim(delim = ";", col_types = readr::cols())

  biome_cols <- c(
    "#993404", "#D95F0E", "#004529", "#238443",
    "#D9F0A3", "#4EB3D3", "#2B8CBE", "#c4e2f4",
    "#FE9929", "#FEC44F", "#FEE391", "#A8DDB5",
    "#E0F3DB", "#F7FCF0", "#c79999", "#0868AC",
    "#FFFFD4", "white", "#dad4d4"
  )
  names(biome_cols) <- biome_mapping$short_name

  # define order of biomes for the legend
  order_legend <- c(
    1, 2, 9, 10, 11, 3, 4, 5, 12, 13, 14, 6, 7, 8, 15, 16, 17, 18, 19
  )
  biome_names_legend <- biome_mapping$short_name[order_legend]

  # create SpatRaster with biome ids
  biomes_lpjml <- to_raster(
    lpjml_array = x$biome_id,
    projection = projection,
    grid_path = grid_path
  )
  levels(biomes_lpjml) <- data.frame(id = order_legend,
                                     biome = biome_names_legend)
  terra::is.factor(biomes_lpjml)

  # define color map
  coltb <- data.frame(value = c(1:19), col = biome_cols)
  terra::coltab(biomes_lpjml) <- coltb

  # get country outlines
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

  p <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = biomes_lpjml) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "#ffffff"),
      panel.grid.major = ggplot2::element_line(linewidth = 0.1,
                                               color = "#8d8b8b"),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5, size = 8),
      plot.margin = ggplot2::margin(0, 0, 0.5, 0, "cm"),
      panel.spacing = ggplot2::unit(0, "lines"),
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::geom_sf(data = world, fill = NA, linewidth = 0.12,
                     color = "#7e7d7d") +
    ggplot2::xlim(terra::ext(biomes_lpjml)[1], terra::ext(biomes_lpjml)[2]) +
    ggplot2::ylim(terra::ext(biomes_lpjml)[3], terra::ext(biomes_lpjml)[4]) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 5, byrow = FALSE))

  if (!is.null(filename)) {
    ggplot2::ggsave(
      filename,
      p,
      width = 21,
      height = 12,
      dpi = 600,
      units = "cm",
      pointsize = 7
    )
  } else {
    # plot maps to screen
    print(p)
    return(p)
  }

}
