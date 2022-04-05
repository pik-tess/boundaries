# calculate continent mask as extended lpjml_grid array
#   based on ./extdata/world_continents.shp
#   returns lpjml grid with continent vals besides lon & lat (names as dimnames)
calc_continents_mask <- function(path_output, eurasia = TRUE) {

  # TO BE REPLACED BY lpjmlKit::read_output -------------------------------- #
  # read grid
  ncell <- 67420
  size <- 2
  grid_file <- file(paste(path_output, "grid.bin", sep = "/"), "rb")
  lpjml_grid <- readBin(grid_file, integer(), n = 2 * ncell, size = size) /
                100
  close(grid_file)
  dim(lpjml_grid) <- c(coordinate = 2, cell = ncell)
  dimnames(lpjml_grid) <- list(coordinate = c("lon", "lat"),
                               cell = seq_len(ncell))

  # grid to simple feature collection (sfc) 
  #   https://r-spatial.github.io/sf/articles/sf1.html
  grid_sf <- lpjml_grid %>%
    reshape2::melt() %>%
    tibble::as_tibble() %>%
    tidyr::spread(coordinate, value) %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = sf::st_crs("epsg:4326"))

  # read continents from shapefile as sfc
  continents <- system.file("extdata",
                            "world_continents.shp",
                            package = "pbfunctions") %>%
    sf::st_read(crs = "epsg:4326", quiet = TRUE)

  # intersect both for coastal points that do not intersect calculate nearest
  #   feature; also assign continent names
  continent_grid <- grid_sf %>% dplyr::mutate(
    intersection = as.integer(sf::st_intersects(geometry, continents)),
    continent = dplyr::if_else(is.na(intersection),
                               as.character(
                                 continents$CONTINENT)[
                                   sf::st_nearest_feature(., continents)
                                 ],
                               as.character(
                                 continents$CONTINENT)[intersection]
                               ))
  # workaround to overwrite intersection continent factor with actual factor
  continent_grid$intersection[
    is.na(continent_grid$intersection)
  ] <- match(continent_grid$continent, continents$CONTINENT)[
    is.na(continent_grid$intersection)
  ]

  # if eurasia merge Asia and Europe to Eurasia (intersection = 1)
  if (eurasia) {
    continent_grid$continent[
      which(continent_grid$continent %in% c("Asia", "Europe"))
    ] <- "Eurasia"
    continent_grid$intersection[
      which(continent_grid$intersection == 3)
    ] <- 1
  }

  # bind to lpjml_grid 
  continent_mask <- rbind(lpjml_grid, continent = continent_grid$intersection)

  # add additional continent "coordinate"
  dimnames(continent_mask) <- list(coordinate = c("lon", "lat", "continent"),
                                   cell = continent_grid$continent)

  return(continent_mask)
}