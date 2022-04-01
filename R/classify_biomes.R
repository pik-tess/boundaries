#' Classify biomes
#'
#' Classify biomes based on foliage protected cover (FPC) LPJmL output
#'
#' @param path_data output directory (character string) of the LPJmL run where
#' binary files (soon with metafiles) are written
#'
#' @param time_span time span to be used, defined as an integer vector,
#' e.g. `1982:2011` (default)
#'
#' @param vegc_proxy logical. Use vegetation carbon (7500 gC/m2) as a proxy
#' threshold to distinguish forests and savannahs

#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#' time series analysis
#'
#' @return array with cell and second dimension if avg_nyear_args is specified
#'
#' @md
#' @export
classify_biomes <- function(path_data,
                            time_span = c(1982, 2011),
                            vegc_proxy = TRUE,
                            avg_nyear_args = list(),
                            # to be replaced by lpjmlKit::read_output
                            start_year = 1901) {

  # TO BE REPLACED BY lpjmlKit::read_output ---------------------------------- #
  #   hardcoded values to be internally replaced
  # read grid
  ncell <- 67420
  size <- 2
  grid_file <- file(paste(path_data, "grid.bin", sep = "/"), "rb")
  lpjml_grid <- readBin(grid_file, integer(), n = 2 * ncell, size = size) /
                100
  close(grid_file)
  dim(lpjml_grid) <- c(coordinate = 2, cell = ncell)
  dimnames(lpjml_grid) <- list(coordinate = c("lon", "lat"),
                               cell = seq_len(ncell))

  # biome_names after biome classification in Ostberg et al. 2013
  # (https://doi.org/10.5194/esd-4-347-2013), Ostberg et al 2015
  # (https://doi.org/10.1088/1748-9326/10/4/044011) and Gerten et al. 2020
  # (https://doi.org/10.1038/s41893-019-0465-1)
  # biome names
  biome_names <- c("Tropical Rainforest", # 1
                   "Tropical Seasonal & Deciduous Forest", # 2
                   "Temperate Broadleaved Evergreen Forest", # 3
                   "Temperate Broadleaved Deciduous Forest", # 4
                   "Mixed Forest", # 6
                   "Temperate Coniferous Forest", # 5
                   "Boreal Evergreen Forest", # 7
                   "Boreal Deciduous Forest", # 8
                   "Warm Woody Savanna, Woodland & Shrubland", # 9
                   "Warm Savanna & Open Shrubland", #10
                   "Warm Grassland", #11
                   "Temperate Woody Savanna, Woodland & Shrubland", # 12
                   "Temperate Savanna & Open Shrubland", #13
                   "Temperate Grassland", #14
                   "Arctic Tundra", #15
                   "Desert", #16
                   "Rocks and Ice",
                   "Water"
  )

  # band names for fpc.bin
  fpc_names <- c("natvegfrac", # 1
                 "Tropical Broadleaved Evergreen Tree", # 2
                 "Tropical Broadleaved Raingreen Tree", # 3
                 "Temperate Needleleaved Evergreen Tree", # 4
                 "Temperate Broadleaved Evergreen Tree", # 5
                 "Temperate Broadleaved Summergreen Tree", # 6
                 "Boreal Needleleaved Evergreen Tree", # 7
                 "Boreal Broadleaved Summergreen Tree", # 8
                 "Boreal Needleleaved Summergreen Tree", # 9
                 "Tropical C4 grass", # 10
                 "Temperate C3 grass", # 11
                 "Polar C3 grass") # 12

  # indices for pft subsets
  fpc_tropical_trees <- c("Tropical Broadleaved Evergreen Tree",
                          "Tropical Broadleaved Raingreen Tree")
  fpc_temperate_trees <- c("Temperate Needleleaved Evergreen Tree",
                           "Temperate Broadleaved Evergreen Tree",
                           "Temperate Broadleaved Summergreen Tree")
  fpc_boreal_trees <- c("Boreal Needleleaved Evergreen Tree",
                        "Boreal Broadleaved Summergreen Tree",
                        "Boreal Needleleaved Summergreen Tree")
  fpc_needle_trees <- c("Temperate Needleleaved Evergreen Tree",
                        "Boreal Needleleaved Evergreen Tree",
                        "Boreal Needleleaved Summergreen Tree")
  fpc_evergreen_trees <- c("Tropical Broadleaved Evergreen Tree",
                           "Temperate Needleleaved Evergreen Tree",
                           "Temperate Broadleaved Evergreen Tree",
                           "Temperate Broadleaved Summergreen Tree")
  fpc_grass <- c("Tropical C4 grass",
                 "Temperate C3 grass",
                 "Polar C3 grass")
  fpc_trees <- c(fpc_tropical_trees,
                 fpc_temperate_trees,
                 fpc_boreal_trees)

  # TO BE REPLACED BY lpjmlKit::read_output ---------------------------------- #
  #   hardcoded values to be internally replaced
  # read foliage protected cover (fpc) output
  fpc %<-% tmp_read_pft_yearly(
    file_name = paste0(path_data,
                       "/",
                       "fpc.bin"),
    time_span = time_span,
    start_year = start_year,
    ncell = ncell,
    nbands = 12,
    band_names = fpc_names,
    size = 4
  )

  # read vegetation carbon output
  vegc %<-% tmp_read_yearly(
    file_name = paste0(path_data,
                       "/",
                       "vegc.bin"),
    time_span = time_span,
    start_year = start_year,
    ncell = ncell,
    nbands = 1,
    size = 4
  )

  # read temperature output
  temp %<-% tmp_read_yearly(
    file_name = paste0(path_data,
                       "/",
                       "temp.bin"),
    time_span = time_span,
    start_year = start_year,
    ncell = ncell,
    nbands = 1,
    size = 4
  )

  # average fpc
  avg_fpc %<-% do.call(average_nyear_window,
                       append(list(x = fpc),
                              avg_nyear_args))


  # average fpc
  avg_vegc %<-% do.call(average_nyear_window,
                        append(list(x = vegc),
                               avg_nyear_args))


  # average fpc
  avg_temp %<-% do.call(average_nyear_window,
                        append(list(x = temp),
                               avg_nyear_args))

  # latitudes (same dimension for vectorized biome classification)
  latitudes <- array(
    lpjmlKit::subset_array(lpjml_grid, list(coordinate = "lat")),
    dim = dim(avg_temp)
  )

  third_dim <- names(dim(avg_fpc))[
        !names(dim(avg_fpc)) %in% c("cell", "band")
      ] %>% {
        if (rlang::is_empty(.)) NULL else .
      }

  fpc_tree_total %<-% apply(
    lpjmlKit::subset_array(avg_fpc, list(band = fpc_trees)),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  fpc_tree_tropical %<-% apply(
    lpjmlKit::subset_array(avg_fpc, list(band = fpc_tropical_trees)),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  fpc_tree_temperate %<-% apply(
    lpjmlKit::subset_array(avg_fpc, list(band = fpc_temperate_trees)),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  fpc_tree_boreal %<-% apply(
    lpjmlKit::subset_array(avg_fpc, list(band = fpc_boreal_trees)),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  fpc_tree_needle %<-% apply(
    lpjmlKit::subset_array(avg_fpc, list(band = fpc_needle_trees)),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  fpc_tree_evergreen %<-% apply(
    lpjmlKit::subset_array(avg_fpc, list(band = fpc_evergreen_trees)),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  fpc_grass_total %<-% apply(
    lpjmlKit::subset_array(avg_fpc, list(band = fpc_grass)),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  fpc_total %<-% apply(
    lpjmlKit::subset_array(avg_fpc,
                           list(band = fpc_names[
                             which(fpc_names != "natvegfrac")
                           ])),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  max_share_trees %<-% apply(
    lpjmlKit::subset_array(avg_fpc, list(band = fpc_trees)),
    c("cell", third_dim),
    max,
    na.rm = TRUE
  )
  max_share %<-% apply(
    lpjmlKit::subset_array(avg_fpc,
                           list(band = fpc_names[
                             which(fpc_names != "natvegfrac")
                           ])),
    c("cell", third_dim),
    max,
    na.rm = TRUE
  )
  fpc_tree_broadleaf <- fpc_tree_total - fpc_tree_needle

  # initiate biome_class array
  biome_class <- array(NA, dim = dim(avg_vegc), dimnames = dimnames(avg_vegc))

  # use vegc 7500 gC/m2 as proxy threshold for forest/savannah "boundary"
  if (vegc_proxy) {
    is_tropical_proxy <- avg_vegc >= 7500
    is_savannah_proxy <- avg_vegc < 7500
  } else {
    is_tropical_proxy <- TRUE
    is_savannah_proxy <- FALSE
  }

  # Arctic Tundra
  is_arctic_tundra %<-% {
    fpc_total <= 0.05 &
    avg_temp < -2
  }
  # Desert
  is_desert %<-% {
    fpc_total <= 0.05 &
    avg_temp >= -2
  }

  # FORESTS ------------------------------------------------------------------ #
  is_forest %<-% {
    fpc_total > 0.05 &
    fpc_tree_total >= 0.6
  }
  # Boreal Evergreen
  is_boreal_evergreen %<-% {
    is_forest &
    lpjmlKit::subset_array(avg_fpc,
                           list(band = "Boreal Needleleaved Evergreen Tree")) ==
      max_share_trees &
    fpc_tree_broadleaf < (0.4 * fpc_tree_total)
  }
  # Boreal Deciduous
  is_boreal_deciduous %<-% {
    is_forest &
    (lpjmlKit::subset_array(avg_fpc,
                           list(band = "Boreal Broadleaved Summergreen Tree")) == # nolint
      max_share_trees |
    lpjmlKit::subset_array(avg_fpc,
                           list(band = "Boreal Needleleaved Summergreen Tree")) == # nolint
      max_share_trees) &
    fpc_tree_evergreen < (0.4 * fpc_tree_total)
  }
  # Temperate Coniferous Forest
  is_temperate_coniferous %<-% {
    is_forest &
    lpjmlKit::subset_array(avg_fpc,
                           list(band = "Temperate Needleleaved Evergreen Tree")) == # nolint
      max_share_trees &
    fpc_tree_broadleaf < (0.4 * fpc_tree_total)
  }
  # Temperate Broadleaved Evergreen Forest
  is_temperate_broadleaved_evergreen %<-% { # nolint
    is_forest &
    lpjmlKit::subset_array(avg_fpc,
                           list(band = "Temperate Broadleaved Evergreen Tree")) == # nolint
      max_share_trees &
    fpc_tree_tropical < (0.4 * fpc_tree_total) &
    fpc_tree_needle < (0.4 * fpc_tree_total)
  }
  # Temperate Broadleaved Deciduous Forest
  is_temperate_broadleaved_deciduous %<-% { # nolint
    is_forest &
    lpjmlKit::subset_array(avg_fpc,
                           list(band = "Temperate Broadleaved Summergreen Tree")) == # nolint
      max_share_trees &
    fpc_tree_tropical < (0.4 * fpc_tree_total) &
    fpc_tree_needle < (0.4 * fpc_tree_total)
  }
  # Tropical Rainforest
  is_tropical_evergreen %<-% {
    is_forest &
    lpjmlKit::subset_array(avg_fpc,
                           list(band = "Tropical Broadleaved Evergreen Tree")) == # nolint
      max_share_trees &
    (fpc_tree_boreal + fpc_tree_temperate) < (0.4 * fpc_tree_total) &
    is_tropical_proxy
  }
  # Tropical Seasonal & Deciduous Forest
  is_tropical_raingreen %<-% {
    is_forest &
    (lpjmlKit::subset_array(avg_fpc,
                           list(band = "Tropical Broadleaved Raingreen Tree")) == # nolint
      max_share_trees) &
    (fpc_tree_boreal + fpc_tree_temperate) < (0.4 * fpc_tree_total) &
    is_tropical_proxy
  }
  # Warm Woody Savanna, Woodland & Shrubland
  is_tropical_forest_savannah %<-% {
    is_forest &
    (lpjmlKit::subset_array(avg_fpc,
                           list(band = "Tropical Broadleaved Evergreen Tree")) == # nolint
      max_share_trees |
    lpjmlKit::subset_array(avg_fpc,
                           list(band = "Tropical Broadleaved Raingreen Tree")) == # nolint
      max_share_trees) &
    (fpc_tree_boreal + fpc_tree_temperate) < (0.4 * fpc_tree_total) &
    is_savannah_proxy
  }
  is_mixed_forest %<-% {
    is_forest &
    !is_boreal_evergreen &
    !is_boreal_deciduous &
    !is_temperate_coniferous &
    !is_temperate_broadleaved_evergreen &
    !is_temperate_broadleaved_deciduous &
    !is_tropical_evergreen &
    !is_tropical_raingreen &
    !is_tropical_forest_savannah
  }

  # WOODY SAVANNAH ----------------------------------------------------------- #
  is_woody_savannah %<-% {
    fpc_total > 0.05 &
    fpc_tree_total < 0.6 &
    fpc_tree_total >= 0.3
  }
  # Temperate Woody Savanna, Woodland & Shrubland
  is_temperate_woody_savannah %<-% {
    is_woody_savannah &
    lpjmlKit::subset_array(avg_fpc, list(band = "Temperate C3 grass")) >
    lpjmlKit::subset_array(avg_fpc, list(band = "Tropical C4 grass")) &
    avg_temp >= -2 &
    latitudes < 55
  }
  # Warm Woody Savanna, Woodland & Shrubland
  is_tropical_woody_savannah %<-% {
    is_woody_savannah &
    lpjmlKit::subset_array(avg_fpc, list(band = "Temperate C3 grass")) <
    lpjmlKit::subset_array(avg_fpc, list(band = "Tropical C4 grass"))
  }
  # Arctic Tundra
  is_woody_arctic_tundra %<-% {
    is_woody_savannah &
    !is_temperate_woody_savannah &
    !is_tropical_woody_savannah
  }

  # OPEN SHRUBLAND / SAVANNAHS ----------------------------------------------- #

  is_shrubbland %<-% {
    fpc_total > 0.05 &
    fpc_tree_total < 0.3 &
    fpc_tree_total >= 0.1
  }
  # Temperate Savanna & Open Shrubland
  is_temperate_shrubland %<-% {
    is_shrubbland &
    lpjmlKit::subset_array(avg_fpc, list(band = "Temperate C3 grass")) >
    lpjmlKit::subset_array(avg_fpc, list(band = "Tropical C4 grass")) &
    avg_temp >= -2 &
    latitudes < 55
  }
  # Warm Savanna & Open Shrubland
  is_tropical_shrubland %<-% {
    is_shrubbland &
    lpjmlKit::subset_array(avg_fpc, list(band = "Temperate C3 grass")) <
    lpjmlKit::subset_array(avg_fpc, list(band = "Tropical C4 grass")) &
    avg_temp >= -2
  }
  # Arctic Tundra
  is_arctic_shrubland %<-% {
    is_shrubbland &
    !is_temperate_shrubland &
    !is_tropical_shrubland
  }

  # GRASSLAND ---------------------------------------------------------------- #
  is_grassland %<-% {
    fpc_total > 0.05 &
    fpc_tree_total < 0.1
  }
  # Temperate Savanna & Open Shrubland
  is_temperate_grassland %<-% {
    is_shrubbland &
    lpjmlKit::subset_array(avg_fpc, list(band = "Temperate C3 grass")) >
    lpjmlKit::subset_array(avg_fpc, list(band = "Tropical C4 grass")) &
    avg_temp >= -2 &
    latitudes < 55
  }
  # Warm Savanna & Open Shrubland
  is_tropical_grassland %<-% {
    is_grassland &
    lpjmlKit::subset_array(avg_fpc, list(band = "Temperate C3 grass")) <
    lpjmlKit::subset_array(avg_fpc, list(band = "Tropical C4 grass")) &
    avg_temp >= -2
  }
  # Arctic Tundra
  is_arctic_grassland %<-% {
    is_grassland &
    (!is_temperate_grassland &
    !is_tropical_grassland)
  }
  # Rocks and Ice
  is_rocks_and_ice %<-% {
    fpc_total == 0 &
    avg_temp < -2
  }
  # Water body
  is_water %<-% {
    lpjmlKit::subset_array(avg_fpc, list(band = "natvegfrac")) == 0
  }

  # CLASSIFY BIOMES ---------------------------------------------------------- #

  biome_class[
    which(is_arctic_tundra)
  ] <- which(biome_names == "Arctic Tundra")

  biome_class[
    which(is_desert)
  ] <- which(biome_names == "Desert")

  biome_class[
    which(is_boreal_evergreen)
  ] <- which(biome_names == "Boreal Evergreen Forest")

  biome_class[
    which(is_boreal_deciduous)
  ] <- which(biome_names == "Boreal Deciduous Forest")

  biome_class[
    which(is_temperate_coniferous)
  ] <- which(biome_names == "Temperate Coniferous Forest")

  biome_class[
    which(is_temperate_broadleaved_evergreen)
  ] <- which(biome_names == "Temperate Broadleaved Evergreen Forest")

  biome_class[
    which(is_temperate_broadleaved_deciduous)
  ] <- which(biome_names == "Temperate Broadleaved Deciduous Forest")

  biome_class[
    which(is_tropical_evergreen)
  ] <- which(biome_names == "Tropical Rainforest")

  biome_class[
    which(is_tropical_raingreen)
  ] <- which(biome_names == "Tropical Seasonal & Deciduous Forest")

  biome_class[
    which(is_tropical_forest_savannah)
  ] <- which(biome_names == "Warm Woody Savanna, Woodland & Shrubland")

  biome_class[
    which(is_mixed_forest)
  ] <- which(biome_names == "Mixed Forest")

  biome_class[
    which(is_temperate_woody_savannah)
  ] <- which(biome_names == "Temperate Woody Savanna, Woodland & Shrubland")

  biome_class[
    which(is_tropical_woody_savannah)
  ] <- which(biome_names == "Warm Woody Savanna, Woodland & Shrubland")

  biome_class[
    which(is_woody_arctic_tundra)
  ] <- which(biome_names == "Arctic Tundra")

  biome_class[
    which(is_temperate_shrubland)
  ] <- which(biome_names == "Temperate Savanna & Open Shrubland")

  biome_class[
    which(is_tropical_shrubland)
  ] <- which(biome_names == "Warm Savanna & Open Shrubland")

  biome_class[
    which(is_arctic_shrubland)
  ] <- which(biome_names == "Arctic Tundra")

  biome_class[
    which(is_temperate_grassland)
  ] <- which(biome_names == "Temperate Grassland")

  biome_class[
    which(is_tropical_grassland)
  ] <- which(biome_names == "Warm Grassland")

  biome_class[
    which(is_arctic_grassland)
  ] <- which(biome_names == "Arctic Tundra")

  biome_class[
    which(is_rocks_and_ice)
  ] <- which(biome_names == "Rocks and Ice")

  biome_class[
    which(is_water)
  ] <- which(biome_names == "Water")

  return(biome_class)
}


tmp_read_yearly <- function(file_name,
                             time_span,
                             start_year,
                             ncell,
                             nbands,
                             size) {
  # scenario runofv
  file_con <- file(file_name,
                   "rb")
  seek(file_con,
       where = (time_span[1] - start_year) *
                nbands * ncell * size,
       origin = "start")
  lpjml_data <- readBin(file_con,
                        double(),
                         n = (ncell * nbands *
                              (time_span[2] -
                               time_span[1] + 1)),
                        size = size)
  close(file_con)
  dim(lpjml_data) <- c(cell = ncell,
                       year = (time_span[2] -
                                time_span[1] + 1))
  dimnames(lpjml_data) <- list(cell = seq_len(ncell),
                               year = seq(time_span[1],
                                           time_span[2]))
  return(lpjml_data)
}


tmp_read_pft_yearly <- function(file_name,
                             time_span,
                             start_year,
                             ncell,
                             nbands,
                             band_names = NULL,
                             size) {
  # scenario runofv
  file_con <- file(file_name,
                   "rb")
  seek(file_con,
       where = (time_span[1] - start_year) *
                nbands * ncell * size,
       origin = "start")
  lpjml_data <- readBin(file_con,
                        double(),
                         n = (ncell * nbands *
                              (time_span[2] -
                               time_span[1] + 1)),
                        size = size)
  close(file_con)
  dim(lpjml_data) <- c(cell = ncell,
                       band = nbands,
                       year = (time_span[2] -
                                time_span[1] + 1))
  dimnames(lpjml_data) <- list(cell = seq_len(ncell),
                               band = band_names,
                               year = seq(time_span[1],
                                           time_span[2]))
  return(lpjml_data)
}