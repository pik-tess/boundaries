#' Classify biomes
#'
#' Classify biomes based on foliage protected cover (FPC) and temperature
#' LPJmL output plus either vegetation carbon or pft_lai depending on
#' the savanna_proxy option and elevation if montane_arctic_proxy requires this
#'
#' @param path_data output directory (character string) of the LPJmL run o read
#'        default outputs from (fpc, grid, vegc, pft_lai, temp)
#' @param timespan time span to be used, defined as an integer vector, e.g.
#'       `c(1982,2011)`, to average over (default over all years, else see
#'       \link[pbfunctions]{average_nyear_window})
#' @param input_files optional `list` containing additional input (!) files,
#'        not in the `path_data`, e.g. if temp was not written out:
#'        `list(grid=..., temp = ..., elevation = ...)`
#' @param diff_output_files optional list for specification of output file names
#'        differing from default, which is list(grid = "grid.bin", fpc = "fpc.bin",# nolint
#'        vegc = "vegc.bin", pft_lai = "pft_lai.bin", temp = "temp.bin")
#' @param file_ending replace default file ending. default: ".bin"
#' @param savanna_proxy `list` with either "pft_lai" or "vegc" as
#'        name/key and value in m2/m2 for pft_lai (default = 6) and gC/m2 for
#'        vegc (default would be 7500), Set to `NULL` if no proxy is used.
#' @param montane_arctic_proxy `list` with either "elevation" or "latitude" as
#'        name/key and value in m for elevation (default 1000) and degree for
#'        latitude (default would be 55), Set to `NULL` if no proxy is used.
#' @param tree_cover_thresholds list with minimum tree cover thresholds for
#'        definition of forest, woodland, savanna and grassland. Only changes to
#'        the default have to be included in the list, for the rest the default
#'        is used. Default values:
#'        "boreal forest" = 0.6
#'        "temperate forest" = 0.6
#'        "temperate woodland" = 0.3
#'        "temperate savanna" = 0.1
#'        "tropical forest" = 0.6
#'        "tropical woodland" = 0.3
#'        "tropical savanna" = 0.1
#' @param avg_nyear_args list of arguments to be passed to
#'        \link[pbfunctions]{average_nyear_window} (see for more info). To be used for # nolint
#'        time series analysis
#' @param read_args list of arguments for reading input/output. only required
#'        for:
#'        nc output: specification of header_size and ncell to read in
#'                   lpjml grid input
#'        raw/clm output: specification of header_size, ncell, firstyear and 
#'                   fpc_nbands (12 or 10)
#' @return list object containing biome_id (main biome per grid cell [dim=c(ncells)]), # nolint
#' and list of respective biome_names[dim=c(nbiomes)]
#'
#' @examples
#' \dontrun{
#' classify_biomes(
#'   path_data = "/p/projects/open/Fabian/runs/Gamma/output/historic_gamma"
#'   timespan = c(1982,2011))
#' }
#'
#' @export
classify_biomes <- function(path_data,
                            timespan,
                            input_files = NULL,
                            diff_output_files = NULL,
                            file_type = "raw",
                            savanna_proxy = list(pft_lai = 6),
                            montane_arctic_proxy = list(elevation = 1000),
                            tree_cover_thresholds = list(),
                            avg_nyear_args = list(), # currently a place holder
                            read_args = list(
                              header_size = 0,
                              ncell = 67420, firstyear = 1901, fpc_nbands = 12,
                              size = 4
                            )) {


  file_extension <- switch(file_type,
                           raw = ".bin",
                           meta = ".bin.json",
                           clm = ".clm",
                           nc = ".nc",
                           nc4 = ".nc4",
                           cdf = ".nc")
  # default output files with defined file_extension
  output_files <- list(grid = "grid",
                       fpc = "fpc",
                       vegc = "vegc",
                       pft_lai = "pft_lai",
                       temp = "temp") %>%
    lapply(paste0, file_extension)


  # replace output files that differ from standard output  id
  replace_idx <- match(names(diff_output_files), names(output_files))
  if (any(is.na(replace_idx))) {
    stop(paste0(
      names(diff_output_files)[which(is.na(replace_idx))],
      " is not valid. Please use a name of: ",
      paste0(names(output_files), collapse = ", ")
    ))
  }
  output_files[c(replace_idx)] <- diff_output_files

  #  concatenate path_data and output_files
  output_files <- lapply(output_files,
                         function(x, path_data) {
                           paste(path_data, x, sep = "/")
                         },
                         path_data = path_data)

  # define default minimum tree cover for forest / woodland / savanna
  min_tree_cover <- list("boreal forest" = 0.6, "temperate forest" = 0.6,
                         "temperate woodland" = 0.3, "temperate savanna" = 0.1,
                         "tropical forest" = 0.6, "tropical woodland" = 0.3,
                         "tropical savanna" = 0.1)

  # replace default values by values defined in tree_cover_thresholds
  # parameter -> won't be applied if not specified
  replace_idx <- match(names(tree_cover_thresholds), names(min_tree_cover))
  if (any(is.na(replace_idx))) {
    stop(paste0(
      names(tree_cover_thresholds)[which(is.na(replace_idx))],
      " is not valid. Please use a name of: ",
      paste0(names(min_tree_cover), collapse = ", ")
    ))
  }
  min_tree_cover[replace_idx] <- tree_cover_thresholds

  # test if forest threshold is always > woodland threshold > savanna threshold
  if (min_tree_cover[["temperate forest"]] <=
        min_tree_cover[["temperate woodland"]] |
      min_tree_cover[["temperate woodland"]] <=
        min_tree_cover[["temperate savanna"]] |
      min_tree_cover[["tropical woodland"]] <=
        min_tree_cover[["tropical savanna"]] |
      min_tree_cover[["tropical forest"]] <=
        min_tree_cover[["tropical woodland"]]) {
    stop(paste0("Tree cover threshold for forest are not always higher than",
                "tree cover thresholds for woodland and savanna. Aborting."))
  }

  # test if provided proxies are valid
  savanna_proxy_name <- match.arg(names(savanna_proxy), c(NA, "vegc", "pft_lai"))
  montane_arctic_proxy_name <- match.arg(names(montane_arctic_proxy),
                                         c(NA, "elevation", "latitude"))

  if (file.exists(output_files$grid) && file_type == "meta") {
      grid <- lpjmliotools::autoReadMetaOutput(
        metaFile = output_files$grid
      )
      ncell <- length(grid) / 2
      lon   <- grid[c(1:ncell) * 2 - 1]
      lat   <- grid[c(1:ncell) * 2]
  } else if (file_type %in% c("nc", "cdf", "nc4")) {
    # if nc output is defined, we need an lpjml grid to convert to 
    # the correct array size, this needs to be given in input_files$grid
    message("Reading of netcdf output is still preliminary. Please specify LPJmL grid input.") # nolint
    grid <- lpjmliotools::readGridInputBin(inFile = input_files$grid,
                             headersize = read_args$header_size,
                             ncells = read_args$ncell)
    lon <- grid$lon
    lat <- grid$lat
    ncell <- length(grid$lon)
  } else if (file.exists(output_files$grid) && file_type %in% c("raw", "clm")) {
    grid <- lpjmliotools::readGridOutputBin(inFile = output_files$grid,
                                            headersize = read_args$header_size, #nolint
                                            ncells = read_args$ncell)
    lon <- grid$lon
    lat <- grid$lat
    ncell <- length(grid$lon)
  } else {
    stop(paste0("Output file ",
                output_files$grid,
                " does not exist. Make sure the specified input path_data is ",
                "correct. If your file names differ from the default, please ",
                "use diff_output_files to specify them. "))
  }
  lpjml_grid <- rbind(lon, lat)
  # TODO: convert to yearly if output is monthly
  if (file_type == "meta") {
    fpc <- lpjmliotools::autoReadMetaOutput(
      metaFile = output_files$fpc,
      getyearstart = timespan[1],
      getyearstop = timespan[2]
    )
    if (file.exists(output_files$temp) && is.null(input_files$temp)) {
      temp <- lpjmliotools::autoReadMetaOutput(
        metaFile = output_files$temp,
        getyearstart = timespan[1],
        getyearstop = timespan[2]
      )
    }
    if (!is.na(savanna_proxy_name)) {
      savanna_proxy_data <- lpjmliotools::autoReadMetaOutput(
        metaFile = output_files[[savanna_proxy_name]],
        getyearstart = timespan[1],
        getyearstop = timespan[2]
      )
    }
  } else if (file_type %in% c("nc", "nc4", "cdf")) {
    fpc <- lpjmliotools::netcdfCFT2lpjarray(
      ncInFile = output_files$fpc,
      var = "FPC",
      lon = lon,
      lat = lat
    )
    if (file.exists(output_files$temp) && is.null(input_files$temp)) {
      temp <- lpjmliotools::netcdfCFT2lpjarray(
        ncInFile = output_files$temp,
        var = "temp",
        lon = lon,
        lat = lat
        )
    }
    if (savanna_proxy_name == "vegc") {
      savanna_proxy_data <- lpjmliotools::netcdfCFT2lpjarray(
        ncInFile = files$vegc,
        var = "VegC",
        lon = lon,
        lat = lat
      )
    } else if (savanna_proxy_name == "pft_lai") {
      fpc <- lpjmliotools::netcdfCFT2lpjarray(
        ncInFile = output_files$pft_lai,
        var = "LAI",
        lon = lon,
        lat = lat
      )
    }
  } else if (file_type %in% c("raw", "clm")) {
    fpc <- lpjmliotools::readCFToutput(inFile = output_files$fpc,
                                      startyear = read_args$firstyear,
                                      stopyear = timespan[2],
                                      size = read_args$size,
                                      headersize = read_args$header_size,
                                      getyearstart = timespan[1],
                                      getyearstop = timespan[2],
                                      ncells = read_args$ncell,
                                      bands = read_args$fpc_nbands)
    if (file.exists(output_files$temp) && is.null(input_files$temp)) {
      temp <- lpjmliotools::readDaily(inFile = output_files$temp,
                                      startyear = read_args$firstyear,
                                      stopyear = timespan[2],
                                      size = read_args$size,
                                      headersize = read_args$header_size,
                                      getyearstart = timespan[1],
                                      getyearstop = timespan[2],
                                      ncells = read_args$ncell)
    }
    if (savanna_proxy_name == "vegc") {
      savanna_proxy_data <- lpjmliotools::readYearly(inFile = output_files$vegc,
                                      startyear = read_args$firstyear,
                                      stopyear = timespan[2],
                                      size = read_args$size,
                                      headersize = read_args$header_size,
                                      getyearstart = timespan[1],
                                      getyearstop = timespan[2],
                                      ncells = read_args$ncell)
    } else if (savanna_proxy_name == "pft_lai") {
      savanna_proxy_data <- lpjmliotools::readCFToutput(inFile = output_files$pft_lai, #nolint
                                      startyear = read_args$firstyear,
                                      stopyear = timespan[2],
                                      size = read_args$size,
                                      headersize = read_args$header_size,
                                      getyearstart = timespan[1],
                                      getyearstop = timespan[2],
                                      ncells = read_args$ncell,
                                      bands = c(read_args$fpc_nbands -
                                                     1 + 16 * 2))
    }
  } else {
      stop(paste0("Unknown file ending (",
                  fpc_ending,
                  "). Aborting."))
  }

  if (!is.null(input_files$temp)) {
    temp <- lpjmliotools::autoReadInput(inFile = input_files$temp,
                        getyearstart = timespan[1], getyearstop = timespan[2])
    # monthly temperature
    # TODO has to be processed to yearly?
  }

  if(!file_type %in% c("raw", "clm")) {
    fpc_nbands <- dim(fpc)[["nbands"]]
  }
  npft <- fpc_nbands - 1

  # average fpc
  avg_fpc %<-% do.call(
    average_nyear_window,
    append(list(x = fix_dimnames(fpc, "fpc", ncell, npft, fpc_nbands)),
           avg_nyear_args)
  )

  # average vegc or pft_lai
  if (!is.na(savanna_proxy_name)) {
    avg_savanna_proxy_data %<-% do.call(
      average_nyear_window,
      append(list(x = fix_dimnames(savanna_proxy_data, savanna_proxy_name, ncell, npft, fpc_nbands)), # nolint
             avg_nyear_args)
    )
  }

  # average temp
  avg_temp %<-% do.call(
    average_nyear_window,
    append(list(x = fix_dimnames(temp, "temp", ncell, npft)),
           avg_nyear_args)
  )


  if (montane_arctic_proxy_name == "elevation") {
      elevation <- lpjmliotools::autoReadInput(
        inFile = input_files$elevation
      )[1, ]
  }

  # biome_names after biome classification in Ostberg et al. 2013
  # (https://doi.org/10.5194/esd-4-347-2013), Ostberg et al 2015
  # (https://doi.org/10.1088/1748-9326/10/4/044011) and Gerten et al. 2020
  # (https://doi.org/10.1038/s41893-019-0465-1)

  # biome names
  biome_mapping <- system.file("extdata",
                               "biomes.csv",
                               package = "pbfunctions") %>%
                   readr::read_delim(col_types = readr::cols(), delim = ";")
  biome_names <- biome_mapping$id
  names(biome_names) <- biome_mapping$name


  pft_categories <- system.file("extdata",
                                "pft_categories.csv",
                                package = "pbfunctions") %>%
    read_pft_categories() %>%
    dplyr::filter(npft_proxy == npft)

  fpc_names <- dplyr::filter(pft_categories, category == "natural")$pft

  # indices (when estimation only via npft possible) or names for pft subsets
  fpc_temperate_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & zone == "temperate" & category == "natural"
  ) %>% {
      if (any(is.na(.$npft_proxy))) .$pft else .$lpjml_index + 1
    }

  fpc_tropical_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & zone == "temperate" & category == "natural"
  ) %>% {
    if (any(is.na(.$npft_proxy))) .$pft else .$lpjml_index + 1
  }

  fpc_boreal_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & zone == "boreal" & category == "natural"
  ) %>% {
    if (any(is.na(.$npft_proxy))) .$pft else .$lpjml_index + 1
  }

  fpc_needle_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & category == "needle"
  ) %>% {
    if (any(is.na(.$npft_proxy))) .$pft else .$lpjml_index + 1
  }

  fpc_grass <- dplyr::filter(
    pft_categories,
    type == "grass" & category == "natural"
  ) %>% {
    if (any(is.na(.$npft_proxy))) .$pft else .$lpjml_index + 1
  }

  fpc_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & category == "natural"
  ) %>% {
    if (any(is.na(.$npft_proxy))) .$pft else .$lpjml_index + 1
  }


  # latitudes (same dimension for vectorized biome classification)
  latitudes <- array(
    lpjmliotools::subset_array(lpjml_grid, list(coordinate = "lat")),
    dim = c(ncell)
  )

  third_dim <- names(dim(avg_fpc))[
    !names(dim(avg_fpc)) %in% c("cell", "band")
  ] %>% {
    if (rlang::is_empty(.)) NULL else .
  }

  fpc_tree_total <- apply(
    lpjmliotools::subset_array(avg_fpc, list(band = fpc_trees)),
    c("cell"),
    sum,
    na.rm = TRUE
  )

  fpc_tree_tropical <- apply(
    lpjmliotools::subset_array(avg_fpc, list(band = fpc_tropical_trees)),
    c("cell"),
    sum,
    na.rm = TRUE
  )
  fpc_tree_temperate <- apply(
    lpjmliotools::subset_array(avg_fpc, list(band = fpc_temperate_trees)),
    c("cell"),
    sum,
    na.rm = TRUE
  )
  fpc_tree_boreal <- apply(
    lpjmliotools::subset_array(avg_fpc, list(band = fpc_boreal_trees)),
    c("cell"),
    sum,
    na.rm = TRUE
  )
  fpc_tree_needle <- apply(
    lpjmliotools::subset_array(avg_fpc, list(band = fpc_needle_trees)),
    c("cell"),
    sum,
    na.rm = TRUE
  )
  fpc_tree_evergreen <- apply(
    lpjmliotools::subset_array(avg_fpc, list(band = fpc_evergreen_trees)),
    c("cell"),
    sum,
    na.rm = TRUE
  )
  fpc_grass_total <- apply(
    lpjmliotools::subset_array(avg_fpc, list(band = fpc_grass)),
    c("cell"),
    sum,
    na.rm = TRUE
  )
  fpc_total <- apply(
    lpjmliotools::subset_array(avg_fpc,
                               list(band = fpc_names[fpc_names != "natvegfrac"])), # nolint
    c("cell"),
    sum,
    na.rm = TRUE
  )
  max_share_trees <- apply(
    lpjmliotools::subset_array(avg_fpc, list(band = fpc_trees)),
    c("cell"),
    max,
    na.rm = TRUE
  )

  fpc_tree_broadleaf <- fpc_tree_total - fpc_tree_needle

  # use vegc 7500 gC/m2 or natLAI 6 as proxy threshold for forest/savanna
  #   "boundary
  if (!is.null(savanna_proxy)) {
    if (savanna_proxy_name == "pft_lai") {
      # TODO: dynamic subsetting
      avg_savanna_proxy_data <- rowSums(
        avg_savanna_proxy_data[, 1:npft] * avg_fpc[, 2:(npft + 1)] * avg_fpc[, 1] # nolint
      )
    }
    is_tropical_proxy <- avg_savanna_proxy_data >= savanna_proxy[[savanna_proxy_name]] # nolint
    is_savannah_proxy <- avg_savanna_proxy_data < savanna_proxy[[savanna_proxy_name]] # nolint
  } else {
    is_tropical_proxy <- rep(TRUE, ncell)
    is_savannah_proxy <- rep(FALSE, ncell)
  }

  # Desert
  is_desert <- {
    fpc_total <= 0.05 &
      avg_temp >= 0 #-2
  }

  # montane (for classification of montane grassland)
  if (montane_arctic_proxy_name == "elevation") {
    is_montane_artic <- elevation > montane_arctic_proxy[[
      montane_arctic_proxy_name
    ]]
  } else if (montane_arctic_proxy_name == "elevation") {
    is_montane_artic <- !(abs(latitudes) > montane_arctic_proxy[[
      montane_arctic_proxy_name
    ]])
  }

  # FORESTS ------------------------------------------------------------------ #
  is_boreal_forest <- {
    fpc_tree_total >= min_tree_cover[["boreal forest"]]
  }
  is_temperate_forest <- {
    fpc_tree_total >= min_tree_cover[["temperate forest"]]
  }
  is_tropical_forest <- {
    fpc_tree_total >= min_tree_cover[["tropical forest"]]
  }
  # Boreal Evergreen
  is_boreal_evergreen <- {
    is_boreal_forest &
    lpjmliotools::subset_array(
      avg_fpc,
      list(band = "Boreal Needleleaved Evergreen Tree")
    ) == max_share_trees &
    fpc_tree_broadleaf < (0.4 * fpc_tree_total)
  }

  if (npft == 9) {
    # Boreal Deciduous
    is_boreal_deciduous <- {
      is_boreal_forest &
      (
        lpjmliotools::subset_array(
          avg_fpc,
          list(band = "Boreal Broadleaved Summergreen Tree")
        ) == max_share_trees
      ) &
      fpc_tree_evergreen < (0.4 * fpc_tree_total)
    }
  }else if (npft == 11) {
    # Boreal Deciduous
    is_boreal_deciduous <- {
      is_boreal_forest &
      (
        lpjmliotools::subset_array(
          avg_fpc,
          list(band = "Boreal Broadleaved Summergreen Tree")
        ) == max_share_trees |
        lpjmliotools::subset_array(
          avg_fpc,
          list(band = "Boreal Needleleaved Summergreen Tree")
        ) == max_share_trees
      ) &
      fpc_tree_evergreen < (0.4 * fpc_tree_total)
    }
  } else {
    stop(paste("Unknown number of pfts:", npft))
  }

  # Temperate Coniferous Forest
  is_temperate_coniferous <- {
    is_temperate_forest &
    lpjmliotools::subset_array(
      avg_fpc,
      list(band = "Temperate Needleleaved Evergreen Tree")
    ) == max_share_trees &
    fpc_tree_broadleaf < (0.4 * fpc_tree_total)
  }
  # Temperate Broadleaved Evergreen Forest
  is_temperate_broadleaved_evergreen <- { # nolint
    is_temperate_forest &
    lpjmliotools::subset_array(
      avg_fpc,
      list(band = "Temperate Broadleaved Evergreen Tree")
    ) == max_share_trees &
    fpc_tree_tropical < (0.4 * fpc_tree_total) &
    fpc_tree_needle < (0.4 * fpc_tree_total)
  }
  # Temperate Broadleaved Deciduous Forest
  is_temperate_broadleaved_deciduous <- { # nolint
    is_temperate_forest &
    lpjmliotools::subset_array(
      avg_fpc,
      list(band = "Temperate Broadleaved Summergreen Tree")
    ) == max_share_trees &
    fpc_tree_tropical < (0.4 * fpc_tree_total) &
    fpc_tree_needle < (0.4 * fpc_tree_total)
  }
  # Tropical Rainforest
  is_tropical_evergreen <- {
    is_tropical_forest &
    lpjmliotools::subset_array(
      avg_fpc,
      list(band = "Tropical Broadleaved Evergreen Tree")
    ) == max_share_trees &
    (fpc_tree_boreal + fpc_tree_temperate) < (0.4 * fpc_tree_total) &
    is_tropical_proxy
  }
  # Tropical Seasonal & Deciduous Forest
  is_tropical_raingreen <- {
    is_tropical_forest &
    (lpjmliotools::subset_array(
      avg_fpc,
      list(band = "Tropical Broadleaved Raingreen Tree")
    ) == max_share_trees) &
    (fpc_tree_boreal + fpc_tree_temperate) < (0.4 * fpc_tree_total) &
    is_tropical_proxy
  }
  # Warm Woody Savanna, Woodland & Shrubland
  is_tropical_forest_savanna <- {
    is_tropical_forest &
    (
      lpjmliotools::subset_array(
        avg_fpc,
        list(band = "Tropical Broadleaved Evergreen Tree")
      ) == max_share_trees |
      lpjmliotools::subset_array(
        avg_fpc,
        list(band = "Tropical Broadleaved Raingreen Tree")
      ) == max_share_trees
    ) &
    (fpc_tree_boreal + fpc_tree_temperate) < (0.4 * fpc_tree_total) &
    is_savanna_proxy
  }
  is_mixed_forest <- {
    is_temperate_forest &
      !is_boreal_evergreen &
      !is_boreal_deciduous &
      !is_temperate_coniferous &
      !is_temperate_broadleaved_evergreen &
      !is_temperate_broadleaved_deciduous &
      !is_tropical_evergreen &
      !is_tropical_raingreen &
      !is_tropical_forest_savanna
  }

  # WOODY savanna ----------------------------------------------------------- #

  # Temperate Woody Savanna, Woodland & Shrubland
  is_temperate_woody_savanna <- {
    fpc_tree_total <= min_tree_cover[["temperate forest"]] &
    fpc_tree_total >= min_tree_cover[["temperate woodland"]] &
    lpjmliotools::subset_array(avg_fpc, list(band = "Temperate C3 Grass")) >
    lpjmliotools::subset_array(avg_fpc, list(band = "Tropical C4 Grass")) &
    avg_temp >= 0 #-2 &
    #latitudes < 55
  }
  # Warm Woody Savanna, Woodland & Shrubland
  is_tropical_woody_savanna <- {
    fpc_tree_total <= min_tree_cover[["tropical forest"]] &
    fpc_tree_total >= min_tree_cover[["tropical woodland"]] &
    lpjmliotools::subset_array(avg_fpc, list(band = "Temperate C3 Grass")) <
    lpjmliotools::subset_array(avg_fpc, list(band = "Tropical C4 Grass"))
  }

  # OPEN SHRUBLAND / SAVANNAS ----------------------------------------------- #

  # Temperate Savanna & Open Shrubland
  is_temperate_shrubland <- {
    fpc_tree_total <= min_tree_cover[["temperate woodland"]] &
    fpc_tree_total >= min_tree_cover[["temperate savanna"]] &
    lpjmliotools::subset_array(avg_fpc, list(band = "Temperate C3 Grass")) >
    lpjmliotools::subset_array(avg_fpc, list(band = "Tropical C4 Grass")) &
    avg_temp >= 0 #-2 &
    #latitudes < 55
  }
  # Warm Savanna & Open Shrubland
  is_tropical_shrubland <- {
    fpc_tree_total <= min_tree_cover[["tropical woodland"]] &
    fpc_tree_total >= min_tree_cover[["tropical savanna"]] &
    lpjmliotools::subset_array(avg_fpc, list(band = "Temperate C3 Grass")) <
    lpjmliotools::subset_array(avg_fpc, list(band = "Tropical C4 Grass")) &
    avg_temp >= 0 #-2
  }

  # GRASSLAND ---------------------------------------------------------------- #

  # Temperate Savanna & Open Shrubland
  is_temperate_grassland <- {
    fpc_total > 0.05 &
    fpc_tree_total <= min_tree_cover[["temperate savanna"]] &
    lpjmliotools::subset_array(avg_fpc, list(band = "Temperate C3 Grass")) >
    lpjmliotools::subset_array(avg_fpc, list(band = "Tropical C4 Grass")) &
    avg_temp >= 0 #-2 &
    #latitudes < 55
  }
  # Warm Savanna & Open Shrubland
  is_tropical_grassland <- {
    fpc_total > 0.05 &
    fpc_tree_total <= min_tree_cover[["tropical savanna"]] &
    lpjmliotools::subset_array(avg_fpc, list(band = "Temperate C3 Grass")) <
    lpjmliotools::subset_array(avg_fpc, list(band = "Tropical C4 Grass")) &
    avg_temp >= 0 #-2
  }

  # Arctic Tundra ------------------------------------------------------------ #
  is_arctic_tundra <- {
    (!is_boreal_evergreen &
    !is_boreal_deciduous &
    !is_temperate_forest &
    (
      avg_temp < 0 |
      lpjmliotools::subset_array(avg_fpc, list(band = "Temperate C3 Grass")) ==
      lpjmliotools::subset_array(avg_fpc, list(band = "Tropical C4 Grass"))) &
      fpc_total > 0.05
    ) |
    (avg_temp < 0 & fpc_total < 0.05)
  }

  # Rocks and Ice
  is_rocks_and_ice <- {
    fpc_total == 0 &
      avg_temp < 0 #-2
  }
  # Water body
  is_water <- {
    lpjmliotools::subset_array(avg_fpc, list(band = "natvegfrac")) == 0
  }

  # CLASSIFY BIOMES ---------------------------------------------------------- #

  # initiate biome_class array
  biome_class <- array(NA, dim = c(ncell), dimnames = dimnames(fpc_total))

  biome_class[is_desert] <- biome_names["Desert"]

  # forests
  biome_class[is_boreal_evergreen] <- biome_names["Boreal Evergreen Forest"]
  biome_class[is_boreal_deciduous] <- biome_names["Boreal Deciduous Forest"]
  biome_class[is_temperate_coniferous] <- biome_names["Temperate Coniferous Forest"] # nolint
  biome_class[is_temperate_broadleaved_evergreen] <- biome_names["Temperate Broadleaved Evergreen Forest"] # nolint
  biome_class[is_temperate_broadleaved_deciduous] <- biome_names["Temperate Broadleaved Deciduous Forest"] # nolint
  biome_class[is_tropical_evergreen] <- biome_names["Tropical Rainforest"]
  biome_class[is_tropical_raingreen] <- biome_names["Tropical Seasonal & Deciduous Forest"] # nolint
  biome_class[is_tropical_forest_savanna] <- biome_names["Warm Woody Savanna, Woodland & Shrubland"] # nolint
  biome_class[is_mixed_forest] <- biome_names["Mixed Forest"]

  # woody savanna
  biome_class[is_temperate_woody_savanna] <- biome_names["Temperate Woody Savanna, Woodland & Shrubland"] # nolint
  biome_class[is_tropical_woody_savanna] <- biome_names["Warm Woody Savanna, Woodland & Shrubland"] # nolint

  # open shrubland / savanna
  biome_class[is_temperate_shrubland] <- biome_names["Temperate Savanna & Open Shrubland"] # nolint
  biome_class[is_tropical_shrubland] <- biome_names["Warm Savanna & Open Shrubland"] # nolint

  # grassland
  biome_class[is_temperate_grassland] <- biome_names["Temperate Grassland"]
  biome_class[is_tropical_grassland] <- biome_names["Warm Grassland"]

  biome_class[is_arctic_tundra] <- biome_names["Arctic Tundra"]
  if (!is.na(is_montane_artic_proxy)) {
    biome_class[
      biome_class == biome_names["Arctic Tundra"] & is_montane_artic
    ] <- biome_names["Montane Grassland"]
  }

  # other
  biome_class[is_rocks_and_ice] <- biome_names["Rocks and Ice"]
  biome_class[is_water] <- biome_names["Water"]

  return(list(biome_id = biome_class, biome_names = names(biome_names)))
}


read_pft_categories <- function(file_path) {
  # read_delim, col_types = readr::cols(), delim = ";")to suppress messages
  readr::read_delim(file_path, col_types = readr::cols(), delim = ";") %>%
    # change 1, 0.5, 0 values to TRUE and NAs (NA's can be dropped)
    dplyr::mutate_at(dplyr::vars(dplyr::starts_with(c("category_", "zone_"))),
                     function(x) ifelse(as.logical(x), TRUE, NA)) %>%
    # filter natural pfts
    dplyr::filter(category_natural) %>%
    # all binary zone columns (tropical, temperate, boreal) in one categorical
    #   zone column
    tidyr::pivot_longer(cols = starts_with("zone_"),
                 names_to = "zone",
                 names_prefix = "zone_",
                 values_to = "zone_value",
                 values_drop_na = TRUE) %>%
    # all binary category columns (natural, needle, evergreen) in one categorical # nolint
    #   category column
    tidyr::pivot_longer(cols = starts_with("category_"),
                 names_to = "category",
                 names_prefix = "category_",
                 values_to = "category_value",
                 values_drop_na = TRUE) %>%
    # delete side product - logical columns
    dplyr::select(-c("category_value", "zone_value")) %>%
    # values to lpjml_index, names to length of npft (convert to numeric)
    tidyr::pivot_longer(cols = starts_with("lpjml_index_npft_"),
                 values_to = "lpjml_index",
                 names_to = "npft_proxy",
                 names_transform = list(npft_proxy = function(x) suppressWarnings(as.numeric(x))), # nolint
                 names_prefix = "lpjml_index_npft_",
                 values_drop_na = TRUE) %>%
    return()
}


# bring in line dimnames from lpjmiotools with function
fix_dimnames <- function(data,
                         name,
                         ncell = 67420,
                         npft = NULL,
                         band_names = NULL) {
  if (name == "grid") {
    dim(data) <- c(coordinate = 2, cell = ncell)
    dimnames(data) <- list(coordinate = c("lon", "lat"),
                                 cell = seq_len(ncell))
  } else if (name == "fpc") {
    # process foliage projected cover (fpc)
    dim(data) <- c(cell = ncell,
                   band = npft + 1,
                   year = length(timespan[1]:timespan[2]))
    dimnames(data) <- list(cell = seq_len(ncell),
                           band = fpc_names,
                           year = c(timespan[1]:timespan[2]))
  } else if (name == "vegc" || name == "temp" ) {
      dim(data) <- c(cell = ncell,  year = length(timespan[1]:timespan[2]))
      dimnames(data) <- list(cell = seq_len(ncell),
                             year = c(timespan[1]:timespan[2]))
  } else if (name == "pft_lai") {
      # process pft_lai input
      pre_dims <- dim(pft_lai)
      dim(data) <- c(cell = pre_dims[1],
                     band = pre_dims[2],
                     year = length(timespan[1]:timespan[2]))
      dimnames(data) <- list(cell = seq_len(pre_dims[1]),
                             band = c(fpc_names[2:(npft + 1)],(npft + 1):pre_dims[2]),
                             year = c(timespan[1]:timespan[2]))
  }
  return(data)
}

