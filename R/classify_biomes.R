#' Classify biomes
#'
#' Classify biomes based on foliage protected cover (FPC) and temperature
#' LPJmL output plus either vegetation carbon or pft_lai depending on
#' the savanna_proxy option and elevation if montane_arctic_proxy requires this
#' information.
#'
#' @param config_reference character string. File path to the LPjmL
#' configuration file (json) of the reference run. The configuration file
#' contains the information about the LPJmL run, e.g. the output
#' directory
#'
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the reference LPJmL run. All needed files are
#' provided as key value pairs, e.g. `list(vegc = "/temp/vegc.bin.json`.
#' If `config_reference` is supplied with all needed files, files
#' reference can be set to `NULL`.
#'
#' @param time_span_reference time span to be used for the classification of
#' biomes, defined as character string, e.g. `as.character(1901:1930)`.
#'
#' @param approach character string indicating which biome classification
#' approach to use. Currently only one is defined ("default").
#'
#' @param time_series_avg integer. Number of years to be used for the moving
#' average calculation. If `NULL`, all years are averaged for one status
#' calculation, for `1` the whole time span is used to calculate a status
#' time series.
#'
#' @param config_args list of arguments to be passed on from the model
#' configuration.
#'
#' @param savanna_proxy `list` with either pft_lai or vegc as
#' key and value in m2/m2 for pft_lai (default: 6) and gC/m2 for
#' vegc (current default: 7500); set to `NULL` if no proxy should be
#' used.
#'
#' @param montane_arctic_proxy `list` with either "elevation" or "latitude" as
#' name/key and value in m for elevation (default: 1000) and degree for
#' latitude (default: 55); set to `NULL` if no proxy is used.
#'
#' @param tree_cover_thresholds list with minimum tree cover thresholds for
#' definition of forest, woodland, savanna and grassland. Only changes to
#' the default have to be included in the list, for the rest the default
#' is used.
#' Default values, based on the IGBP land cover classification system:
#' "boreal forest" = 0.6
#' "temperate forest" = 0.6
#' "temperate woodland" = 0.3
#' "temperate savanna" = 0.1
#' "tropical forest" = 0.6
#' "tropical woodland" = 0.3
#' "tropical savanna" = 0.1
#' In the boreal zone, there is no woodland, everything below the
#' boreal forest threshold will be classified as boreal tundra.
#'
#' @return list object containing biome_id (main biome per grid
#'  `cell[dim=c(ncells)]`), and list of respective
#'  `biome_names[dim=c(nbiomes)]`
#'
#' @examples
#' \dontrun{
#' classify_biomes(
#'   config_reference = "./outputs/config.json",
#'   time_span_reference = 1982:2011
#' )
#' }
#'
#' @md
#' @export
classify_biomes <- function(config_reference = NULL, # nolint:cyclocomp_linter
                            files_reference = NULL,
                            time_span_reference,
                            savanna_proxy = list(vegc = 7500),
                            montane_arctic_proxy = list(elevation = 1000),
                            tree_cover_thresholds = list(),
                            approach = "default",
                            time_series_avg = NULL,
                            config_args = list()) {

  if (is.null(files_reference) && is.null(config_reference)) {
    stop("files_reference or path_reference must be provided")

  } else if (!is.null(config_reference) && is.null(files_reference)) {
    # Get main file type (meta, clm)
    config_reference <- lpjmlkit::read_config(config_reference)

    if (!all(time_span_reference %in% get_sim_time(config_reference))) {
      stop("Time span not available in reference run.")
    }
    # List required output files for each boundary

    output_files <- list_outputs(
      "biome",
      approach = list("biome" = approach),
      spatial_scale = "subglobal",
      only_first_filename = FALSE
    )

    files_reference <- get_filenames(
      config = config_reference,
      output_files = output_files
    )
  }

  # test if provided proxies are valid
  savanna_proxy_name <- match.arg(
    names(savanna_proxy),
    c(NA, "vegc", "pft_lai")
  )
  montane_arctic_proxy_name <- match.arg(names(montane_arctic_proxy),
                                         c(NA, "elevation", "latitude"))

  # define default minimum tree cover for forest / woodland / savanna
  min_tree_cover <- list("boreal forest" = 0.6,
                         "temperate forest" = 0.6,
                         "temperate woodland" = 0.3,
                         "temperate savanna" = 0.1,
                         "tropical forest" = 0.6,
                         "tropical woodland" = 0.3,
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
      min_tree_cover[["temperate woodland"]] | # nolint
      min_tree_cover[["temperate woodland"]] <=
        min_tree_cover[["temperate savanna"]] | # nolint
      min_tree_cover[["tropical woodland"]] <=
        min_tree_cover[["tropical savanna"]] | # nolint
      min_tree_cover[["tropical forest"]] <=
        min_tree_cover[["tropical woodland"]]
  ) {
    stop(paste0("Tree cover threshold for forest are not always higher than",
                "tree cover thresholds for woodland and savanna. Aborting."))
  }

  # Assumptions on fpc and temperature thresholds for biome classification
  # (e.g. distinction between desert and grassland)
  low_fpc <- 0.05
  low_temp_threshold <- 0

  # -------------------------------------------------------------------------- #
  # read in relevant data
  grid <- lpjmlkit::read_io(
    files_reference$grid,
    silent = TRUE
  ) %>%
    conditional_subset(config_args$spatial_subset)

  # bands currently not named in grid file, therefore band 2 manually selected
  lat <- lpjmlkit::as_array(grid, subset = list(band = 2)) %>%
    drop()

  # please R CMD check for use of future operator
  fpc <- temp <- NULL
  fpc %<-% read_io_format(
    files_reference$fpc,
    time_span_reference,
    spatial_subset = config_args$spatial_subset
  )
  temp %<-% {
    lpjmlkit::read_io(
      files_reference$temp,
      subset = list(year = time_span_reference),
      silent = TRUE
    ) %>%
      conditional_subset(config_args$spatial_subset) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      lpjmlkit::as_array(aggregate =
                           list(month = mean, day = mean, band = mean)) %>%
      suppressWarnings()
  }

  if (!is.na(savanna_proxy_name)) {
    # please R CMD check for use of future operator
    savanna_proxy_data <- NULL
    savanna_proxy_data %<-% read_io_format(
      files_reference[[savanna_proxy_name]],
      time_span_reference,
      spatial_subset = config_args$spatial_subset,
      aggregate = list(band = sum)
    )
  }

  if (!is.na(montane_arctic_proxy_name)) {
    if (montane_arctic_proxy_name == "elevation") {
      elevation <- lpjmlkit::read_io(
        files_reference$elevation,
        silent = TRUE
      ) %>%
        conditional_subset(config_args$spatial_subset) %>%
        lpjmlkit::as_array() %>%
        drop()
    }
  }

  fpc_nbands <- dim(fpc)[["band"]]
  npft <- fpc_nbands - 1

  # please R CMD check for use of future operator
  avg_fpc <- NULL
  # average fpc
  avg_fpc %<-% aggregate_time(
    x = fpc,
    time_series_avg = time_series_avg
  )

  # average vegc or pft_lai
  if (!is.na(savanna_proxy_name)) {
    # please R CMD check for use of future operator
    avg_savanna_proxy_data <- NULL
    avg_savanna_proxy_data %<-% aggregate_time(
      x = savanna_proxy_data,
      time_series_avg = time_series_avg
    )
  }
  # please R CMD check for use of future operator
  avg_temp <- NULL
  # average temp
  avg_temp %<-% aggregate_time(
    x = temp,
    time_series_avg = time_series_avg
  )

  # biome_names after biome classification in Ostberg et al. 2013
  # (https://doi.org/10.5194/esd-4-347-2013), Ostberg et al 2015
  # (https://doi.org/10.1088/1748-9326/10/4/044011) and Gerten et al. 2020
  # (https://doi.org/10.1038/s41893-019-0465-1)

  # read in biome names
  biome_mapping <- system.file("extdata",
                               "biomes.csv",
                               package = "boundaries") %>%
    readr::read_delim(col_types = readr::cols(), delim = ";")
  biome_names <- biome_mapping$id
  names(biome_names) <- biome_mapping$name

  # please R CMD check for use of dplyr syntax
  npft_proxy <- zone <- category <- type <- NULL
  # read in pft categories
  pft_categories <- system.file("extdata",
                                "pft_categories.csv",
                                package = "boundaries") %>%
    read_pft_categories() %>%
    dplyr::filter(., npft_proxy == npft)

  # define pft categories
  fpc_temperate_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & zone == "temperate" & category == "natural"
  )$pft

  fpc_tropical_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & zone == "tropical" & category == "natural"
  )$pft

  fpc_boreal_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & zone == "boreal" & category == "natural"
  )$pft

  fpc_needle_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & category == "needle"
  )$pft

  fpc_evergreen_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & category == "evergreen"
  )$pft

  fpc_grass <- dplyr::filter(
    pft_categories,
    type == "grass" & category == "natural"
  )$pft

  fpc_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & category == "natural"
  )$pft

  # sum up fpc over defined pft categories
  third_dim <- names(dim(avg_fpc))[
    !names(dim(avg_fpc)) %in% c("cell", "band") # can be "year"
  ] %>% {
    if (rlang::is_empty(.)) NULL else .
  }
  # please R CMD check for use of future operator
  fpc_tree_total <- NULL
  fpc_tree_total %<-% apply(
    lpjmlkit::asub(avg_fpc, band = fpc_trees, drop = FALSE),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  # please R CMD check for use of future operator
  fpc_tree_tropical <- NULL
  fpc_tree_tropical %<-% apply(
    lpjmlkit::asub(avg_fpc, band = fpc_tropical_trees, drop = FALSE),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  # please R CMD check for use of future operator
  fpc_tree_temperate <- NULL
  fpc_tree_temperate %<-% apply(
    lpjmlkit::asub(avg_fpc, band = fpc_temperate_trees, drop = FALSE),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  # please R CMD check for use of future operator
  fpc_tree_boreal <- NULL
  fpc_tree_boreal %<-% apply(
    lpjmlkit::asub(avg_fpc, band = fpc_boreal_trees, drop = FALSE),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  # please R CMD check for use of future operator
  fpc_tree_needle <- NULL
  fpc_tree_needle %<-% apply(
    lpjmlkit::asub(avg_fpc, band = fpc_needle_trees, drop = FALSE),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  # please R CMD check for use of future operator
  fpc_tree_evergreen <- NULL
  fpc_tree_evergreen %<-% apply(
    lpjmlkit::asub(avg_fpc, band = fpc_evergreen_trees, drop = FALSE),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  # please R CMD check for use of future operator
  fpc_grass_total <- NULL
  fpc_grass_total %<-% apply(
    lpjmlkit::asub(avg_fpc, band = fpc_grass, drop = FALSE),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )
  # please R CMD check for use of future operator
  fpc_total <- NULL
  fpc_total %<-% apply(
    lpjmlkit::asub(avg_fpc, band = -1, drop = FALSE),
    c("cell", third_dim),
    sum,
    na.rm = TRUE
  )

  # please R CMD check for use of future operator
  max_share_trees <- NULL
  # define maximum share of trees among all tree PFTs per grid cell
  max_share_trees %<-% apply(
    lpjmlkit::asub(avg_fpc, band = fpc_trees, drop = FALSE),
    c("cell", third_dim),
    max,
    na.rm = TRUE
  )

  # use vegc 7500 gC/m2 or natLAI 6 as proxy threshold for forest/savanna
  #   "boundary
  if (!is.null(savanna_proxy)) {
    if (savanna_proxy_name == "pft_lai") {
      # please R CMD check for use of future operator
      avg_savanna_proxy_data <- NULL
      avg_savanna_proxy_data %<-% apply(
        lpjmlkit::asub(avg_savanna_proxy_data, band = 1:npft, drop = FALSE) * # nolint
          lpjmlkit::asub(avg_fpc, band = 2: (npft + 1), drop = FALSE) *
          lpjmlkit::asub(avg_fpc, band = 1, drop = FALSE),
        c("cell", third_dim),
        sum
      )
    }
    is_tropical_proxy <- avg_savanna_proxy_data >= savanna_proxy[[savanna_proxy_name]] # nolint
    is_savanna_proxy <- avg_savanna_proxy_data < savanna_proxy[[savanna_proxy_name]] # nolint
  } else {
    is_tropical_proxy <- array(TRUE,
                               dim = dim(avg_temp),
                               dimnames = dimnames(avg_temp))
    is_savanna_proxy <- array(FALSE,
                              dim = dim(avg_temp),
                              dimnames = dimnames(avg_temp))
  }

  # Desert
  is_desert <- {
    fpc_total <= low_fpc &
      avg_temp >= low_temp_threshold
  }

  # montane (for classification of montane grassland)
  if (!is.na(montane_arctic_proxy_name)) {
    if (montane_arctic_proxy_name == "elevation") {
      is_montane_artic <- elevation > montane_arctic_proxy[[
        montane_arctic_proxy_name
      ]]
    } else if (montane_arctic_proxy_name == "latitude") {
      is_montane_artic <- !(abs(lat) > montane_arctic_proxy[[
        montane_arctic_proxy_name
      ]])
    }
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
      (abind::adrop(
        lpjmlkit::asub(
          avg_fpc, band = "boreal needleleaved evergreen tree", drop = FALSE
        ),
        drop = "band"
      ) == max_share_trees)
  }


  if (npft == 9) {
    # Boreal Broadleaved Deciduous
    # no simulation of boreal needleleaved summergreen trees
    is_boreal_broad_deciduous <- {
      is_boreal_forest &
        (abind::adrop(
          lpjmlkit::asub(
            avg_fpc, band = "boreal broadleaved summergreen tree",
            drop = FALSE
          ),
          drop = "band"
        ) == max_share_trees)
    }
  } else {
    # Boreal Deciduous
    is_boreal_broad_deciduous <- {
      is_boreal_forest &
        (abind::adrop(
          lpjmlkit::asub(
            avg_fpc,
            band = "boreal broadleaved summergreen tree",
            drop = FALSE
          ),
          drop = "band"
        ) == max_share_trees)
    }

    is_boreal_needle_deciduous <- {
      is_boreal_forest &
        (abind::adrop(
          lpjmlkit::asub(
            avg_fpc,
            band = "boreal needleleaved summergreen tree",
            drop = FALSE
          ),
          drop = "band"
        ) == max_share_trees)
    }
  }

  # Temperate Coniferous Forest
  is_temperate_coniferous <- {
    is_temperate_forest &
      (abind::adrop(
        lpjmlkit::asub(
          avg_fpc,
          band = "temperate needleleaved evergreen tree",
          drop = FALSE
        ),
        drop = "band"
      ) == max_share_trees)
  }
  # Temperate Broadleaved Evergreen Forest
  is_temperate_broadleaved_evergreen <- { # nolint
    is_temperate_forest &
      (abind::adrop(
        lpjmlkit::asub(
          avg_fpc,
          band = "temperate broadleaved evergreen tree",
          drop = FALSE
        ),
        drop = "band"
      ) == max_share_trees)
  }
  # Temperate Broadleaved Deciduous Forest
  is_temperate_broadleaved_deciduous <- { # nolint
    is_temperate_forest &
      (abind::adrop(
        lpjmlkit::asub(
          avg_fpc,
          band = "temperate broadleaved summergreen tree",
          drop = FALSE
        ),
        drop = "band"
      ) == max_share_trees)
  }

  # Tropical Rainforest
  is_tropical_evergreen <- {
    is_tropical_forest &
      (abind::adrop(
        lpjmlkit::asub(
          avg_fpc,
          band = "tropical broadleaved evergreen tree",
          drop = FALSE
        ),
        drop = "band"
      ) == max_share_trees) &
      is_tropical_proxy
  }

  # Tropical Seasonal & Deciduous Forest
  is_tropical_raingreen <- {
    is_tropical_forest &
      (abind::adrop(
        lpjmlkit::asub(
          avg_fpc,
          band = "tropical broadleaved raingreen tree",
          drop = FALSE
        ),
        drop = "band"
      ) == max_share_trees) &
      is_tropical_proxy
  }
  # Warm Woody Savanna, Woodland & Shrubland
  is_tropical_forest_savanna <- {
    is_tropical_forest &
      (abind::adrop(
        lpjmlkit::asub(
          avg_fpc,
          band = "tropical broadleaved evergreen tree",
          drop = FALSE
        ),
        drop = "band"
      ) == max_share_trees |
        (abind::adrop(
          lpjmlkit::asub(
            avg_fpc,
            band = "tropical broadleaved raingreen tree",
            drop = FALSE
          ),
          drop = "band"
        ) == max_share_trees)
      ) &
      is_savanna_proxy
  }

  # WOODY savanna ----------------------------------------------------------- #

  # Temperate Woody Savanna, Woodland & Shrubland
  is_temperate_woody_savanna <- {
    fpc_tree_total <= min_tree_cover[["temperate forest"]] &
      fpc_tree_total >= min_tree_cover[["temperate woodland"]] &
      abind::adrop(lpjmlkit::asub(avg_fpc, band = "temperate c3 grass",
                                  drop = FALSE), drop = "band") >
        abind::adrop(lpjmlkit::asub(avg_fpc, band = "tropical c4 grass",
                                    drop = FALSE), drop = "band") &
      avg_temp >= low_temp_threshold
  }
  # Warm Woody Savanna, Woodland & Shrubland
  is_tropical_woody_savanna <- {
    fpc_tree_total <= min_tree_cover[["tropical forest"]] &
      fpc_tree_total >= min_tree_cover[["tropical woodland"]] &
      abind::adrop(lpjmlkit::asub(avg_fpc, band = "temperate c3 grass",
                                  drop = FALSE), drop = "band") <
        abind::adrop(lpjmlkit::asub(avg_fpc, band = "tropical c4 grass",
                                    drop = FALSE), drop = "band")
  }

  # OPEN SHRUBLAND / SAVANNAS ----------------------------------------------- #

  # Temperate Savanna & Open Shrubland
  is_temperate_shrubland <- {
    fpc_tree_total <= min_tree_cover[["temperate woodland"]] &
      fpc_tree_total >= min_tree_cover[["temperate savanna"]] &
      abind::adrop(lpjmlkit::asub(avg_fpc, band = "temperate c3 grass",
                                  drop = FALSE), drop = "band") >
        abind::adrop(lpjmlkit::asub(avg_fpc, band = "tropical c4 grass",
                                    drop = FALSE), drop = "band") &
      avg_temp >= low_temp_threshold
  }
  # Warm Savanna & Open Shrubland
  is_tropical_shrubland <- {
    fpc_tree_total <= min_tree_cover[["tropical woodland"]] &
      fpc_tree_total >= min_tree_cover[["tropical savanna"]] &
      abind::adrop(lpjmlkit::asub(avg_fpc, band = "temperate c3 grass",
                                  drop = FALSE), drop = "band") <
        abind::adrop(lpjmlkit::asub(avg_fpc, band = "tropical c4 grass",
                                    drop = FALSE), drop = "band") &
      avg_temp >= low_temp_threshold
  }

  # GRASSLAND ---------------------------------------------------------------- #

  # Temperate grassland
  is_temperate_grassland <- {
    fpc_total > low_fpc &
      fpc_tree_total <= min_tree_cover[["temperate savanna"]] &
      abind::adrop(lpjmlkit::asub(avg_fpc, band = "temperate c3 grass",
                                  drop = FALSE), drop = "band") >
        abind::adrop(lpjmlkit::asub(avg_fpc, band = "tropical c4 grass",
                                    drop = FALSE), drop = "band") &
      avg_temp >= low_temp_threshold
  }
  # Warm grassland
  is_tropical_grassland <- {
    fpc_total > low_fpc &
      fpc_tree_total <= min_tree_cover[["tropical savanna"]] &
      abind::adrop(lpjmlkit::asub(avg_fpc, band = "temperate c3 grass",
                                  drop = FALSE), drop = "band") <
        abind::adrop(lpjmlkit::asub(avg_fpc, band = "tropical c4 grass",
                                    drop = FALSE), drop = "band") &
      avg_temp >= low_temp_threshold
  }

  # Arctic Tundra ------------------------------------------------------------ #
  is_arctic_tundra <- {
    (!is_boreal_forest &
      !is_temperate_forest &
      (
        avg_temp < low_temp_threshold |
          abind::adrop(lpjmlkit::asub(avg_fpc, band = "temperate c3 grass",
                                      drop = FALSE), drop = "band") ==
            abind::adrop(lpjmlkit::asub(avg_fpc, band = "tropical c4 grass",
                                        drop = FALSE), drop = "band")
      ) &
      fpc_total > low_fpc
    ) |
      (avg_temp < low_temp_threshold & fpc_total < low_fpc)
  }

  # Rocks and Ice
  is_rocks_and_ice <- {
    fpc_total == 0 &
      avg_temp < low_temp_threshold
  }
  # Water body
  is_water <- {
    abind::adrop(lpjmlkit::asub(avg_fpc, band = 1, drop = FALSE),
                 drop = "band") == 0
  }

  # CLASSIFY BIOMES ---------------------------------------------------------- #
  # initiate biome_class array
  biome_class <- array(NA,
                       dim = dim(fpc_total),
                       dimnames = dimnames(fpc_total))

  biome_class[is_desert] <- biome_names["Desert"]

  # forests
  biome_class[is_boreal_evergreen] <- biome_names["Boreal Needleleaved Evergreen Forest"] # nolint
  biome_class[is_boreal_broad_deciduous] <- biome_names["Boreal Broadleaved Deciduous Forest"] # nolint
  biome_class[is_boreal_needle_deciduous] <- biome_names["Boreal Needleleaved Deciduous Forest"] # nolint
  biome_class[is_temperate_coniferous] <- biome_names["Temperate Needleleaved Evergreen Forest"] # nolint
  biome_class[is_temperate_broadleaved_evergreen] <- biome_names["Temperate Broadleaved Evergreen Forest"] # nolint
  biome_class[is_temperate_broadleaved_deciduous] <- biome_names["Temperate Broadleaved Deciduous Forest"] # nolint
  biome_class[is_tropical_evergreen] <- biome_names["Tropical Rainforest"]
  biome_class[is_tropical_raingreen] <- biome_names["Tropical Seasonal & Deciduous Forest"] # nolint
  biome_class[is_tropical_forest_savanna] <- biome_names["Warm Woody Savanna, Woodland & Shrubland"] # nolint

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
  if (!is.na(montane_arctic_proxy_name)) {
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
  # please R CMD check for use of dplyr syntax
  category_natural <- NULL
  # read_delim, col_types = readr::cols(), delim = ";")to suppress messages
  readr::read_delim(file_path, col_types = readr::cols(), delim = ";") %>%
    # change 1, 0.5, 0 values to TRUE and NAs (NA's can be dropped)
    dplyr::mutate_at(
      dplyr::vars(tidyselect::starts_with(c("category_", "zone_"))),
      function(x) ifelse(as.logical(x), TRUE, NA)
    ) %>%
    # filter natural pfts
    dplyr::filter(category_natural) %>%
    # all binary zone columns (tropical, temperate, boreal) in one categorical
    #   zone column
    tidyr::pivot_longer(cols = tidyselect::starts_with("zone_"),
                        names_to = "zone",
                        names_prefix = "zone_",
                        values_to = "zone_value",
                        values_drop_na = TRUE) %>%
    # all binary category columns (natural, needle, evergreen) in one
    #   category column
    tidyr::pivot_longer(cols = tidyselect::starts_with("category_"),
                        names_to = "category",
                        names_prefix = "category_",
                        values_to = "category_value",
                        values_drop_na = TRUE) %>%
    # delete side product - logical columns
    dplyr::select(-c("category_value", "zone_value")) %>%
    # values to lpjml_index, names to length of npft (convert to numeric)
    tidyr::pivot_longer(cols = tidyselect::starts_with("lpjml_index_npft_"),
                        values_to = "lpjml_index",
                        names_to = "npft_proxy",
                        names_transform = list(npft_proxy = function(x) suppressWarnings(as.numeric(x))), # nolint
                        names_prefix = "lpjml_index_npft_") %>%
    return()
}
