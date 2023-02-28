#' Calculate the planetary boundary status for the land-system change boundary
#'
#' Calculate the PB status for the LSC (land-system change) boundary based
#' on a scenario LPJmL run and a reference LPJmL run.
#'
#' @param files_scenario list with variable names and corresponding file paths
#' (character string) of the scenario LPJmL run. All needed files are
#' provided in XXX. E.g.: list(leaching = "/temp/leaching.bin.json")
#'
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the reference LPJmL run. All needed files are
#' provided in XXX. E.g.: list(leaching = "/temp/leaching.bin.json"). If not
#' needed for the applied method, set to NULL.
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#' as an integer vector, e.g. `1982:2011` (default)
#'
#' @param time_span_reference time span to be used for the scenario run, defined
#' as an integer vector, e.g. `1901:1930`. Can differ in offset and length from
#' `time_span_scenario`! If `NULL` value of `time_span_scenario` is used
#'
#' @param lsc_spatial_resolution character. Spatial resolution, available options
#'        are `"biome"` (default) and `"grid"`
#'
#' @param eurasia logical. If `lsc_spatial_resolution` = `"biome"` merge continents
#'        Europe and Asia to avoid arbitrary biome cut at europe/asia border.
#'        Defaults to `TRUE`
#'
#' @param lsc_thresholds list with deforestation thresholds for defining safe,
#' increasing risk and high risk zone for temperate/temperate/boreal forest
#' biomes, Default based on Steffen et al. 2015
#' (https://doi.org/10.1126/science.1259855):
#' list(boreal = list(lower = 0.5, upper = 0.7),
#'      temperate = list(lower = 0.15, upper = 0.4),
#'      tropical = list(lower = 0.15, upper = 0.4))
#' lower = threshold between safe zone and increasing risk zone (e.g. 50% for
#'         boreal forest with default value)
#' upper = threshold between increasing risk and high risk zone
#'
#' @param avg_nyear_args list of arguments to be passed to
#'        \link[pbfunctions]{average_nyear_window} (see for more info).
#'        To be used for time series analysis
#'
#' @param ... arguments forwarded to \link[boundaries](classify_biomes)
#'
#' @examples
#' \dontrun{
#'  calc_lsc_status(files_scenario, files_reference)
#' }
#'
#' @md
#' @export
calc_lsc_status <- function(files_scenario,
                            files_reference,
                            time_span_scenario = c(1982:2011),
                            time_span_reference = NULL,
                            lsc_spatial_resolution = "biome",
                            eurasia = TRUE,
                            lsc_thresholds = list(
                              temperate = list(lower = 0.5, upper = 0.7),
                              boreal = list(lower = 0.15, upper = 0.4),
                              tropical = list(lower = 0.15, upper = 0.4)),
                            avg_nyear_args = list(),
                            ...
                            ) {

  # verify available temporal resolution
  lsc_spatial_resolution <- match.arg(lsc_spatial_resolution, c("biome",
                                                        "grid"))

  # check time_spans of scenario and reference runs
  if (is.null(time_span_reference)) {
    time_span_reference <- time_span_scenario
    nyear_ref <- NULL
  } else {
    if (diff(time_span_reference) > diff(time_span_scenario)) {
      stop(paste0("time_span_reference is longer than time_span_scenario.",
                  "Define a time_span_reference that is shorter than",
                  "time_span_scenario"))
    } else if (diff(time_span_reference) < diff(time_span_scenario)) {
      nyear_ref <- length(time_span_scenario[1]:time_span_scenario[2])
    } else {
      nyear_ref <- NULL
    }
  }

  # classify biomes based on foliage projected cover (FPC) output
  biome_classes <- classify_biomes(
    files_reference = files_reference,
    time_span_reference = time_span_reference,
    avg_nyear_args = avg_nyear_args,
    montane_arctic_proxy = NULL,
    ...
    )
  biome_classes <- biome_classes$biome_id

  if (lsc_spatial_resolution == "biome") {
    # get continents mask - pass arg of whether to merge europe and asia
    continent_grid <- calc_continents_mask(files_reference$grid, eurasia = eurasia)
  }

  # read in biome mapping
  biome_mapping <- read_biome_mapping(system.file("extdata",
                               "biomes.csv",
                               package = "boundaries"))

  # read grid
  grid <- lpjmlkit::read_io(
      files_reference$grid,
      silent = TRUE
      )$data %>% drop()

  # calculate cell area
  cell_area <- lpjmlkit::calc_cellarea(grid[, 2])

  # read fpc
  fpc_scenario <- lpjmlkit::read_io(
      files_scenario$fpc,
      subset = list(year = as.character(time_span_scenario)),
      silent = TRUE
      ) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      as_array()
  fpc_reference <- lpjmlkit::read_io(
      files_reference$fpc,
      subset = list(year = as.character(time_span_reference)),
      silent = TRUE
      ) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      as_array()

  fpc_nbands <- dim(fpc_scenario)[["band"]]
  npft <- fpc_nbands - 1

  pft_categories <- system.file("extdata",
                                "pft_categories.csv",
                                package = "boundaries") %>%
    read_pft_categories() %>%
    dplyr::filter(., npft_proxy == npft)

  fpc_names <- dplyr::filter(pft_categories, category == "natural")$pft

  # only needed for header files: 
  dimnames(fpc_scenario)$band <- c("natural stand fraction", fpc_names)
  dimnames(fpc_reference)$band <- c("natural stand fraction", fpc_names)

  # subset only tree pfts
  fpc_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & category == "natural"
  )$pft

  tree_share_scenario <- lpjmlkit::asub(
    fpc_scenario,
    band = fpc_trees
  )
  tree_share_reference <- lpjmlkit::asub(
    fpc_reference,
    band = fpc_trees
  )
  # calculate actual tree cover area
  tree_cover_scenario <- (
    tree_share_scenario *
    array(lpjmlkit::asub(fpc_scenario,
                         band = c("natural stand fraction")),
         dim = dim(tree_share_scenario)) * cell_area
  )
  tree_cover_reference <- (
    tree_share_reference *
    array(lpjmlkit::asub(fpc_reference,
                         band = c("natural stand fraction")),
         dim = dim(tree_share_reference)) * cell_area
  )
  # sum tree pfts for forest cover
  all_tree_cover_scenario <- apply(tree_cover_scenario,
                                     c("cell", "year"),
                                     sum,
                                     na.rm = TRUE)
  all_tree_cover_reference <- apply(tree_cover_reference,
                                      c("cell", "year"),
                                      sum,
                                      na.rm = TRUE)
  # average forest over time
  avg_trees_scenario <- do.call(average_nyear_window,
                                  append(list(x = all_tree_cover_scenario),
                                         avg_nyear_args))
  if (!is.null(nyear_ref)) {
    avg_nyear_args["nyear_reference"] <- nyear_ref
  }

  # average forest over time
  avg_trees_reference <- do.call(average_nyear_window,
                                   append(list(x = all_tree_cover_reference),
                                          avg_nyear_args))

  # binary is forest biome - mask
  is_forest <- array(0,
                     dim = dim(avg_trees_reference),
                     dimnames = dimnames(avg_trees_reference))

  # binary is tropical forest biome - mask
  is_tropical_forest <- is_forest
  # binary is temperate forest biome - mask
  is_temperate_forest <- is_forest
  # binary is boreal forest biome - mask
  is_boreal_forest <- is_forest
  # init forest type mask (tropical:1, temperate:2, boreal:3)
  forest_type <- is_forest

  # 8 forest biomes
  is_forest[biome_classes %in% dplyr::filter(biome_mapping,
                                             category == "forest"
                                            )$id] <- 1

  # 2 tropical forest biomes
  is_tropical_forest[biome_classes %in% dplyr::filter(biome_mapping,
                                             category == "forest" &
                                             zone == "tropical"
                                            )$id] <- 1
  forest_type[which(is_tropical_forest == 1)] <- 1
  # 4 temperate forest biomes
  is_temperate_forest[biome_classes %in% dplyr::filter(biome_mapping,
                                             category == "forest" &
                                             zone == "temperate"
                                            )$id] <- 1
  forest_type[which(is_temperate_forest == 1)] <- 2

  # 2 boreal forest biomes
  is_boreal_forest[biome_classes %in% dplyr::filter(biome_mapping,
                                             category == "forest" &
                                             zone == "boreal"
                                            )$id] <- 1
  forest_type[which(is_boreal_forest == 1)] <- 3

  # calculate deforestation status by share of scenario and reference run
  deforestation <- avg_trees_scenario / (avg_trees_reference + 1e-9)
  deforestation[deforestation > 1] <- 1
  deforestation[deforestation > 0] <- 1 - deforestation[deforestation > 0]

  # get flexibly named time dimension
  third_dim <- names(dim(deforestation))[
    !names(dim(deforestation)) %in% c("cell")
  ]

  if (lsc_spatial_resolution == "biome") {
    # create space of combinations to loop over (even though not all make sense)
    comb <- expand.grid(
      continent = sort(unique(lpjmlkit::asub(continent_grid,
                                             band = "continent"))),
      forest = sort(unique(factor(forest_type)))[-1]
    )
    for (idx in seq_len(nrow(comb))) {
      sub_deforest <- deforestation
      # replace for every combination a subset of cells with mean of each
      sub_cells <- {
        # match forest type for cells
        forest_type == comb$forest[idx] &
        # match continent for cells
        array(
          lpjmlkit::asub(continent_grid, band = "continent", drop = FALSE) == comb$continent[idx], # nolint
          dim = dim(deforestation),
          dimnames = dimnames(deforestation)
        )
      }
      if (!any(sub_cells)) {
        next
      }
      sub_deforest[!sub_cells] <- NA
      sub_deforest <- apply(
        sub_deforest,
        third_dim,
        function(x) {
          # mean over cell subset of forest type and continent
          # weighted by cell area
          return(rep(weighted.mean(x, cell_area, na.rm = TRUE), length(x)))
        }
      )
      deforestation[sub_cells] <- sub_deforest[sub_cells]
    }
  }
  # init pb_status array with 0 = no data and initial dimensions
  pb_status <- array(0,
                     dim = dim(deforestation),
                     dimnames = dimnames(deforestation))
  ## tropical forest
  # high risk
  pb_status[
    is_tropical_forest &
    deforestation >= lsc_thresholds$tropical$upper
  ] <- 3
  # uncertainty zone
  pb_status[
    is_tropical_forest &
    deforestation < lsc_thresholds$tropical$upper &
    deforestation > lsc_thresholds$tropical$lower
  ] <- 2
  # safe space
  pb_status[
    is_tropical_forest &
    deforestation <= lsc_thresholds$tropical$lower
  ] <- 1

  ## temperate forest
  # high risk
  pb_status[
    is_temperate_forest &
    deforestation >= lsc_thresholds$temperate$upper
  ] <- 3
  # uncertainty zone
  pb_status[
    is_temperate_forest &
    deforestation < lsc_thresholds$temperate$upper &
    deforestation > lsc_thresholds$temperate$lower
  ] <- 2
  # safe space
  pb_status[
    is_temperate_forest &
    deforestation <= lsc_thresholds$temperate$lower
  ] <- 1

  ## boreal forest
  # high risk
  pb_status[
    is_boreal_forest &
    deforestation >= lsc_thresholds$boreal$upper
  ] <- 3
  # uncertainty zone
  pb_status[
    is_boreal_forest &
    deforestation < lsc_thresholds$boreal$upper &
    deforestation > lsc_thresholds$boreal$lower
  ] <- 2
  # safe space
  pb_status[
    is_boreal_forest &
    deforestation <= lsc_thresholds$boreal$lower
  ] <- 1

  return(pb_status)
}

read_biome_mapping <- function(file_path) {
  # read_delim, col_types = readr::cols(), delim = ";")to suppress messages
  readr::read_delim(file_path, col_types = readr::cols(), delim = ";") %>%
    # change 1, 0.5, 0 values to TRUE and NAs (NA's can be dropped)
    dplyr::mutate_at(dplyr::vars(dplyr::starts_with(c("category_", "zone_"))),
                     function(x) ifelse(as.logical(x), TRUE, NA)) %>%
    # all binary zone columns (tropical, temperate, boreal) in one categorical
    #   zone column
    tidyr::pivot_longer(cols = starts_with("zone_"),
                 names_to = "zone",
                 names_prefix = "zone_",
                 values_to = "zone_value",
                 values_drop_na = TRUE) %>%
    # all binary category columns (forest, xx) in one categorical
    #   category column
    tidyr::pivot_longer(cols = starts_with("category_"),
                 names_to = "category",
                 names_prefix = "category_",
                 values_to = "category_value",
                 values_drop_na = TRUE) %>%
    # delete side product - logical columns
    dplyr::select(-c("category_value", "zone_value")) %>%
    return()
}