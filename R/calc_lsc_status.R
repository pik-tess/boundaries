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
#' as a character string, e.g. `as.character(1982:2011)` (default)
#'
#' @param time_span_reference time span to be used for the scenario run, defined
#' as an integer vector, e.g. `as.character(1901:1930)`. Can differ in offset
#' and length from `time_span_scenario`! If `NULL` value of `time_span_scenario`
#' is used
#'
#' @param spatial_resolution character. Spatial resolution, available options
#'        are `"subglobal"` (at the biome level, default), `"global"` and
#'        `"grid"`
#'
#' @param eurasia logical. If `spatial_resolution` = `"biome"` merge continents
#'        Europe and Asia to avoid arbitrary biome cut at europe/asia border.
#'        Defaults to `TRUE`
#'
#' @param thresholds list with deforestation thresholds for defining safe,
#' increasing risk and high risk zone for temperate/temperate/boreal forest
#' biomes, Default based on Steffen et al. 2015
#' (https://doi.org/10.1126/science.1259855):
#' list(boreal = list(pb = 0.5, highrisk = 0.7),
#'      temperate = list(pb = 0.15, highrisk = 0.4),
#'      tropical = list(pb = 0.15, highrisk = 0.4))
#' pb = threshold between safe zone and increasing risk zone (e.g. 50% for
#'         boreal forest with default value)
#' highrisk = threshold between increasing risk and high risk zone
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
                            time_span_scenario = as.character(1982:2011),
                            time_span_reference = NULL,
                            spatial_resolution = "subglobal",
                            eurasia = TRUE,
                            thresholds = list(holocene = list(temperate = 0,
                                                              tropical = 0,
                                                              boreal = 0),
                                              pb = list(temperate = 0.5,
                                                        tropical = 0.15,
                                                        boreal = 0.15),
                                              highrisk = list(temperate = 0.7,
                                                              tropical = 0.4,
                                                              boreal = 0.4)),
                            avg_nyear_args = list(),
                            ...
                            ) {




  # Filter out method and thresholds arguments from ellipsis
  ellipsis_filtered <- list(...)
  ellipsis_filtered$method <- NULL
  ellipsis_filtered$thresholds <- NULL

  # classify biomes based on foliage projected cover (FPC) output
  do.call(classify_biomes,
          append(list(files_reference = files_reference,
                      time_span_reference = time_span_reference,
                      avg_nyear_args = avg_nyear_args),
                 ellipsis_filtered))

  biome_classes <- biome_classes$biome_id

  if (spatial_resolution == "subglobal") {
    # get continents mask - pass arg of whether to merge europe and asia
    continent_grid <- calc_continents_mask(files_reference$grid,
                                           eurasia = eurasia)
  }

  # read in biome mapping
  biome_mapping <- read_biome_mapping(system.file("extdata",
                               "biomes.csv",
                               package = "boundaries"))

  # calculate cell area
  cell_area <- lpjmlkit::read_io(
      files_reference$grid,
      silent = TRUE
      ) %>% 
      lpjmlkit::calc_cellarea()

  # read fpc
  fpc_scenario <- lpjmlkit::read_io(
      files_scenario$fpc,
      subset = list(year = time_span_scenario),
      silent = TRUE
      ) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      as_array()
  fpc_reference <- lpjmlkit::read_io(
      files_reference$fpc,
      subset = list(year = time_span_reference),
      silent = TRUE
      ) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      as_array()

  fpc_nbands <- dim(fpc_scenario)[["band"]]
  npft <- fpc_nbands - 1

  # read in pft mapping
  pft_categories <- system.file("extdata",
                                "pft_categories.csv",
                                package = "boundaries") %>%
    read_pft_categories() %>%
    dplyr::filter(., npft_proxy == npft)

  # add band names for clm outputs
  if(lpjmlkit::detect_io_type(files_reference$fpc) == "clm") {
    fpc_names <- dplyr::filter(pft_categories, category == "natural")$pft
    dimnames(fpc_scenario)$band <- c("natural stand fraction", fpc_names)
    dimnames(fpc_reference)$band <- c("natural stand fraction", fpc_names)
  }

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


  if (length(time_span_reference) < length(time_span_scenario)) {
    avg_nyear_args["nyear_reference"] <- length(time_span_scenario)
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

  if (spatial_resolution == "subglobal") {
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

  if (spatial_resolution %in% c("grid", "subglobal")) {
    # init threshold array with NA = no data and initial dimensions + threshold
    # dimension
    dim_names <- dimnames(deforestation)
    dim_names$thresholds <- names(thresholds)
    threshold_attr <- array(NA,
                       dim = c(dim(deforestation), length(names(thresholds))),
                       dimnames = dim_names)
    ## tropical forest
    # high risk
    threshold_attr[
      is_tropical_forest, "highrisk"
    ] <- thresholds$highrisk$tropical
    # pb
    threshold_attr[
      is_tropical_forest, "pb"
    ] <- thresholds$pb$tropical
    # holocene
    threshold_attr[
      is_tropical_forest, "holocene"
    ] <- thresholds$holocene$tropical

    ## temperate forest
    # high risk
      # high risk
    threshold_attr[
      is_temperate_forest, "highrisk"
    ] <- thresholds$highrisk$temperate
    # pb
    threshold_attr[
      is_temperate_forest, "pb"
    ] <- thresholds$pb$temperate
    # holocene
    threshold_attr[
      is_temperate_forest, "holocene"
    ] <- thresholds$holocene$temperate

    ## boreal forest
    # high risk
    threshold_attr[
      is_boreal_forest, "highrisk"
    ] <- thresholds$highrisk$boreal
    # pb
    threshold_attr[
      is_boreal_forest, "pb"
    ] <- thresholds$pb$boreal
    # holocene
    threshold_attr[
      is_boreal_forest, "holocene"
    ] <- thresholds$holocene$boreal

    attr(deforestation, "thresholds") <- threshold_attr

  } else if (spatial_resolution == "global") {
    #TODO: thresholds to be read in from yml file if not explicitely defined

    dim_remain <- names(dim(deforestation))[names(dim(deforestation)) != "cell"]
    deforestation <- apply(deforestation[is_forest], dim_remain, mean)
    attr(deforestation, "thresholds") <- thresholds

  }

  return(deforestation)
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