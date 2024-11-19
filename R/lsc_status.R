#' Status calculation of the land-system change boundary
#'
#' Planetary Boundary status calculation of the LSC (land-system change) boundary based
#' on a scenario LPJmL run and a reference LPJmL run.
#'
#' @param files_scenario list with variable names and corresponding file paths
#' (character string) of the scenario LPJmL run. Handled automatically via
#' [`calc_status()`].
#'
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the reference LPJmL run. Handled automatically via
#' [`calc_status()`].
#'
#' @param spatial_scale character. Spatial resolution, available options
#' are `"regional"` (at the biome level, default), `"global"` and `"grid"`
#'
#' @param time_span_scenario time span to use output from the scenario run,
#' e.g. `1982:2011`.
#'
#' @param time_span_reference time span use output from the reference run,
#' e.g. `1901:1930`.
#'
#' @param approach approach (character string) to be used , currently available
#' approach is `"steffen2015"`
#'
#' @param thresholds list with deforestation thresholds for defining safe,
#' increasing risk and high risk zone. Default based on Steffen et al. 2015 if
#' thresholds set to NULL (https://doi.org/10.1126/science.1259855):
#' for gridded and biome scale application:
#' `list(pb = list(temperate = 0.5, tropical = 0.15, boreal = 0.15),
#'      highrisk = list(temperate = 0.7, tropical = 0.4, boreal = 0.4))
#' for global scale application:
#' list(holocene = 0, pb = 0.25, highrisk = 0.46)
#' pb = threshold between safe zone and increasing risk zone
#'      (e.g. 50% for boreal forest with default value)
#' highrisk = threshold between increasing risk and high risk zone
#'
#' @param time_series_avg integer. Number of years to be used for the moving
#' average calculation. If `NULL`, all years are averaged for one status
#' calculation, for `1` the whole time span is used to calculate a status time
#' series.
#'
#' @param config_args list of arguments to be passed on from the model
#' configuration.
#'
#' @param eurasia logical. If `spatial_scale` = `"regional"` merge continents
#' Europe and Asia to avoid arbitrary biome cut at europe/asia border.
#' Defaults to `TRUE`
#'
#' @param ... arguments forwarded to [`classify_biomes`]
#'
#'@return Object of class `control_variable` with the boundary status of the
#' lsc boundary.
#'
#' @examples
#' \dontrun{
#' boundary_status <- calc_status(
#'   boundary = "lsc",
#'   config_scenario = "path/to/config_scenario.json",
#'   config_reference = "path/to/config_reference.json",
#'   spatial_scale = "global",
#'   time_span_scenario = 1901:2019,
#'   time_span_reference = 1901:1930
#' )
#' }
#'
#' @md
#' @export
lsc_status <- function(
  files_scenario,
  files_reference,
  spatial_scale = "regional",
  time_span_scenario = as.character(1982:2011),
  time_span_reference = time_span_scenario,
  approach = "steffen2015",
  time_series_avg = NULL,
  config_args = list(),
  thresholds = NULL,
  eurasia = TRUE,
  ...
) {


  # Filter out approach and thresholds arguments from ellipsis
  ellipsis_filtered <- list(...)
  ellipsis_filtered$approach <- NULL
  ellipsis_filtered$thresholds <- NULL

  # please R CMD check for use of future operator
  biome_classes <- NULL
  # please R CMD check for dplyr syntax
  npft_proxy <- category <- type <- NULL

  # classify biomes based on foliage projected cover (FPC) output
  biome_classes %<-% do.call(
    classify_biomes,
    append(list(files_reference = files_reference,
                time_span_reference = time_span_reference,
                time_series_avg = NULL,
                config_args = config_args,
                montane_arctic_proxy = NULL),
           ellipsis_filtered)
  )

  if (spatial_scale == "regional") {
    # please R CMD check for use of future operator
    continent_grid <- NULL
    # get continents mask - pass arg of whether to merge europe and asia
    continent_grid %<-% mask_continents(
      files_reference$grid,
      eurasia = eurasia
    )
  }

  # read in biome mapping
  biome_mapping <- read_biome_mapping(
    system.file("extdata", "biomes.csv", package = "boundaries")
  )

  # please R CMD check for use of future operator
  terr_area <- NULL
  # calculate cell area
  terr_area %<-% {
    lpjmlkit::read_io(
      files_reference$terr_area,
      silent = TRUE
    )$data %>% drop()
  }

  # please R CMD check for use of future operator
  fpc_scenario <- fpc_reference <- NULL
  # read fpc
  fpc_scenario %<-% {
    lpjmlkit::read_io(
      files_scenario$fpc,
      subset = list(year = time_span_scenario),
      silent = TRUE
    ) %>%
      conditional_subset(config_args$spatial_subset) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      lpjmlkit::as_array()
  }

  fpc_reference %<-% {
    lpjmlkit::read_io(
      files_reference$fpc,
      subset = list(year = time_span_reference),
      silent = TRUE
    ) %>%
      conditional_subset(config_args$spatial_subset) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      lpjmlkit::as_array()
  }

  fpc_nbands <- dim(fpc_scenario)[["band"]]
  npft <- fpc_nbands - 1

  # read in pft mapping
  pft_categories <- system.file("extdata",
                                "pft_categories.csv",
                                package = "boundaries") %>%
    read_pft_categories() %>%
    dplyr::filter(., `npft_proxy` == npft)

  # add band names for clm outputs
  if (lpjmlkit::detect_io_type(files_reference$fpc) == "clm") {
    fpc_names <- dplyr::filter(pft_categories, `category` == "natural")$pft
    dimnames(fpc_scenario)$band <- c("natural stand fraction", fpc_names)
    dimnames(fpc_reference)$band <- c("natural stand fraction", fpc_names)
  }

  # subset only tree pfts
  fpc_trees <- dplyr::filter(
    pft_categories,
    `type` == "tree" & `category` == "natural"
  )$pft

  tree_share_scenario <- lpjmlkit::asub(
    fpc_scenario,
    band = fpc_trees,
    drop = FALSE
  )

  tree_share_reference <- lpjmlkit::asub(
    fpc_reference,
    band = fpc_trees,
    drop = FALSE
  )

  # calculate actual tree cover area
  tree_cover_scenario <- (
    tree_share_scenario *
      array(lpjmlkit::asub(fpc_scenario,
                           band = c("natural stand fraction")),
            dim = dim(tree_share_scenario)) * terr_area
  )

  tree_cover_reference <- (
    tree_share_reference *
      array(lpjmlkit::asub(fpc_reference,
                           band = c("natural stand fraction")),
            dim = dim(tree_share_reference)) * terr_area
  )

  # please R CMD check for use of future operator
  all_tree_cover_scenario <- all_tree_cover_reference <- NULL
  # sum tree pfts for forest cover
  all_tree_cover_scenario %<-% apply(
    tree_cover_scenario,
    c("cell", "year"),
    sum,
    na.rm = TRUE
  )

  all_tree_cover_reference %<-% apply(
    tree_cover_reference,
    c("cell", "year"),
    sum,
    na.rm = TRUE
  )

  # please R CMD check for use of future operator
  avg_trees_scenario <- avg_trees_reference <- NULL
  # average forest over time
  avg_trees_scenario %<-% aggregate_time(
    x = all_tree_cover_scenario,
    time_series_avg = time_series_avg
  )

  # average forest over time

  avg_trees_reference %<-% apply(
    all_tree_cover_reference,
    "cell",
    mean
  )

  # binary is forest biome - mask
  is_forest <- array(0,
                     dim = dim(biome_classes$biome_id),
                     dimnames = dimnames(biome_classes$biome_id))

  # binary is tropical forest biome - mask
  is_tropical_forest <- is_forest
  # binary is temperate forest biome - mask
  is_temperate_forest <- is_forest
  # binary is boreal forest biome - mask
  is_boreal_forest <- is_forest
  # init forest type mask (tropical:1, temperate:2, boreal:3)
  forest_type <- is_forest

  # forest biomes
  is_forest[
    biome_classes$biome_id %in% dplyr::filter(
      biome_mapping, `category` == "forest"
    )$id
  ] <- 1

  # please R CMD check for use of dplyr syntax
  zone <- NULL
  # tropical forest biomes
  is_tropical_forest[
    biome_classes$biome_id %in% dplyr::filter(
      biome_mapping, `category` == "forest" & `zone` == "tropical"
    )$id
  ] <- 1
  forest_type[which(is_tropical_forest == 1)] <- 1

  # temperate forest biomes
  is_temperate_forest[
    biome_classes$biome_id %in% dplyr::filter(
      biome_mapping, `category` == "forest" & `zone` == "temperate"
    )$id
  ] <- 1
  forest_type[which(is_temperate_forest == 1)] <- 2

  # boreal forest biomes
  is_boreal_forest[
    biome_classes$biome_id %in% dplyr::filter(
      biome_mapping, `category` == "forest" & `zone` == "boreal"
    )$id
  ] <- 1
  forest_type[which(is_boreal_forest == 1)] <- 3

  # calculate deforestation status by share of scenario and reference run
  deforestation <- avg_trees_scenario / (avg_trees_reference + 1e-9)
  deforestation[deforestation > 1] <- 1
  deforestation[deforestation > 0] <- 1 - deforestation[deforestation > 0]
  deforestation[!is_forest] <- NA
  deforestation <- deforestation * 100 # in percent

  if (spatial_scale == "regional") {

    # create space of combinations to loop over (even though not all make sense)
    comb <- expand.grid(
      continent = sort(
        unique(lpjmlkit::asub(continent_grid, band = "continent"))
      ),
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
      # get flexibly named time dimension
      third_dim <- names(dim(deforestation))[
        !names(dim(deforestation)) %in% c("cell")
      ]
      sub_deforest[!sub_cells] <- NA
      sub_deforest <- apply(
        sub_deforest,
        third_dim,
        function(x) {
          # mean over cell subset of forest type and continent
          # weighted by cell area
          return(
            rep(stats::weighted.mean(x, terr_area, na.rm = TRUE), length(x))
          )
        }
      )
      deforestation[sub_cells] <- sub_deforest[sub_cells]
    }
  }

  if (spatial_scale %in% c("grid", "regional")) {
    # init threshold array with NA = no data and initial dimensions + threshold
    # dimension

    assign_thresholds <- function(forest_type, thresholds) {
      threshold_array <- array(
        NA,
        dim = dim(deforestation),
        dimnames = dimnames(deforestation)
      )
      threshold_array[forest_type == 1] <- thresholds$tropical
      threshold_array[forest_type == 2] <- thresholds$temperate
      threshold_array[forest_type == 3] <- thresholds$boreal
      return(threshold_array)
    }

    pb_thresholds <- assign_thresholds(forest_type, thresholds$pb)
    highrisk_thresholds <- assign_thresholds(forest_type, thresholds$highrisk)
    holocene_thresholds <- assign_thresholds(forest_type, thresholds$holocene)

    dim_names <- dimnames(deforestation)
    dim_names$thresholds <- names(thresholds)
    threshold_attr <- array(
      NA,
      dim = c(dim(deforestation), length(names(thresholds))),
      dimnames = dim_names
    )
    threshold_attr[, , "pb"] <- pb_thresholds
    threshold_attr[, , "highrisk"] <- highrisk_thresholds
    threshold_attr[, , "holocene"] <- holocene_thresholds

    control_variable <- deforestation
    attr(control_variable, "thresholds") <- threshold_attr

  } else if (spatial_scale == "global") {
    dim_remain <- names(dim(deforestation))[names(dim(deforestation)) != "cell"]
    deforestation[as.vector(is_forest == 0), ] <- NA
    control_variable <- apply(deforestation, dim_remain, mean, na.rm = TRUE)
    attr(control_variable, "thresholds") <- thresholds

  }

  attr(control_variable, "control_variable") <- "deforestation"
  attr(control_variable, "spatial_scale") <- spatial_scale
  attr(control_variable, "unit") <- list_unit("lsc", approach, spatial_scale)
  attr(control_variable, "long_name") <- list_long_name("lsc")
  class(control_variable) <- c("control_variable")
  return(control_variable)
}


read_biome_mapping <- function(file_path) {
  # read_delim, col_types = readr::cols(), delim = ";")to suppress messages
  readr::read_delim(file_path, col_types = readr::cols(), delim = ";") %>%
    # change 1, 0.5, 0 values to TRUE and NAs (NA's can be dropped)
    dplyr::mutate_at(
      dplyr::vars(tidyselect::starts_with(c("category_", "zone_"))),
      function(x) ifelse(as.logical(x), TRUE, NA)
    ) %>%
    # all binary zone columns (tropical, temperate, boreal) in one categorical
    #   zone column
    tidyr::pivot_longer(
      cols = tidyselect::starts_with("zone_"),
      names_to = "zone",
      names_prefix = "zone_",
      values_to = "zone_value",
      values_drop_na = TRUE
    ) %>%
    # all binary category columns (forest, xx) in one categorical
    #   category column
    tidyr::pivot_longer(
      cols = tidyselect::starts_with("category_"),
      names_to = "category",
      names_prefix = "category_",
      values_to = "category_value",
      values_drop_na = TRUE
    ) %>%
    # delete side product - logical columns
    dplyr::select(-c("category_value", "zone_value")) %>%
    return()
}
