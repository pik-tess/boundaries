#' Calculate the planetary boundary status for the land-system change boundary
#'
#' Calculate the PB status for the LSC (land-system change) boundary based
#' on a scenario LPJmL run and a reference LPJmL run.
#'
#' @param path_scenario output directory (character string) of the scenario
#' LPJmL run where binary files (soon with metafiles) are written
#'
#' @param path_reference output directory (character string) of the reference
#' LPJmL run where binary files (soon with metafiles) are written
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#' as an integer vector, e.g. `1982:2011` (default)
#'
#' @param time_span_reference time span to be used for the scenario run, defined
#' as an integer vector, e.g. `1901:1930`. Can differ in offset and length from
#' `time_span_scenario`! If `NULL` value of `time_span_scenario` is used
#'
#' @param spatial_resolution character. Spatial resolution, available options
#' are `"biome"` (default) and `"grid"`
#'
#' @param eurasia logical. If `spatial_resolution` = `"biome"` merge continents
#' Europe and Asia to avoid arbitrary biome cut at europe/asia border. Defaults
#' to `TRUE`
#'
#' @param vegc_proxy logical. If `TRUE` vegetation carbon (vegc) threshold of
#' 7500 gC/m² is used to distinguish between forest and savannahs
#'
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#' time series analysis
#'
#' @examples
#' \dontrun{
#'  calc_lsc_status(path_scenario, path_reference)
#' }
#'
#' @md
#' @export
calc_lsc_status <- function(path_scenario,
                            path_reference,
                            time_span_scenario = c(1982, 2011),
                            time_span_reference = NULL,
                            spatial_resolution = "biome",
                            eurasia = TRUE,
                            vegc_proxy = FALSE,
                            avg_nyear_args = list(),
                            # to be replaced by lpjmlKit::read_output
                            start_year = 1901) {

  # verify available temporal resolution
  spatial_resolution <- match.arg(spatial_resolution, c("biome",
                                                        "grid"))

  if (.Platform$OS.type == "windows") {
    future_plan <- future::plan("multisession")
  } else {
    future_plan <- future::plan("multicore")
  }
  on.exit(future::plan(future_plan))

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
  biome_classes %<-% classify_biomes(
      path_data = path_reference,
      time_span = time_span_reference,
      vegc_proxy = vegc_proxy,
      avg_nyear_args = avg_nyear_args,
      # to be replaced by lpjmlKit::read_output
      start_year = start_year)

  if (spatial_resolution == "biome") {
    # get continents mask - pass arg of whether to merge europe and asia
    continent_grid %<-% calc_continents_mask(path_reference, eurasia = eurasia)
  }

  # TO BE REPLACED BY lpjmlKit::read_output -------------------------------- #
  # read grid
  ncell <- 67420
  size <- 2
  grid_file <- file(paste(path_reference, "grid.bin", sep = "/"), "rb")
  lpjml_grid <- readBin(grid_file, integer(), n = 2 * ncell, size = size) /
                100
  close(grid_file)
  dim(lpjml_grid) <- c(coordinate = 2, cell = ncell)
  dimnames(lpjml_grid) <- list(coordinate = c("lon", "lat"),
                               cell = seq_len(ncell))
  # ------------------------------------------------------------------------ #
  cell_area <-  lpjmlKit::calc_cellarea(
    lpjmlKit::subset_array(lpjml_grid, list(coordinate = "lat"))
  )

  # TO BE REPLACED BY lpjmlKit::read_output ---------------------------------- #
  #   hardcoded values to be internally replaced
  # read foliage projected cover (fpc) output
  ncell <- 67420
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

  fpc_scenario %<-% tmp_read_pft_yearly(
    file_name = paste0(path_scenario,
                       "/",
                       "fpc.bin"),
    time_span = time_span_scenario,
    start_year = start_year,
    ncell = ncell,
    nbands = 12,
    band_names = fpc_names,
    size = 4
  )

  fpc_reference %<-% tmp_read_pft_yearly(
    file_name = paste0(path_reference,
                       "/",
                       "fpc.bin"),
    time_span = time_span_reference,
    start_year = start_year,
    ncell = ncell,
    nbands = 12,
    band_names = fpc_names,
    size = 4
  )
  # subset only tree pfts
  tree_share_scenario <- lpjmlKit::subset_array(
    fpc_scenario,
    list(band = c("Tropical Broadleaved Evergreen Tree", # 2
                  "Tropical Broadleaved Raingreen Tree", # 3
                  "Temperate Needleleaved Evergreen Tree", # 4
                  "Temperate Broadleaved Evergreen Tree", # 5
                  "Temperate Broadleaved Summergreen Tree", # 6
                  "Boreal Needleleaved Evergreen Tree", # 7
                  "Boreal Broadleaved Summergreen Tree", # 8
                  "Boreal Needleleaved Summergreen Tree"))
  )
  tree_share_reference <- lpjmlKit::subset_array(
    fpc_reference,
    list(band = c("Tropical Broadleaved Evergreen Tree", # 2
                  "Tropical Broadleaved Raingreen Tree", # 3
                  "Temperate Needleleaved Evergreen Tree", # 4
                  "Temperate Broadleaved Evergreen Tree", # 5
                  "Temperate Broadleaved Summergreen Tree", # 6
                  "Boreal Needleleaved Evergreen Tree", # 7
                  "Boreal Broadleaved Summergreen Tree", # 8
                  "Boreal Needleleaved Summergreen Tree"))
  )
  # calculate actual tree cover area
  tree_cover_scenario <- (
    tree_share_scenario *
    array(lpjmlKit::subset_array(fpc_scenario,
                    list(band = c("natvegfrac"))),
         dim = dim(tree_share_scenario)) * cell_area
  )
  tree_cover_reference <- (
    tree_share_reference *
    array(lpjmlKit::subset_array(fpc_reference,
                    list(band = c("natvegfrac"))),
         dim = dim(tree_share_reference)) * cell_area
  )
  # sum tree pfts for forest cover
  all_tree_cover_scenario %<-% apply(tree_cover_scenario,
                                     c("cell", "year"),
                                     sum,
                                     na.rm = TRUE)
  all_tree_cover_reference %<-% apply(tree_cover_reference,
                                      c("cell", "year"),
                                      sum,
                                      na.rm = TRUE)
  # average forest over time
  avg_trees_scenario %<-% do.call(average_nyear_window,
                                  append(list(x = all_tree_cover_scenario),
                                         avg_nyear_args))
  if (!is.null(nyear_ref)) {
    avg_nyear_args["nyear_reference"] <- nyear_ref
  }

  # average forest over time
  avg_trees_reference %<-% do.call(average_nyear_window,
                                   append(list(x = all_tree_cover_reference),
                                          avg_nyear_args))

  # binary is forest biome - mask
  is_forest <- array(0,
                     dim = dim(biome_classes),
                     dimnames = dimnames(biome_classes))

  # binary is tropical forest biome - mask
  is_tropical_forest <- is_forest
  # binary is temperate forest biome - mask
  is_temperate_forest <- is_forest
  # binary is boreal forest biome - mask
  is_boreal_forest <- is_forest
  # init forest type mask (tropical:1, temperate:2, boreal:3)
  forest_type <- is_forest

  # 8 forest biomes
  is_forest[biome_classes %in% seq_len(8)] <- 1

  # 2 tropical forest biomes
  is_tropical_forest[biome_classes %in% seq_len(2)] <- 1
  forest_type[which(is_tropical_forest == 1)] <- 1
  # 4 temperate forest biomes
  is_temperate_forest[biome_classes %in% seq(3, 6)] <- 1
  forest_type[which(is_temperate_forest == 1)] <- 2

  # 2 boreal forest biomes
  is_boreal_forest[biome_classes %in% seq(7, 8)] <- 1
  forest_type[which(is_boreal_forest == 1)] <- 3

  # calculate deforestation status by share of scenario and reference run
  deforestation <- avg_trees_scenario / (avg_trees_reference + 1e-9)
  deforestation[deforestation > 1] <- 1
  deforestation[deforestation > 0] <- 1 - deforestation[deforestation > 0]

  # get flexibly named time dimension
  third_dim <- names(dim(deforestation))[
    !names(dim(deforestation)) %in% c("cell")
  ]

  if (spatial_resolution == "biome") {
    # create space of combinations to loop over (even though not all make sense)
    comb <- expand.grid(
      continent = sort(unique(lpjmlKit::subset_array(continent_grid, list(coordinate = "continent")))), # nolint
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
          lpjmlKit::subset_array(continent_grid, list(coordinate = "continent"), drop = FALSE) == comb$continent[idx], # nolint
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
  # boundaries for tropical forest and boreal forest after Steffen et al. 2015
  #   (https://doi.org/10.1126/science.1259855)
  #   share for safe space <= 0.15, uncertainty zone < 0.4 & >0.15 and
  #   for high risk >= 0.4
  # high risk
  pb_status[
    (is_tropical_forest |
    is_boreal_forest) &
    deforestation >= 0.4
  ] <- 3
  # uncertainty zone
  pb_status[
    (is_tropical_forest |
    is_boreal_forest) &
    deforestation < 0.4 &
    deforestation > 0.15
  ] <- 2
  # safe space
  pb_status[
    (is_tropical_forest |
    is_boreal_forest) &
    deforestation <= 0.15
  ] <- 1

  # boundaries for temperate forest after Steffen et al. 2015
  #   (https://doi.org/10.1126/science.1259855)
  #   share for safe space <= 0.5, uncertainty zone < 0.7 & >0.5 and
  #   for high risk >= 0.7
  # high risk
  pb_status[
    is_temperate_forest &
    deforestation >= 0.7
  ] <- 3
  # uncertainty zone
  pb_status[
    is_temperate_forest &
    deforestation < 0.7 &
    deforestation > 0.5
  ] <- 2
  # safe space
  pb_status[
    is_temperate_forest &
    deforestation <= 0.5
  ] <- 1

  return(pb_status)
}