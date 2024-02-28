#' Calculate biosphere status based on BioCol (HANPP) from a PNV run (reference)
#' and LU run (scenario) of LPJmL, both using time_span_scenario. Additionally
#' a separate reference NPP file (e.g. from a Holocene run) can be supplied as
#' reference_npp_file, which will use time_span_reference, or file index years
#' 3:32 if time_span_reference is not supplied
#'
#' Calculate biosphere status based on BioCol from a PNV run and LU run of LPJmL
#'
#' @param files_scenario list with variable names and corresponding file paths
#' (character string) of the scenario LPJmL run. All needed files need to be
#' provided. E.g.: list(grid = "/temp/grid.bin.json",
#'                      npp = "/temp/npp.bin.json")
#'
<<<<<<< HEAD
=======
#' @param path_baseline character string with path to outputs for the baseline
#' run, file names are taken from files scenario.
#'
>>>>>>> development
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the reference NPP, HANPP should be compared against.
#' In this case only NPP is required. list(npp = "/temp/npp.bin.json").
#'
#' @param time_span_scenario time span to be used for the scenario run and
#' parallel PNV run, defined as a character string,
#' e.g. `as.character(1982:2011)`
#'
#' @param time_span_reference time span to read reference_npp_file from, using
#' index years 3:32 if set to NULL (default: NULL)
#'
#' @param spatial_scale character string indicating spatial resolution
#' either "grid", "subglobal" or "global"
#'
#' @param thresholds named character string with thresholds to be used to
#' define the lower end of safe, increasing risk and high risk zone,
#' e.g. c(holocene = 0.0, pb = 0.1, highrisk = 0.2). If set to NULL, default
#' values from metric_files.yml will be used.
#'
#' @param time_aggregation_args list of arguments to be passed to
#' \link[boundaries]{aggregate_time} (see for more info).
#' To be used for time series analysis
#'
#' @param config_args list of arguments to be passed on from the model
#' configuration.
#'
#' @param path_baseline character string with path to outputs for the baseline
#' run, file names are taken from files scenario.
#'
#' @param time_span_baseline time span to be used for the baseline run, defined
#' as a character vector, e.g. `as.character(1901:1930)`. Can differ in offset
#' and length from `time_span_scenario`! If `NULL` value of `time_span_scenario`
#' is used
#'
#' @param gridbased logical; are pft outputs from LPJmL gridbased or pft-based?
#'
#' @param npp_threshold lower threshold for npp (to mask out non-lu areas
#' according to Haberl et al. 2007). Below BioCol will be set to 0.
#' (default: 20 gC/m2)
#'
<<<<<<< HEAD
=======
#' @param avg_nyear_args list of arguments to be passed to
#' \link[boundaries]{average_nyear_window} (see for more info).
#' To be used for time series analysis
#'
#' @param biocol_option which biocol values to use for aggregation. options:
#' netsum, only_above_zero, abs - TODO: finish
#'
>>>>>>> development
#' @param ... arguments forwarded to \link[boundaries](classify_biomes)
#'
#' @return pb_status list data object
#'
#' @examples
#' \dontrun{
#' }
#' @export
calc_biosphere_status <- function(
  files_scenario,
  files_reference,
  time_span_scenario = as.character(1982:2011),
  time_span_reference = time_span_scenario,
  spatial_scale = "subglobal",
  thresholds = NULL,
  method = "stenzel2023",
  time_aggregation_args = list(),
  config_args = list(),
  path_baseline,
  time_span_baseline = time_span_reference,
  gridbased = TRUE,
  npp_threshold = 20,
<<<<<<< HEAD
=======
  avg_nyear_args = list(),
  biocol_option = "only_above_zero",
>>>>>>> development
  ...
) {

  files_baseline <- lapply(
    files_scenario,
    function(x) {
      file.path(path_baseline, basename(x))
    }
  )

  # workaround to filter input files (files_scenario and files_reference match)
  files_baseline <- mapply(
    function(x, y, z) {
      if (x == y) {
        z <- y
      } else {
        z
      }
      z
    },
    files_scenario,
    files_reference,
    files_baseline,
    SIMPLIFY = FALSE
  )

  biocol <- biospheremetrics::read_calc_biocol(
    files_scenario = files_scenario,
    files_baseline = files_baseline,
    files_reference = files_reference,
    time_span_scenario = time_span_scenario,
    time_span_baseline = time_span_baseline,
    time_span_reference = time_span_reference,
    gridbased = gridbased,
    npp_threshold = npp_threshold,
    # further options that can be defined in biospheremetrics:: read_calc_biocol
    # are diregarded here:
    read_saved_data = FALSE,
    save_data = FALSE,
    data_file = NULL,
    include_fire = FALSE,
    external_fire = FALSE,
    external_wood_harvest = FALSE,
    grass_scaling = FALSE,
    grass_harvest_file = NULL,
    external_fire_file = NULL,
    external_wood_harvest_file = NULL
  ) %>% suppressMessages()

  if (spatial_scale == "grid") {
    control_variable_raw <- abs(biocol$biocol_frac_piref)

  } else if (spatial_scale == "subglobal") {

    # classify biomes
    # Filter out method and thresholds arguments from ellipsis (thresholds
    # is not used in classify_biomes and method should not be passed on to use
    # the default method defined for classify biomes)
    ellipsis_filtered <- list(...)
    ellipsis_filtered$method <- NULL
    ellipsis_filtered$thresholds <- NULL

    # classify biomes based on foliage projected cover (FPC) output
    biome_classes <- do.call(
      classify_biomes,
      append(list(files_reference = files_baseline,
               time_span_reference = time_span_baseline,
               time_aggregation_args = time_aggregation_args
             ),
             ellipsis_filtered)
    )


    # TODO: also for global
    # initialize control variable vector
    if (biocol_option == "abs"){
      control_variable_raw <- abs(biocol$biocol)
    }else if (biocol_option == "only_above_zero"){
      control_variable_raw <- biocol$biocol
      control_variable_raw[control_variable_raw<0] <- 0
    }else if (biocol_option == "netsum"){
      control_variable_raw <- biocol$biocol
    }else { # TODO: do with matcharg
      stop("Not defined option for biocol_option.")
    }

    # get continents mask - pass arg of whether to merge europe and asia
    continent_grid %<-% calc_continents_mask(files_reference$grid) # implicit eurasia = TRUE

<<<<<<< HEAD
    terr_area <- lpjmlkit::read_io(
      files_scenario$terr_area
    ) %>%
      conditional_subset(config_args$spatial_subset) %>%
      lpjmlkit::as_array()

    # TODO: this is averaged over all cells in a biomes, but
    # (i) could be weighted by area of each cell
    # (ii) the calculation should be done for each continent seperately (?)
    for (year_i in seq_len(dim(biome_classes$biome_id)["year"])) {
      biome_id <- biome_classes$biome_id[, year_i]
      for (b in sort(unique(biome_id))){
        biome_cells <- which(biome_id == b)
        control_variable_raw[biome_cells, year_i] <- (
          weighted.mean(
            abs(biocol$biocol)[biome_cells, year_i],
            w = terr_area[biome_cells]
          ) /
            weighted.mean(
              biocol$npp_ref[biome_cells, year_i],
              w = terr_area[biome_cells]
            )
        )
=======
    # create space of combinations to loop over (even though not all make sense)
    comb <- expand.grid(
      continent = sort(
        unique(lpjmlkit::asub(continent_grid, band = "continent"))
      ),
      biome = sort(unique(factor(biome_classes$biome_id)))
    )
    # get flexibly named time dimension
    third_dim <- names(dim(control_variable_raw))[
      !names(dim(control_variable_raw)) %in% c("cell")
    ]
    for (idx in seq_len(nrow(comb))) {
      subset_biocol <- control_variable_raw
      subset_npp <- biocol$npp_ref

      # replace for every combination a subset of cells with mean of each
      sub_cells <- {
        # match biome type for cells
        array(
          drop(biome_classes$biome_id) == comb$biome[idx],
          dim = dim(control_variable_raw),
          dimnames = dimnames(control_variable_raw)
        ) &
          # match continent for cells
          array(
            lpjmlkit::asub(continent_grid, band = "continent", drop = FALSE) ==
                                            comb$continent[idx],
            dim = dim(control_variable_raw),
            dimnames = dimnames(control_variable_raw)
          )
      }
      if (!any(sub_cells)) {
        next
>>>>>>> development
      }
      subset_biocol[!sub_cells] <- NA
      summed_biocol <- apply(
        subset_biocol,
        third_dim,
        function(x) {
          # mean over cell subset of forest type and continent
          # weighted by cell area
          return(rep(sum(x, na.rm = TRUE), length(x)))
        }
      )
      subset_npp[!sub_cells] <- NA
      summed_npp <- apply(
        subset_npp,
        third_dim,
        function(x) {
          # mean over cell subset of forest type and continent
          # weighted by cell area
          return(rep(sum(x, na.rm = TRUE), length(x)))
        }
      )
      control_variable_raw[sub_cells] <- summed_biocol[sub_cells] /
                                         summed_npp[sub_cells]
    }

  } else if (spatial_scale == "global") {
    # TODO: also for global
    # initialize control variable vector
    if (biocol_option == "abs"){
      control_variable_raw <- biocol$biocol_overtime_abs_frac_piref
    }else if (biocol_option == "only_above_zero"){
      stop("Missing")
      # TODO: recompute, does not exist as overtime yet
    }else if (biocol_option == "netsum"){
      control_variable_raw <- biocol$biocol_overtime_frac
    }else { # TODO: do with matcharg
      stop("Not defined option for biocol_option.")
    }
  }

  # average
  control_variable <- do.call(aggregate_time,
                              append(list(x = control_variable_raw),
                                     time_aggregation_args))

  attr(control_variable, "spatial scale") <- spatial_scale
  attr(control_variable, "thresholds") <- thresholds
  attr(control_variable, "control_variable") <- "BioCol (in fraction of NPPref)"

  class(control_variable) <- c("control_variable")
  return(control_variable)

} # end of calc_biosphere_status
