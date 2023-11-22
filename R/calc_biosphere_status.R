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
#' @param path_baseline character string with path to outputs for the baseline
#' run, file names are taken from files scenario. 
#'
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the reference NPP, HANPP should be compared against.
#' In this case only NPP is required. list(npp = "/temp/npp.bin.json").
#'
#' @param time_span_scenario time span to be used for the scenario run and
#' parallel PNV run, defined as a character string,
#' e.g. `as.character(1982:2011)`
#'
#' @param time_span_baseline time span to be used for the baseline run, defined
#' as a character vector, e.g. `as.character(1901:1930)`. Can differ in offset
#' and length from `time_span_scenario`! If `NULL` value of `time_span_scenario`
#' is used
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
#' @param gridbased logical; are pft outputs from LPJmL gridbased or pft-based?
#'
#' @param npp_threshold lower threshold for npp (to mask out non-lu areas
#' according to Haberl et al. 2007). Below BioCol will be set to 0.
#' (default: 20 gC/m2)
#'
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info).
#' To be used for time series analysis
#'
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
  path_baseline,
  time_span_scenario = as.character(1982:2011),
  time_span_baseline = time_span_scenario,
  time_span_reference = time_span_scenario,
  spatial_scale = "subglobal",
  thresholds = NULL,
  method = "stenzel2023",
  gridbased = TRUE,
  npp_threshold = 20,
  avg_nyear_args = list(),
  ...
) {

  files_baseline <- lapply(
    files_scenario,
    function(x) {
      file.path(path_baseline, basename(x))
    }
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
  )

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
               avg_nyear_args = avg_nyear_args
             ),
             ellipsis_filtered)
    )

    # initialize control variable vector
    control_variable_raw <- biocol$biocol * 0

    # TODO: this is averaged over all cells in a biomes, but
    # (i) could be weighted by area of each cell
    # (ii) the calculation should be done for each continent seperately (?)
    for (b in sort(unique(biome_classes$biome_id))){
      biome_cells <- which(biome_classes$biome_id == b)
      if (length(biome_cells) > 1) {
        control_variable_raw[biome_cells, ] <-
          colSums(abs(biocol$biocol)[biome_cells, ]) /
          sum(rowMeans(biocol$npp_ref[biome_cells, ]))
      } else if (length(biome_cells) == 1) {
        control_variable_raw[biome_cells, ] <-
          abs(biocol$biocol)[biome_cells, ] /
          mean(biocol$npp_ref[biome_cells, ])
      }
    }
  } else if (spatial_scale == "global") {
    # add dummy cell dimension
    # TODO can be removed once avg_nyear_args is compatible with it
    control_variable_raw <- array(biocol$biocol_overtime_abs_frac_piref,
                                  dim = c(cell = 1,
                                          year = length(biocol$biocol_overtime_abs_frac_piref) #nolint
                                  ),
                                  dimnames = list(cell = 1,
                                                  year = names(biocol$biocol_overtime_abs_frac_piref) #nolint
                                  ))
  } else {
    stop(paste("Unknown value for spatial_scale: ", spatial_scale))
  }
  # average
  control_variable <- do.call(average_nyear_window,
                              append(list(x = control_variable_raw),
                                     avg_nyear_args))

  if (spatial_scale == "global") {
    # remove dummy cell dimension
    control_variable <- control_variable[1, ]
  }
  attr(control_variable, "thresholds") <- thresholds
  attr(control_variable, "control variable") <- "BioCol (in fraction of NPPref)"
  return(control_variable)

} # end of calc_biosphere_status
