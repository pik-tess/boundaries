#' Calculate biosphere status based on BioCol (HANPP) from a PNV run (reference) and
#' LU run (scenario) of LPJmL, both using time_span_scenario. Additionally
#' a separate reference NPP file (e.g. from a Holocene run) can be supplied as
#' reference_npp_file, which will use time_span_reference, or file index years
#' 3:32 if time_span_reference is not supplied
#'
#' Calculate biosphere status based on BioCol from a PNV run and LU run of LPJmL
#' @param files_scenario list with variable names and corresponding file paths
#' (character string) of the scenario LPJmL run. All needed files need to be
#' provided. E.g.: list(grid = "/temp/grid.bin.json",
#'                      npp = "/temp/npp.bin.json")
#' @param files_baseline list with variable names and corresponding file paths
#' (character string) of the baseline LPJmL run. All needed files are
#' provided in XXX. E.g.: list(npp = "/temp/npp.bin.json"). If not
#' needed for the applied method, set to NULL.
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the reference NPP, HANPP should be compared against.
#' In this case only NPP is required. list(npp = "/temp/npp.bin.json").
#' @param time_span_scenario time span to be used for the scenario run and
#' parallel PNV run, defined as a character string,
#' e.g. `as.character(1982:2011)`
#' @param time_span_baseline time span to be used for the baseline run, defined
#' as a character vector, e.g. `as.character(1901:1930)`. Can differ in offset
#' and length from `time_span_scenario`! If `NULL` value of `time_span_scenario`
#' is used
#' @param time_span_reference time span to read reference_npp_file from, using
#' index years 3:32 if set to NULL (default: NULL)
#' @param spatial_resolution character string indicating spatial resolution
#' either "grid" or "global"
#' @param thresholds named character string with thresholds to be used to
#' define the lower end of safe, increasing risk and high risk zone,
#' e.g. c(holocene = 0.0, pb = 0.1, highrisk = 0.2).
#' @param gridbased logical are pft outputs gridbased or pft-based?
#' @param npp_threshold lower threshold for npp (to mask out non-lu areas
#' according to Haberl et al. 2007). Below BioCol will be set to 0.
#' (default: 20 gC/m2)
#' @param avg_nyear_args list of arguments to be passed to
#'        \link[pbfunctions]{average_nyear_window} (see for more info).
#'        To be used for time series analysis
#' @param ... arguments forwarded to \link[boundaries](classify_biomes)
#'
#' @return pb_status list data object
#'
#' @examples
#' \dontrun{
#' }
#' @export
calc_biosphere_status <- function(files_scenario,
                                  files_baseline,
                                  files_reference = NULL,
                                  time_span_scenario,
                                  time_span_baseline = NULL,
                                  time_span_reference = NULL,
                                  spatial_resolution = "subglobal",
                                  thresholds = NULL,
                                  method = "stenzel2023",
                                  gridbased = T,
                                  npp_threshold = 20,
                                  avg_nyear_args = list(),
                                  ...
                                  ) {
  if (is.null(files_reference))
    files_reference <- list(npp = baseline_npp_file)
  if (is.null(time_span_baseline))
    time_span_baseline <- time_span_scenario
  if (is.null(time_span_reference))
    time_span_reference <- time_span_scenario[3:12]
  if (is.null(thresholds)) {
    thresholds <- c(holocene = 0,
                    pb = 0.1,
                    highrisk = 0.2)
  }
  biocol <- biospheremetrics::read_calc_biocol(
                    files_scenario = files_scenario,
                    files_baseline = files_baseline,
                    files_reference = files_reference,
                    time_span_scenario = as.character(time_span_scenario),
                    time_span_baseline = as.character(time_span_baseline),
                    time_span_reference = as.character(time_span_reference),
                    gridbased = gridbased,
                    npp_threshold = npp_threshold,
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

  if (spatial_resolution == "grid") {
    control_variable_raw <- abs(biocol$biocol_frac_piref)
  } else if (spatial_resolution == "subglobal") {
    # Filter out method and thresholds arguments from ellipsis
    ellipsis_filtered <- list(...)
    ellipsis_filtered$method <- NULL
    ellipsis_filtered$thresholds <- NULL
    # classify biomes based on foliage projected cover (FPC) output
    biome_classes <- do.call(classify_biomes,
                             append(list(files_reference = files_baseline,
                                         time_span_reference = time_span_baseline,
                                         avg_nyear_args = avg_nyear_args,
                                         montane_arctic_proxy = NULL # todo: ideally this line should be removed, but has to be fixed first - does not work yet with elevation
                                         ),
                                    ellipsis_filtered))

    # initialize control variable vector
    control_variable_raw <- biocol$biocol*0

    for (b in sort(unique(biome_classes$biome_id))){
      biome_cells <- which(biome_classes$biome_id == b)
      if (length(biome_cells) > 1){
        control_variable_raw[biome_cells,] <- colSums(abs(biocol$biocol)[biome_cells,])/
          sum(rowMeans(biocol$npp_ref[biome_cells,]))
      } else if (length(biome_cells) == 1){
        control_variable_raw[biome_cells,] <- abs(biocol$biocol)[biome_cells,]/
          mean(biocol$npp_ref[biome_cells,])
      }
    }
  } else if (spatial_resolution == "global") {
    # add dummy cell dimension
    control_variable_raw <- array(biocol$biocol_overtime_abs_frac_piref,
                               dim = c(cell = 1,
                                       year = length(biocol$biocol_overtime_abs_frac_piref)
                               ),
                               dimnames = list(cell = 1,
                                               year = names(biocol$biocol_overtime_abs_frac_piref)
                               ))
  }else{
    stop(paste("Unknown value for spatial_resolution: ", spatial_resolution))
  }
  # average
  control_variable <- do.call(average_nyear_window,
                        append(list(x = control_variable_raw),
                               avg_nyear_args))
  if (spatial_resolution == "global") {
    # remove dummy cell dimension
    control_variable <- control_variable[1,]
  }
  attr(control_variable, "thresholds") <- thresholds
  return(control_variable)

} # end of calc_biosphere_status
