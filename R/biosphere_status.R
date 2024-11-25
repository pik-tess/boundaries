#' Status calculation of the biosphere integrity boundary.
#'
#' Biosphere status calculation based on BioCol (HANPP) from a baseline run
#' (with potential natural vegetation) and a scenario run (actual land use)
#' of LPJmL, both within the `time_span_scenario`.
#' Additionally a separate reference NPP file (e.g. from a Holocene run) can be
#' supplied with `files_reference` = list(npp = "path/to/npp.bin.json"),
#' which will use time_span_reference, or file index years 3:32 if
#' time_span_reference is not supplied.
#'
#' @param files_scenario list with variable names and corresponding file paths
#' (character string) of the scenario LPJmL run. All needed files need to be
#' provided. E.g.: list(grid = "/temp/grid.bin.json",
#' npp = "/temp/npp.bin.json"). Handled via [`calc_status`].
#'
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the reference NPP, HANPP should be compared against.
#' In this case only NPP is required. list(npp = "/temp/npp.bin.json").
#'
#' @param spatial_scale character string indicating spatial resolution
#' either "grid", "subglobal" or "global"
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#' as character string
#'
#' @param time_span_reference time span to be used for the reference run,
#' defined as character string, e.g. `as.character(1901:1930)`.
#'
#' @param approach approach (character string) to be used , currently available
#' approach is `"stenzel2023"`
#'
#' @param time_series_avg integer. Number of years to be used for the moving
#' average calculation. If `NULL`, all years are averaged for one status
#' calculation, for `1` the whole time span is used to calculate a status time
#' series.
#'
#' @param config_args list of arguments to be passed on from the model
#' configuration.
#'
#' @param thresholds named character string with thresholds to be used to
#' define the lower end of safe, increasing risk and high risk zone,
#' e.g. c(holocene = 0.0, pb = 0.1, highrisk = 0.2). If set to NULL, default
#' values from metric_files.yml will be used.
#'
#' @param path_baseline character string with path to outputs for the baseline
#' run, file names are taken from files scenario.
#'
#' @param time_span_baseline time span to be used for the baseline run, defined
#' as a character vector, e.g. `as.character(1901:1930)`. Can differ in offset
#' and length from `time_span_scenario`! If `NULL` value of `time_span_scenario`
#' is used
#'
#' @param npp_threshold lower threshold for npp (to mask out non-lu areas
#' according to Haberl et al. 2007). Below BioCol will be set to 0.
#' (default: 20 gC/m2)
#'
#' @param biocol_option which biocol values to use for aggregation. options:
#' netsum, only_above_zero, abs
#'
#' @param eurasia logical. If `spatial_scale` = `"subglobal"` merge continents
#' Europe and Asia to avoid arbitrary biome cut at europe/asia border.
#' Defaults to `TRUE`
#'
#' @param ... arguments forwarded to [`classify_biomes`]
#'
#'@return Object of class `control_variable` with the boundary status of the
#' biosphere integrity boundary.
#'
#' @examples
#' \dontrun{
#' boundary_status <- calc_status(
#'   boundary = "biosphere",
#'   config_scenario = "path/to/config_scenario.json",
#'   config_reference = "path/to/config_reference.json",
#'   spatial_scale = "global",
#'   time_span_scenario = 1901:2019,
#'   time_span_reference = 1901:1930,
#'   approach = "stenzel2023",
#'   path_baseline = "path/to/baseline_outputs"
#' )
#' }
#' @export
biosphere_status <- function(
  files_scenario,
  files_reference,
  spatial_scale = "subglobal",
  time_span_scenario = as.character(1982:2011),
  time_span_reference = NULL,
  approach = "stenzel2023",
  time_series_avg = NULL,
  config_args = list(),
  thresholds = NULL,
  path_baseline,
  time_span_baseline = time_span_scenario,
  npp_threshold = 20,
  biocol_option = "only_above_zero",
  eurasia = TRUE,
  ...
) {

  biocol_option <- match.arg(biocol_option,
                             c("only_above_zero", "netsum", "abs"))

  files_baseline <- lapply(
    files_scenario,
    function(x) {
      file.path(path_baseline, basename(x))
    }
  )

  # workaround to filter input files (files_scenario and files_reference match)
  files_baseline <- mapply( # nolint:undesirable_function_linter
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

  gridbased <- config_args$gridbased

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
    # are disregarded here:
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
    # Filter out approach and thresholds arguments from ellipsis (thresholds
    # is not used in classify_biomes and approach should not be passed on to use
    # the default approach defined for classify biomes)
    ellipsis_filtered <- list(...)
    ellipsis_filtered$approach <- NULL
    ellipsis_filtered$thresholds <- NULL

    # classify biomes based on foliage projected cover (FPC) output
    biome_classes <- do.call(
      classify_biomes,
      append(list(files_reference = files_baseline,
               time_span_reference = time_span_baseline,
               time_series_avg = NULL
             ),
             ellipsis_filtered)
    )


    # initialize control variable vector
    if (biocol_option == "abs") {
      control_variable_raw <- abs(biocol$biocol)
    }else if (biocol_option == "only_above_zero") {
      control_variable_raw <- biocol$biocol
      control_variable_raw[control_variable_raw < 0] <- 0
    }else if (biocol_option == "netsum") {
      control_variable_raw <- biocol$biocol
    }

    # please R CMD check for use of future operator
    continent_grid <- NULL
    # get continents mask - pass arg of whether to merge europe and asia
    continent_grid %<-% mask_continents(
      files_reference$grid,
      eurasia = eurasia
    )

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
      di_biocol <- dim(subset_biocol)
      subset_npp <- rep(rowMeans(biocol$npp_ref), times = di_biocol["year"])
      dim(subset_npp) <- dim(subset_biocol)
      dimnames(subset_npp) <- dimnames(subset_biocol)

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
            lpjmlkit::asub(
              continent_grid,
              band = "continent",
              drop = FALSE
            ) == comb$continent[idx],
            dim = dim(control_variable_raw),
            dimnames = dimnames(control_variable_raw)
          )
      }
      if (!any(sub_cells)) {
        next
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
    # initialize control variable vector
    if (biocol_option == "abs") {
      control_variable_raw <- biocol$biocol_overtime_abs_frac_piref
    }else if (biocol_option == "only_above_zero") {
      control_variable_raw <- biocol$biocol_overtime_pos_frac_piref
    }else if (biocol_option == "netsum") {
      control_variable_raw <- biocol$biocol_overtime_frac
    }
  } # end if spatial_scale == "global"

  rm(biocol)
  # average
  control_variable <- aggregate_time(
    x = control_variable_raw * 100,
    time_series_avg = time_series_avg
  )

  control_variable <- set_attributes(
    control_variable,
    approach,
    "biosphere",
    spatial_scale,
    thresholds
  )

  return(control_variable)

} # end of biosphere_status
