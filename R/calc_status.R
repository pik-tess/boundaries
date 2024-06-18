#' Calculate the planetary boundary status
#'
#' Calculate the PB status for a defined planetary boundary based
#' on a scenario LPJmL run and a reference LPJmL run. For boundary function
#' specific arguments to be passed (via `...`) see the respective function
#' documentation of [`biosphere_status()`], [`nitrogen_status()`],
#' [`greenwater_status()`], [`bluewater_status()`] or [`lsc_status()`]
#'
#' @param boundary character vector, boundary for which status is calculated.
#' Available terrestrial boundaries are `c("bluewater", "greenwater", "lsc",
#' "nitrogen", "biosphere")`.
#'
#' @param config_scenario character string. File path to the LPjmL configuration
#' file (json) of the scenario run. The configuration file contains the
#' information about the LPJmL run, e.g. the output directory
#'
#' @param config_reference character string. See config_scenario. For the
#' reference run
#'
#' @param spatial_scale character string indicating spatial resolution
#' options: "global", "subglobal", "grid";
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#' as an integer (or character) vector, e.g. `1982:2011` (default)
#'
#' @param time_span_reference time span to be used for the scenario run, defined
#' as an integer (or character) vector, e.g. `1901:1930`. Can differ in offset
#' and length from `time_span_scenario`! If `NULL` value of `time_span_scenario`
#' is used
#'
#' @param time_series_avg integer. Number of years to be used for the moving
#' average calculation. If `NULL`, all years are averaged for one status
#' calculation, for `1` the whole time span is used to calculate a status time
#' series.
#'
#' @param approach list of methods to be used for each boundary. If `NULL` the
#' default approach is used
#'
#' @param thresholds list of thresholds to be used for each boundary. If `NULL`
#' the default thresholds are used
#'
#' @param in_parallel logical, if `TRUE` the function uses parallelization
#' (default) based on the `future` package (asynchronous execution). If `FALSE`
#' no parallelization is used
#'
#' @param ... further arguments to be passed to each calc_* function
#'
#' @return list with objects of class `control_variable`. To directly get the
#' `boundary_status` use [`as_risk_level()`].
#'
#' @examples
#'
#' \dontrun{
#' boundary_status <- calc_status(
#'   boundary = c("biosphere","nitrogen", "greenwater", "bluewater", "lsc")
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
calc_status <- function(
  boundary,
  config_scenario,
  config_reference,
  spatial_scale,
  time_span_scenario = 1982:2011,
  time_span_reference = time_span_scenario,
  time_series_avg = NULL,
  approach = list(),
  thresholds = list(),
  in_parallel = TRUE,
  ...
) {

  # If in_parallel use future package for asynchronous parallelization
  if (in_parallel) {
    rlang::local_options(
      future.globals.maxSize = 8000 * 1024^2,
      future.rng.onMisuse = "ignore"
    )
    if (.Platform$OS.type == "windows") {
      future_plan <- future::plan("multisession")
    } else {
      future_plan <- future::plan("multicore")
    }
    on.exit(future::plan(future_plan)) # nolint:undesirable_function_linter
  }

  # verify available spatial resolution
  spatial_scale <- match.arg(
    spatial_scale,
    c("global", "subglobal", "grid")
  )

  config_scenario <- lpjmlkit::read_config(config_scenario)
  config_reference <- lpjmlkit::read_config(config_reference)

  # check time span and convert to character if necessary
  time_span_scenario <- check_time_span(
    time_span = time_span_scenario,
    config = config_scenario,
    type = "scenario"
  )
  time_span_reference <- check_time_span(
    time_span = time_span_reference,
    config = config_reference,
    type = "reference"
  )

  # List required output files for each boundary
  output_files <- list_outputs(
    metric = boundary,
    spatial_scale = spatial_scale,
    approach = approach,
    only_first_filename = FALSE
  )

  # Get filenames for scenario and reference
  files_scenario <- get_filenames(
    config = config_scenario,
    output_files = output_files
  )
  files_reference <- get_filenames(
    config = config_reference,
    output_files = output_files
  )

  config_args <- list()
  config_args$spatial_subset <- get_spatial_subset(
    config_scenario,
    config_reference
  )

  config_args$gridbased <- config_scenario$grid_scaled

  # Get arguments for each boundary function
  fun_args <- list_function_args(boundary)

  dot_args <- check_args <- list(...)
  all_status <- list()

  # Loop over boundaries and calculate status
  for (bound in boundary) {
    fun_name <- paste0(bound, "_status")

    # Get arguments for each boundary function
    sub_dots <- get_dots(fun_name, fun_args, dot_args)
    check_args[names(sub_dots)] <- NULL

    inner_args <- list(
      files_scenario = files_scenario,
      files_reference = files_reference,
      time_span_scenario = time_span_scenario,
      time_span_reference = time_span_reference,
      spatial_scale = spatial_scale,
      time_series_avg = time_series_avg,
      config_args = config_args
    )

    if (length(approach[[bound]]) > 0) {
      inner_args$approach <- method_i <- approach[[bound]]
    } else {
      method_i <- formals(get(fun_name))$approach
    }

    if (length(thresholds[[bound]]) > 0) {
      inner_args$thresholds <- thresholds[[bound]]
    } else {
      inner_args$thresholds <- list_thresholds(bound, method_i, spatial_scale)
    }

    # Calculate status
    message(bound)
    all_status[[bound]] <- do.call(
      fun_name,
      args = c(
        inner_args,
        sub_dots
      )
    )
  }

  # Check if all arguments were used
  if (length(check_args) != 0) {
    warning(paste0("The following arguments were not used: ",
                   paste0("`", names(check_args), "`", collapse = ", ")))
  }

  return(all_status)
}


get_filenames <- function(config,
                          output_files) {

  file_names <- list()
  sim_data <- get_sim_data(config)
  # Iterate over required outputs
  for (data_file in names(output_files)) {

    # Get required max. temporal resolution and convert to nstep
    resolution <- output_files[[data_file]]$resolution
    nstep <- switch(
      resolution,
      annual = 1,
      monthly = 12,
      daily = 365,
      stop(paste0("Not supported time resolution: ", dQuote(nstep), "."))
    )

    if (data_file %in% names(sim_data)) {
      file_name <- sim_data[[data_file]]
    } else {
      message("Output file ", dQuote(data_file), " not found in config.")
      file_name <- NULL
    }

    if (!is.null(file_name) && endsWith(file_name, "json")) {

      # Check if data could be read in
      meta <- lpjmlkit::read_meta(file_name)

      # Then check if temporal resolution of file matches required nstep
      if (nstep > meta$nstep && nstep > meta$nbands) {
        stop(
          paste0(
            "Required temporal resolution (nstep = ", nstep, ") ",
            "not supported by file ", dQuote(file_name),
            " (", meta$nstep, ")"
          )
        )
      }
    }
    file_names[[data_file]] <- file_name
  }
  file_names
}

get_dots <- function(fun_name, fun_args, dot_args) {
  sub_dots <- list()
  for (dot_i in seq_along(dot_args)) {
    name_arg <- names(dot_args[dot_i])
    if (name_arg %in% fun_args[[fun_name]]) {
      sub_dots <- c(sub_dots, dot_args[dot_i])
    }
  }
  sub_dots
}


# Avoid note for "."...
utils::globalVariables(".") # nolint:undesirable_function_linter

# Utility function to check time span matches simulation config and convert to
#   character if necessary
check_time_span <- function(time_span, config = NULL, type = "scenario") {

  if (!is.null(config) && !all(time_span %in% get_sim_time(config))) {
    stop("Time span not available in ", type, " run.")
  }
  if (is.numeric(time_span)) {
    time_span <- as.character(time_span)
  }
  time_span
}
