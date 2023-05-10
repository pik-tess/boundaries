#' Calculate the planetary boundary status for the greenwater boundary
#'
#' Calculate the PB status for the greenwater (former freshwater) boundary based
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
#' @param method method (character string) to be used , currently available
#' method is `c("wang-erlandsson2022")` based on
#' [Wang-Erlandsson et al. 2022](https://doi.org/10.1038/s43017-022-00287-8).
#'
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#' time series analysis
#'
#' @examples
#' \dontrun{
#'  calc_greenwater_status(path_scenario, path_reference)
#' }
#'
#' @md
#' @export
calc_greenwater_status <- function(files_scenario,
                                   files_reference,
                                   time_span_scenario = as.character(1982:2011),
                                   time_span_reference = NULL,
                                   method = "wang-erlandsson2022",
                                   avg_nyear_args = list(),
                                   spatial_resolution
                                   ) {
   # verify available methods
  method <- match.arg(method, c("wang-erlandsson2022", "porkka_2023"))
   # verify available spatial resolution
  spatial_resolution <- match.arg(spatial_resolution, c("grid", "global"))
  # TODO not yet compatible with avg_nyear_args

  if (length(time_span_reference) < length(time_span_scenario)) {
    nyear_ref <- length(time_span_scenario)
  } else {
    nyear_ref <- NULL
  }

  # -------------------------------------------------------------------------- #
  # calc deviations for rootmoisture
  deviations_rootmoisture <- calc_deviations(
    file_scenario = files_scenario$rootmoist,
    file_reference = files_reference$rootmoist,
    grid_path = files_reference$grid,
    time_span_scenario = time_span_scenario,
    time_span_reference =  time_span_reference,
    method = method,
    avg_nyear_args = avg_nyear_args,
    spatial_resolution = spatial_resolution
  )

  if (spatial_resolution == "grid") {
    #prop.test to test for significance in departure increases
    # TODO make function that can also be used for bluewater
    p_dry <- p_wet <- array(NA, length(grid$data[, , 2]))
    if (method == "wang-erlandsson2022")
      trials <- c(length(time_span_reference), length(time_span_scenario))
    else if (method == "porkka_2023") {
      trials <- c(length(time_span_reference) * 12,
                   length(time_span_scenario) * 12)
    }
    for (i in seq_len(length(grid$data[, , 2]))) {
      test_dry <- prop.test(x = c(ref_n_depart$dry[i], scen_n_depart$dry[i]),
                              n = trials,
                              alternative = "less") #TODO verify!
      p_dry[i] <- test_dry$p.value

      test_wet <- prop.test(x = c(ref_n_depart$wet[i], scen_n_depart$wet[i]),
                              n = trials,
                              alternative = "less") #TODO verify!
      p_wet[i] <- test_wet$p.value
    }
    #TODO understand where NA values come from
    p_wet[is.na(p_wet)] <- 1
    p_dry[is.na(p_dry)] <- 1
    #TODO translation into PB status only prelimary. Thresholds should be
    # definable in function parameters
    pb_status <- array(0,
                       dim = dim(p_dry))
    pb_status[p_wet < 0.05 | p_dry < 0.05] <- 2
    pb_status[p_wet < 0.05 & p_dry < 0.05] <- 3
    pb_status[p_wet >= 0.05 & p_dry >= 0.05] <- 1
  } else if (spatial_resolution == "global") {
    # TODO translate into pb status
    pb_status <- deviations_rootmoisture$scenario
  }
  return(pb_status)
}