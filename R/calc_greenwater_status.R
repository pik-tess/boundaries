#' Calculate the planetary boundary status for the greenwater boundary
#'
#' Calculate the PB status for the greenwater (former freshwater) boundary based
#' on a scenario LPJmL run and a reference LPJmL run.
#'
#' @param files_scenario list with variable names and corresponding file paths
#'        (character string) of the scenario LPJmL run. All needed files are
#'        provided in metric_files.yml
#'
#' @param files_reference list with variable names and corresponding file paths
#'        (character string) of the reference LPJmL run. All needed files are
#'        provided in metric_files.yml
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#'        as a character string, e.g. `as.character(1982:2011)` (default)
#'
#' @param time_span_reference time span to be used for the reference run, defined
#'        defined as a character string (e.g. `as.character(1901:1930)`).
#'        Can differ in offset and length from `time_span_scenario`!
#'        If `NULL` value of `time_span_scenario` is used
#'
#' @param method method (character string) to be used , currently available
#'        method is `c("wang-erlandsson2022")` based on
#'        [Wang-Erlandsson et al. 2022](https://doi.org/10.1038/s43017-022-00287-8)
#'        (referring only to the driest/wettest month of each year) or
#'        `porkka_2023` based on
#'        [Porkka et al. 2023](https://eartharxiv.org/repository/view/3438/)
#'        (referring to each month of a year)
#' @param spatial_resolution character string indicating spatial resolution
#'        either "grid" for calculation of number of years with transgression
#'        (for wang-erlandsson2022: dim(ncell, nyears); 
#'         for porkka_2023: dim(ncell, nyears, months)) or
#'        "global" for calculation of the share (%) of total global area with
#'        deviations (either one value per year (wang-erlandsson2022) or one
#'        value per year and month (porkka_2023))
#'
#' @param avg_nyear_args list of arguments to be passed to
#'        \link[pbfunctions]{average_nyear_window} (see for more info).
#'        To be used for time series analysis
#'
#' @examples
#' \dontrun{
#'  calc_deviations(file_scenario, file_reference, grid_path,
#'                 time_span_reference, spatial_resolution)
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


  # check time_spans of scenario and reference runs
  if (is.null(time_span_reference)) {
    time_span_reference <- time_span_scenario
    nyear_ref <- NULL
  } else {
    if (length(time_span_reference) > length(time_span_scenario)) {
      stop(paste0("time_span_reference is longer than time_span_scenario.",
                  "Define a time_span_reference that is shorter than",
                  "time_span_scenario"))
    } else if (length(time_span_reference) < length(time_span_scenario)) {
      nyear_ref <- length(time_span_scenario)
    } else {
      nyear_ref <- NULL
    }
  }

  # -------------------------------------------------------------------------- #
  # calc deviations for rootmoisture
  deviations_rootmoisture <- calc_water_status(
    file_scenario = files_scenario$rootmoist,
    file_reference = files_reference$rootmoist,
    grid_path = files_reference$grid,
    time_span_scenario = time_span_scenario,
    time_span_reference =  time_span_reference,
    method = method,
    avg_nyear_args = avg_nyear_args,
    spatial_resolution = spatial_resolution
  )

  return(deviations_rootmoisture)
}