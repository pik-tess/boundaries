#' Status calculation of the greenwater boundary
#'
#' Planetary Boundary status calculation of the greenwater boundary based on
#' rootmoisture in a scenario LPJmL run and a reference LPJmL run.
#'
#' @param files_scenario list with variable names and corresponding file paths
#' (character string) of the scenario LPJmL run. Handled automatically via
#' [`calc_status()`].
#'
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the files_reference LPJmL run. Handled automatically via
#' [`calc_status()`].
#'
#' @param spatial_scale character string indicating spatial resolution
#' either "grid", "subglobal" or "global" for calculation of the share (%)
#' of total global area with deviations
#'
#' @param time_span_scenario time span to use output from the scenario run,
#' e.g. `1982:2011`.
#'
#' @param time_span_reference time span use output from the reference run,
#' e.g. `1901:1930`.
#'
#' @param approach approach (character string) to be used , currently available
#' approach is `c("wang-erlandsson2022")` based on
#' [Wang-Erlandsson et al. 2022](https://doi.org/10.1038/s43017-022-00287-8)
#' (referring only to the driest/wettest month of each year) or
#' `porkka2024` based on
#' [Porkka et al. 2023](https://eartharxiv.org/repository/view/3438/)
#' (referring to each month of a year)
#'
#' @param time_series_avg integer. Number of years to be used for the moving
#' average calculation. If `NULL`, all years are averaged for one status
#' calculation, for `1` the whole time span is used to calculate a status time
#' series.
#'
#' @param config_args list of arguments to be passed on from the model
#' configuration.

#' @param thresholds named character string with thresholds to be used to
#' define the safe, increasing risk and high risk zone,
#' e.g. c(holocene = 0.5, pb = 0.95, highrisk = 0.99).
#' For spatial resolution = "grid", this refers to the p value
#' (significance level of increases in deviations) with the default:
#' c(holocene = 1, pb = 0.05, highrisk = 0.01).
#' For spatial resolution = "global", this refers to the quantiles of
#' the global area with deviations in the reference period. The dafault
#' for global resolution is: c(holocene = 0.5, pb = 0.95,
#' highrisk = 0.99).
#' If set to NULL, the respective default is taken (see above; matching
#' the spatial_scale, defined in metric_files.yml).
#'
#'@return Object of class `control_variable` with the boundary status of the
#' greenwater boundary.
#'
#' @examples
#' \dontrun{
#' boundary_status <- calc_status(
#'   boundary = "greenwater",
#'   config_scenario = "path/to/config_scenario.json",
#'   config_reference = "path/to/config_reference.json",
#'   spatial_scale = "global",
#'   time_span_scenario = 1901:2019,
#'   time_span_reference = 1901:1930,
#'   approach = "porkka2024"
#' )
#' }
#'
#' @md
#' @export
greenwater_status <- function(
  files_scenario,
  files_reference,
  spatial_scale = "global",
  time_span_scenario = as.character(1982:2011),
  time_span_reference = time_span_scenario,
  approach = "wang-erlandsson2022",
  time_series_avg = NULL,
  config_args = list(),
  thresholds = NULL
) {

  # verify available methods
  approach <- match.arg(approach, c("wang-erlandsson2022", "porkka2024"))

  # verify available spatial resolution
  spatial_scale <- match.arg(spatial_scale, c("global", "subglobal"))

  # -------------------------------------------------------------------------- #
  # calc deviations for rootmoisture
  control_variable <- calc_water_deviations(
    files_scenario = files_scenario,
    files_reference = files_reference,
    spatial_scale = spatial_scale,
    time_span_scenario = time_span_scenario,
    time_span_reference =  time_span_reference,
    approach = approach,
    time_series_avg = time_series_avg,
    config_args = config_args,
    thresholds = thresholds,
    variable = "rootmoist"
  )

  control_variable <- set_attributes(
    control_variable,
    approach,
    "greenwater",
    spatial_scale,
    thresholds
  )
  return(control_variable)
}
