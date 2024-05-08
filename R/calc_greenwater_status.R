#' Calculate the planetary boundary status for the greenwater boundary
#'
#' Calculate the PB status for the greenwater (former freshwater) boundary based
#' on rootmoisture in a scenario LPJmL run and a reference LPJmL run.
#'
#' @param files_scenario list with variable names and corresponding file paths
#' (character string) of the scenario LPJmL run. All needed files are
#' provided in metric_files.yml
#'
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the reference LPJmL run. All needed files are
#' provided in metric_files.yml
#'
#' @param spatial_scale character string indicating spatial resolution
#' either "grid", "subglobal" or "global" for calculation of the share (%)
#' of total global area with deviations
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#' as a character string, e.g. `as.character(1982:2011)` (default)
#'
#' @param time_span_reference time span to be used for the reference run,
#' defined as a character string (e.g. `as.character(1901:1930)`).
#' Can differ in offset and length from `time_span_scenario`!
#' If `NULL` value of `time_span_scenario` is used
#'
#' @param time_aggregation_args list of arguments to be passed to
#' [`aggregate_time`] (see for more info).
#' To be used for time series analysis
#'
#' @param config_args list of arguments to be passed on from the model
#' configuration.
#'
#' @param approach approach (character string) to be used , currently available
#' approach is `c("wang-erlandsson2022")` based on
#' [Wang-Erlandsson et al. 2022](https://doi.org/10.1038/s43017-022-00287-8)
#' (referring only to the driest/wettest month of each year) or
#' `porkka2024` based on
#' [Porkka et al. 2023](https://eartharxiv.org/repository/view/3438/)
#' (referring to each month of a year)
#'
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
#' @examples
#' \dontrun{
#'  calc_greenwater_status(files_scenario, files_reference,
#'                 time_span_reference, spatial_scale = "grid")
#' }
#'
#' @md
#' @export
calc_greenwater_status <- function(
  files_scenario,
  files_reference,
  spatial_scale = "global",
  time_span_scenario = as.character(1982:2011),
  time_span_reference = time_span_scenario,
  approach = "wang-erlandsson2022",
  time_aggregation_args = list(),
  config_args = list(),
  thresholds = NULL
) {

  # verify available methods
  approach <- match.arg(approach, c("wang-erlandsson2022", "porkka2024"))

  # verify available spatial resolution
  spatial_scale <- match.arg(spatial_scale, c("global", "subglobal"))

  # -------------------------------------------------------------------------- #
  # calc deviations for rootmoisture
  deviations_rootmoisture <- calc_water_deviations(
    files_scenario = files_scenario,
    files_reference = files_reference,
    spatial_scale = spatial_scale,
    time_span_scenario = time_span_scenario,
    time_span_reference =  time_span_reference,
    approach = approach,
    time_aggregation_args = time_aggregation_args,
    config_args = config_args,
    thresholds = thresholds,
    variable = "rootmoist"
  )
  attr(deviations_rootmoisture, "long_name") <- list_long_name("greenwater")
  return(deviations_rootmoisture)
}
