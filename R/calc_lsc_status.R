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
#' @param method method (character string) to be used , currently available
#' method is `c("gerten2020")` based on
#' [Gerten et al. 2020](https://doi.org/10.1038/s41893-019-0465-1).
#'
#' @param temporal_resolution character. Temporal resolution, available options
#' are `"annual"` (default) and `"monthly"`
#'
#' @param cut_min double. Exclude boundary calculations for Q < cut_min
#'
#' @param prefix_monthly_output character. Provide a prefix if required for
#' monthly LPJmL output files, e.g. `"m"` for `"mdischarge.bin"` instead of
#' `"discharge.bin"` (default is `""`)
#'
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#' time series analysis
#'
#' @examples
#' \dontrun{
#'  calc_bluewater_status(path_scenario, path_reference)
#' }
#'
#' @md
#' @export
calc_lsc_status <- function(path_scenario,
                            path_reference,
                            time_span_scenario = c(1982, 2011),
                            time_span_reference = NULL,
                            method = "gerten2020",
                            temporal_resolution = "annual",
                            # ???
                            cut_min = 0.0864,
                            prefix_monthly_output = "",
                            avg_nyear_args = list(),
                            # to be replaced by lpjmlKit::read_output
                            start_year = 1901) {
  
}