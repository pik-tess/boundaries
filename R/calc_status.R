#' Calculate the planetary boundary status
#'
#' Calculate the PB status for a defined planetary boundary based
#' on a scenario LPJmL run and a reference LPJmL run.
#'
#' @param boundary character vector, boundary for which status is calculated. 
#' Available terrestrial boundaries are c("bluewater", "greenwater", "lsc",
#' "nitrogen")
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
#' @param ... further arguments to be passed to each calc_* function
#'
#' @return list with data array for each `boundary`
#'
#' @examples
#' \dontrun{ 
#'  boundary_status <- calc_status(
#'    boundary = 
#'    path_scenario = "./my_scenario/output",
#'    path_reference = "./my_reference/output")
#' }
#'
#' @md
#' @export
calc_status <- function(boundary,
                        path_scenario,
                        path_reference,
                        time_span_scenario = c(1982, 2011),
                        time_span_reference = NULL,
                        ...) {

  all_status <- list()
  # utility functions
  for (bound in boundary) {
    all_status[[boundary]] <- do.call(
      paste0("calc_", boundary, "_status"),
      args = list(path_scenario = path_scenario,
                  path_reference = path_reference,
                  time_span_scenario = time_span_scenario,
                  time_span_reference = time_span_reference)
    )
  }

  return(all_status)
}