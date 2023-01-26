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
                        path_reference = NULL,
                        time_span_scenario = c(1982, 2011),
                        time_span_reference = NULL,
                        avg_nyear_args = list(),
                        input_files = list(),
                        diff_output_files = list(),
                        read_args = list(),
                        in_parallel = TRUE,
                        # args that use all get prefix like lsc.method or
                        #   lsc.threshold
                        ...) {

  file_type <- get_file_type(path_scenario)

  # if in_parallel use future package for asynchronous parallelization
  if (in_parallel) {
    if (.Platform$OS.type == "windows") {
      future_plan <- future::plan("multisession")
    } else {
      future_plan <- future::plan("multicore")
    }
    on.exit(future::plan(future_plan))
  }

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


get_file_type <- function(path) {
  # get all files in path
  all_files <- list.files(
    path,
    full.names = TRUE
  )

  # get file extensions
  all_file_types <- all_files %>%
  strsplit("^([^\\.]+)") %>%
    sapply(function(x) {
      y <- x[2]
      return(y)
    }) %>%
    substr(2, nchar(.))

  # get most frequent file types
  most_frequent <- all_file_types %>%
    factor() %>%
    table() %>%
    names() %>%
    .[1:5]

  # 5 exemplaric files to detect type
  files_to_check <- sapply(
    most_frequent,
    function(x, y, z) {
      y[which(z == x)[1]]
    },
    y = all_files,
    z = all_file_types)

  # detect actual LPJmL data type
  types <- sapply(
    files_to_check,
    lpjmlkit:::detect_type
  )

  # assign file type after ranking which is available
  #   first preferable: "meta", second: "clm", last: "raw"
  if ("meta" %in% types) {
   file_type <- "meta"
  } else if ("clm" %in% types){
   file_type <- "clm"
  } else if ("raw" %in% types){
   file_type <- "raw"
  }
  return(file_type)
}