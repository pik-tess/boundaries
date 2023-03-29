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
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for # nolint
#' time series analysis
#'
#' @param input_files list of required file(s) using ID (e.g. `prec`, `runoff`)
#' and an absolute file path to the corresponding input file
#'
#' @param diff_output_files list of required file(s) using ID 
#' (e.g. prec, runoff) and the alternative writing (e.g. `"my_runoff"`)
#'
#' @param in_parallel if parallel (asynchronous) execution of code should be
#' used(future package) to speed up boundary status calculation
#'
#' @param ... further arguments to be passed to each calc_* function
#'
#' @return list with data array for each `boundary`
#'
#' @examples
#'
#' \dontrun{
#'  boundary_status <- calc_status(
#'    boundary = "nitrogen"
#'    path_scenario = "./my_scenario/output")
#' }
#'
#' @md
#' @export
calc_status <- function(boundary,
                        path_scenario,
                        path_reference = NULL,
                        time_span_scenario = as.character(1982:2011),
                        time_span_reference = time_span_scenario,
                        avg_nyear_args = list(),
                        input_files = list(),
                        diff_output_files = list(),
                        in_parallel = TRUE,
                        ...) {

  # If in_parallel use future package for asynchronous parallelization
  if (in_parallel) {
    if (.Platform$OS.type == "windows") {
      future_plan <- future::plan("multisession")
    } else {
      future_plan <- future::plan("multicore")
    }
    on.exit(future::plan(future_plan))
  }

  # Get main file type (meta, clm)
  file_ext <- get_file_ext(path_scenario)

  # List required output files for each boundary
  output_files <- list_outputs(boundary,
                               only_first_filename = FALSE)

  # Get filenames for scenario and reference
  files_scenario <- get_filenames(
    path = path_scenario,
    output_files = output_files,
    diff_output_files = diff_output_files,
    input_files = input_files,
    file_ext = file_ext
  )
  files_reference <- get_filenames(
    path = path_scenario,
    output_files = output_files,
    diff_output_files = diff_output_files,
    input_files = input_files,
    file_ext = file_ext
  )

  # Get arguments for each boundary function
  fun_args <- list_function_args(boundary)

  dot_args <- check_args <- list(...)
  all_status <- list()

  # Loop over boundaries and calculate status
  for (bound in boundary) {
    fun_name <- paste0("calc_", bound, "_status")

    # Get arguments for each boundary function
    sub_dots <- get_dots(fun_name, fun_args, dot_args)
    check_args[names(sub_dots)] <- NULL

    # Calculate status
    all_status[[bound]] <- do.call(
      fun_name,
      args = c(
        list(
          files_scenario = files_scenario,
          files_reference = files_reference,
          time_span_scenario = time_span_scenario,
          time_span_reference = time_span_reference
        ),
        sub_dots
      )
    )
    # Check if all arguments were used
    if (length(check_args) != 0) {
      warning(paste0("The following arguments were not used: ",
                     paste0("`", names(check_args), "`", collapse = ", ")))
    }
  }

  return(all_status)
}


get_file_ext <- function(path) {
  # Get all files in path
  all_files <- list.files(
    path,
    full.names = TRUE
  )

  # Get file extensions
  all_file_types <- all_files %>%
  strsplit("^([^\\.]+)") %>%
    sapply(function(x) {
      y <- x[2]
      return(y)
    }) %>%
    substr(2, nchar(.))

  # Get most frequent file types
  #TODO not yet working
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

  # Detect actual LPJmL data type
  types <- sapply(
    files_to_check,
    lpjmlkit:::detect_io_type
  ) %>%
  setNames(names(.), .)

  # Assign file type after ranking which is available
  # first preferable: "meta", second: "clm", last: "raw"
  if ("meta" %in% names(types)) {
   file_type <- types["meta"]
  } else if ("clm" %in% names(types)) {
   file_type <- types["clm"]
  } else if ("raw" %in% names(types)) {
   file_type <- types["raw"]
  }
  return(file_type)
}


get_filenames <- function(path,
                          output_files,
                          diff_output_files,
                          input_files,
                          file_ext) {

  file_names <- list()
  # Iterate over required outputs
  for (ofile in names(output_files)) {

  # Get required max. temporal resolution and convert to nstep
    resolution <- output_files[[ofile]]$resolution
    nstep <- switch(
      resolution,
      annual = 1,
      monthly = 12,
      daily = 365,
      stop(paste0("Not supported time resolution: ", dQuote(nstep), "."))
    )

    # If input file supplied use it as first priority
    if (ofile %in% names(input_files)) {
      file_name <- input_files[[ofile]]

    } else if (ofile %in% names(diff_output_files)) {

      # If different output file should be used - as second priority
      file_name <- paste0(
        path, "/",
        diff_output_files[[ofile]], ".",
        file_ext
      )
    } else {
      file_name <- NULL
    }

    if (!is.null(file_name)) {

      # Check if data could be read in
      meta <- lpjmlkit::read_meta(file_name)

      # Then check if temporal resultion of file matches required nstep
      if (nstep != meta$nstep && nstep != meta$nbands) {
        stop(
          paste0(
            "Required temporal resolution (nstep = ", nstep, ") ",
            "not supported by file ", dQuote(file_name),
            " (", meta$nstep, ")"
          )
        )
      }

    # If nothing specified try to read required files from provided path
    } else {

      # Iterate over different used file name options (e.g. runoff, mrunoff, ...) # nolint
      for (cfile in seq_along(output_files[[ofile]]$file_name)) {
        file_name <- paste0(
          path, "/",
          output_files[[ofile]]$file_name[cfile], ".",
          file_ext
        )

        # Check if file exists and if so check required temporal resolution
        # else next
        if (file.exists(file_name)) {
          meta <- lpjmlkit::read_meta(file_name)
          if (nstep <= meta$nstep || nstep == meta$nbands) {
            # Matching file found, break and use current file_name
            break
          }
        }

        # At end of iteraton raise error that no matching file_name was found
        if (cfile == length(output_files[[ofile]]$file_name) &&
            !output_files[[ofile]]$optional) {
          stop(
            paste0(
              "No matching output for ", dQuote(ofile),
              " with required temporal resolution (nstep = ", nstep, ") ",
              "found at path ", dQuote(path), "."
            )
          )
        }
      }
    }
    file_names[[ofile]] <- file_name
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
