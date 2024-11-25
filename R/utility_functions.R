#' List required LPJmL outputs and temporal resolution
#'
#' Function to return a list of output IDs with required resolution and
#' file names for a given metric. The list is based on the metric_files.yml file
#' in the boundaries package (`"./inst/metric_files.yml"`).
#'
#' @param metric Character string containing name of metric to get
#'        required outputs. Available options are `c("biome",
#'        "nitrogen", "lsc", "bluewater", "greenwater", "biosphere")`
#'        or just `"all"` or `"benchmark"`. Default is `"all"`.
#'
#' @param approach List of character strings containing the approach to
#'       calculate the metric. Or `"all"` to get all approaches (default).
#'
#' @param spatial_scale character. Spatial resolution, available options
#'        are `"subglobal"` (at the biome level), `"global"` and
#'        `"grid"` or `"all"` (default).
#'
#' @param only_first_filename Logical. If TRUE, only the first file name will be
#'        returned for each output. If FALSE, all file names will be returned.
#'
#' @return List of output IDs with required resolution and
#'         file names for a given metric
#'
#' @examples
#' \dontrun{
#' list_outputs(
#'   "biome",
#'   approach = list("biome" = approach),
#'   spatial_scale = "subglobal",
#'   only_first_filename = FALSE
#' )
#' }
#' @export
list_outputs <- function(
  metric = "all",
  spatial_scale = "all",
  approach = "all",
  only_first_filename = TRUE
) {
  metric <- process_metric(metric = metric)

  system.file(
    "extdata",
    "metric_files.yml",
    package = "boundaries"
  ) %>%
    yaml::read_yaml() %>%
    get_outputs(metric, spatial_scale, approach, only_first_filename)

}

# List arguments of functions used in metrics from metric_files.yml
list_function_args <- function(metric = "all") {

  metric <- process_metric(metric = metric)
  system.file(
    "extdata",
    "metric_files.yml",
    package = "boundaries"
  ) %>%
    yaml::read_yaml() %>%
    get_function_args(metric)

}

# Translate metric options into internal metric names
process_metric <- function(metric = "all") {
  all_metrics <- c(
    "biome", "nitrogen", "lsc",
    "bluewater", "greenwater", "biosphere", "benchmark"
  )

  if ("all" %in% metric) {
    metric <- all_metrics
  }

  metric <- match.arg(
    arg = metric,
    choices = all_metrics,
    several.ok = TRUE
  )

  metric
}


# for input list a, all duplicate keys are unified, taking the value with
#     highest temporal resolution (daily>monthly>annual)
get_outputs <- function( # nolint
  x,
  metric_name,
  spatial_scale,
  approach,
  only_first_filename
) {
  outputs <- list()

  # Iterate over all metrics
  for (metric_string in names(x$metric[metric_name])) {
    metric <- x$metric[[metric_string]]

    # Iterate over all spatial scales
    for (scale_string in names(metric$spatial_scale)) { #nolint
      # Check if spatial scale is defined or all scales
      if (spatial_scale != "all" && !scale_string %in% spatial_scale) {
        next
      }
      scale <- metric$spatial_scale[[scale_string]]

      # Iterate over all approaches
      for (method_string in names(scale)) {
        # Check if approach is in list or all approaches
        if (all(approach != "all") && !is.null(approach[[metric_string]]) &&
              method_string != approach[[metric_string]]) { # nolint
          next
        }
        method <- scale[[method_string]]

        # Iterate over all outputs
        for (item in names(method$output)) {
          # Check if output is already in list or if it has higher resolution
          if (!item %in% names(outputs) ||
              (item %in% names(outputs) &&
                 higher_res(metric$output[[item]]$resolution,
                            outputs[[item]]$resolution))
          ) {
            # Assign output resolution from metric file
            outputs[[item]]$resolution <- method$output[[item]]$resolution #nolint
            outputs[[item]]$optional <- method$output[[item]]$optional #nolint
            # Assign output file name from metric file
            if (only_first_filename) {
              outputs[[item]]$file_name <- x$file_name[[item]][1]
            } else {
              outputs[[item]]$file_name <- x$file_name[[item]]
            }
          }
        }
      }
    }
  }
  outputs
}


# Get arguments of functions used in metrics
get_function_args <- function(x, metric_name) {
  # List functions of metrics (metric_name)
  funs <- list()

  for (metric in x$metric[metric_name]) {
    funs[[metric$fun_name]] <- metric$funs
  }

  # Get arguments of functions
  funs %>%
    lapply(function(x) {
      unlist(
        lapply(mget(x, inherits = TRUE), methods::formalArgs),
        use.names = FALSE
      )
    })
}


# Check if resolution of x is higher than resolution of y
higher_res <- function(x, y) {
  levels <- c("annual", "monthly", "daily")
  resolution_x <- match(match.arg(x, levels), levels)
  resolution_y <- match(match.arg(y, levels), levels)

  if (resolution_x > resolution_y) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

list_thresholds <- function(metric, approach, spatial_scale) {
  metric <- process_metric(metric = metric)

  yaml_data <- system.file(
    "extdata",
    "metric_files.yml",
    package = "boundaries"
  ) %>%
    yaml::read_yaml()

  return(yaml_data$metric[[metric]]$spatial_scale[[spatial_scale]][[approach]]$threshold) # nolint:line_length_linter

}

# set attributes for a control variable based on metric_files.yml and
#     user-defined thresholds. If thresholds is not defined, and x does not have
#     a thresholds attribute yet or overwrite is TRUE, the thresholds attribute
#     is set to the user-defined thresholds.
set_attributes <- function(
  x,
  approach = NULL,
  metric = NULL, # "bluewater", "greenwater", "biosphere", "lsc", or "nitrogen"
  spatial_scale = NULL, # "global", "regional", or "grid"
  thresholds = NULL,
  overwrite = FALSE
) {
  if (all(!is.null(c(approach, metric, spatial_scale)))) {
    yaml_data <- system.file(
      "extdata",
      "metric_files.yml",
      package = "boundaries"
    ) %>%
      yaml::read_yaml()

    relevant_data <- yaml_data$metric[[metric]]$spatial_scale[[spatial_scale]][[approach]] # nolint:line_length_linter

    # define the name of the control variable
    attr(x, "control_variable") <- relevant_data$control_variable

    # define the unit of the control variable
    attr(x, "unit") <- relevant_data$unit

    # define long name of the boundary
    attr(x, "long_name") <- yaml_data$metric[[metric]]$long_name
  }

  if (!is.null(spatial_scale)) {
    # define the spatial scale of the control variable
    attr(x, "spatial_scale") <- spatial_scale
  }

  if (!is.null(thresholds)) {
    # define the thresholds of the control variable
    if (is.null(attr(x, "thresholds")) ||
          overwrite == TRUE) {
      attr(x, "thresholds") <- thresholds
    }
  }

  class(x) <- c("control_variable")
  return(x)
}

# set user-defined attributes for a control variable, that was not computed with
#     the boundaries package (to make external data compatible with the
#     control_variable class needed for applying the plotting functions)
set_attributes_ext <- function(
  x,
  spatial_scale = NULL,
  thresholds = NULL,
  long_name = NULL,
  control_variable = NULL,
  unit = NULL
) {

  if (!is.null(control_variable)) {
    # define the name of the control variable
    attr(control_variable, "control_variable") <- control_variable
  }
  if (!is.null(unit)) {
    # define the unit of the control variable
    attr(control_variable, "unit") <- unit
  }
  if (!is.null(long_name)) {
    # define long name of the boundary
    attr(control_variable, "long_name") <- long_name
  }
  if (!is.null(spatial_scale)) {
    # define the spatial scale of the control variable
    attr(control_variable, "spatial_scale") <- spatial_scale
  }
  if (!is.null(thresholds)) {
    # define the thresholds of the control variable
    attr(control_variable, "thresholds") <- thresholds
  }

  class(control_variable) <- c("control_variable")
  return(control_variable)
}
