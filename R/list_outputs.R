#' List required LPJmL outputs and temporal resolution
#'
#' Function to return a list of output IDs with required resolution and
#' file names for a given metric. The list is based on the metric_files.yml file
#' in the boundaries package (`"./inst/metric_files.yml"`).
#'
#' @param metric Character string containing name of metric to get
#'        required outputs. Available options are `c("meco", "mcol", "biome",
#'        "nitrogen", "lsc", "bluewater", "greenwater", "water", "biosphere")`
#'        or just `"all"` or `"benchmark"`.
#'
#' @param only_first_filename Logical. If TRUE, only the first file name will be
#'        returned for each output. If FALSE, all file names will be returned.
#'
#' @return List of output IDs with required resolution and
#'         file names for a given metric
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
list_outputs <- function(metric = "all", method, spatial_scale,
                         only_first_filename = TRUE) {

  metric <- process_metric(metric = metric)

  system.file(
    "extdata",
    "metric_files.yml",
    package = "boundaries"
  ) %>%
    yaml::read_yaml() %>%
    get_outputs(metric, method, spatial_scale, only_first_filename)

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
    "bluewater", "greenwater", "water", "biosphere"
  )

  if ("all" %in% metric) {
    metric <- all_metrics
  }

  if ("benchmark" %in% metric) {
    metric <- "benchmark"
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
get_outputs <- function(x, metric_name, method, spatial_scale,
                        only_first_filename) {

  outputs <- list()
  # Iterate over all metrics
  i <- 0
  for (metric in x$metric[metric_name]) {
    i <- i + 1
    boundary_name <- metric_name[i]
    if (length(method[[boundary_name]]) > 0) {
      method_i <- method[[boundary_name]]
    } else {
      method_i <- formals(get(paste0("calc_", boundary_name, "_status")))$method
    }
    if (length(spatial_scale) > 0) {
      spatial_scale_i <- spatial_scale
    } else {
      spatial_scale_i <- formals(get(paste0("calc_", boundary_name,
                                                 "_status")))$spatial_scale
    }
    # Iterate over all unique keys
    for (item in names(metric$spatial_scale[[spatial_scale_i]][[method_i]]$output)) { #nolint

      # Check if output is already in list or if it has higher resolution
      if (!item %in% names(outputs) ||
          (item %in% names(outputs) &&
           higher_res(metric$output[[item]]$resolution,
                      outputs[[item]]$resolution))
      ) {
        # Assign output resolution from metric file
        outputs[[item]]$resolution <- metric$spatial_scale[[spatial_scale_i]][[method_i]]$output[[item]]$resolution #nolint
        outputs[[item]]$optional <- metric$spatial_scale[[spatial_scale_i]][[method_i]]$output[[item]]$optional #nolint
        # Assign output file name from metric file
        if (only_first_filename) {
          outputs[[item]]$file_name <- x$file_name[[item]][1]
        } else {
          outputs[[item]]$file_name <- x$file_name[[item]]
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
      unlist(lapply(mget(x, inherits = TRUE), formalArgs), use.names = FALSE)
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

list_thresholds <- function(metric, method, spatial_scale) {
    metric <- process_metric(metric = metric)

    yaml_data <- system.file(
      "extdata",
      "metric_files.yml",
      package = "boundaries"
    ) %>%
      yaml::read_yaml()

    return(yaml_data$metric[[metric]]$spatial_scale[[spatial_scale]][[method]]$threshold)

}