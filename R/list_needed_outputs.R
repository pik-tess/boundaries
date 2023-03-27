# output timesteps have to be added (monthly / yearly)
# differentiation between required and optional outputs (prec, temp)

#' Returns LPJmL outputs required for given metric
#'
#' Function to return a list of strings with the LPJmL output names required for
#' the computation of the metrics M-COL, M-ECO, the biome classification
#' and/or planetary boundary calculations
#'
#' @param metric string/list of strings, containing name of metric to get
#'        required outputs for. Pick from "meco", "mcol", "biome", "nitrogen",
#'      "lsc", "bluewater", "greenwater", "water", "biosphere", "all" (default)
#'
#' @param withNitrogen logical: include nitrogen outputs? default: TRUE
#'
#' @return List of strings with LPJmL variable names
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
list_needed_outputs <- function(metric = "all",
                                with_nitrogen = TRUE,
                                only_first_filename = TRUE) {

  optional_metrics <- c("meco", "mcol", "biome", "nitrogen", "lsc", "benchmark",
                        "bluewater", "greenwater", "water", "biosphere", "all")
  notin <- metric[!metric %in% optional_metrics]
  if (length(notin) > 0) {
    stop(paste("Metrics not available:", notin, collapse = ", "))
  }

  requirements <- system.file("extdata",
                              "metric_files.yaml",
                              package = "boundaries") %>%
    yaml::read_yaml()

  outs <- c()
  if ("meco" %in% metric) {
    if (with_nitrogen) {
      outs <- c(outs, requirements[["meco"]], requirements[["meco_nitrogen"]])
    } else {
      outs <- c(outs, requirements[["meco"]])
    }
  }
  if ("mcol" %in% metric) {
    outs <- c(outs, requirements[["mcol"]])
  }
  if ("biome" %in% metric) {
    outs <- c(outs, requirements[["biome"]])
  }
  if ("nitrogen" %in% metric) {
    if (!with_nitrogen)
      stop("You requested the nitrogen boundary without nitrogen?! Aborting.")
    outs <- c(outs, requirements[["nitrogen"]])
  }
  if ("lsc" %in% metric) {
    outs <- c(outs, requirements[["lsc"]])
  }
  if ("bluewater" %in% metric) {
    outs <- c(outs, requirements[["bluewater"]])
  }
  if ("greenwater" %in% metric) {
    outs <- c(outs, requirements[["greenwater"]])
  }
  if ("water" %in% metric) {
    outs <- c(outs, requirements[["bluewater"]], requirements[["greenwater"]])
  }
  if ("biosphere" %in% metric) {
    outs <- c(outs, requirements[["mcol"]], requirements[["meco"]])
  }
  if ("all" %in% metric) {
    if (with_nitrogen) {
      outs <- c(outs, requirements[["meco"]], requirements[["mcol"]],
                requirements[["meco_nitrogen"]], requirements[["nitrogen"]],
                requirements[["biome"]], requirements[["bluewater"]],
                requirements[["greenwater"]], requirements[["lsc"]])
    } else {
      outs <- c(outs, requirements[["meco"]], requirements[["mcol"]],
                requirements[["biome"]], requirements[["bluewater"]],
                requirements[["greenwater"]], requirements[["lsc"]])
    }
  }
  if ("benchmark" %in% metric) {
    outs <- c(requirements[["benchmark"]])
  }
  out <- unify_list(outs)

  if (only_first_filename) {
    output_filenames <- lapply(out, function(x) x[["file_name"]][1])
  } else{
    output_filenames <- lapply(out, function(x) x[["file_name"]])
  }
  return(list(outputs = output_filenames,
              timesteps = sapply(out, function(x) x[["resolution"]])))
}


list_outputs <- function(metric = "all",
                        with_nitrogen = TRUE,
                        only_first_filename = TRUE) {

  metric <- process_metric(metric = metric,
                           with_nitrogen = TRUE)

  system.file(
    "extdata",
    "metric_filefun.yml",
    package = "boundaries"
  ) %>%
    yaml::read_yaml() %>%
    get_outputs(metric, only_first_filename)

}


list_functions <- function(metric = "all",
                           with_nitrogen = TRUE) {

  metric <- process_metric(metric = metric,
                           with_nitrogen = TRUE)

  system.file(
    "extdata",
    "metric_filefun.yml",
    package = "boundaries"
  ) %>%
    yaml::read_yaml() %>%
    get_function_args(metric)

}

process_metric <- function(metric = "all",
                           with_nitrogen = TRUE) {
  all_metrics <- c(
    "meco", "mcol", "biome", "nitrogen", "lsc",
    "bluewater", "greenwater", "water", "biosphere"
  )

  if ("all" %in% metric) {
    metric <- all_metrics
  }

  metric <- match.arg(
    arg = metric,
    choices = all_metrics,
    several.ok = TRUE
  )

  if (with_nitrogen && "meco" %in% metric) {
    metric[metric == "meco"] <- "meco_nitrogen"
  }

  metric
}


# for input list a, all duplicate keys are unified, taking the value with
#     highest temporal resolution (daily>monthly>annual)
get_outputs <- function(x, metric_name, only_first_filename) {

  outputs <- list()
  for (metric in x$metric[metric_name]) {

    # Iterate over all unique keys
    for (item in names(metric$output)) {

      # Check if output is already in list or if it has higher resolution
      if (!item %in% names(outputs) ||
          (item %in% names(outputs) &&
          higher_res(metric$output[[item]], outputs[[item]]$resolution))) {

        outputs[[item]]$resolution <- metric$output[[item]]

        if (only_first_filename) {
          outputs[[item]]$file_name <- x$file_name[[item]][1]
        } else {
          outputs[[item]]$file_name <- x$file_name[[item]]
        }
      }
    }
  }
  return(outputs)
}


get_function_args <- function(x, metric_name) {

  # List functions of metrics (metric_name)
  funs <- c()
  for (metric in x$metric[metric_name]) {
    funs <- append(funs, unlist(metric$fun))
  }
  # Get arguments of functions
  funs %>%
    unique() %>%
    sapply(function(x) formalArgs(get(x)))
}



# for input list a, all duplicate keys are unified, taking the value with
#     highest temporal resolution (daily>monthly>annual)
unify_list <- function(a) {
  merged_list <- list()
  # Iterate over all unique keys
  for (item in unique(names(a))){
    indices <- which(names(a) == item)
    # Get the values for the current key
    values <- sapply(a[indices], function(x) x[["resolution"]])
    outval <- a[indices[1]]
    outval[[item]]$resolution <- highest_temporal_res(values)
    merged_list <- c(merged_list, outval)
  }
  return(merged_list)
}

# Among list of input values a, return that with highest temporal resolution
# (daily>monthly>annual)
highest_temporal_res <- function(a) {
  if ("daily" %in% a) {
    return("daily")
  } else if ("monthly" %in% a) {
    return("monthly")
  } else {
    return("annual")
  }
}


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
