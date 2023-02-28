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
list_needed_outputs <- function(metric = "all", with_nitrogen = TRUE) {

  optional_metrics <- c("meco", "mcol", "biome", "nitrogen", "lsc",
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
    outs <- c(outs,requirements[["biome"]])
  }
  if ("nitrogen" %in% metric) {
    if (!with_nitrogen) stop("You requested the nitrogen boundary without nitrogen?! Aborting.")
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
  out <- unify_list(outs)
  return(list(outputs = sapply(out,function(x) x[["file_name"]][1]), timesteps = sapply(out,function(x) x[["resolution"]])))
}

# for input list a, all duplicate keys are unified, taking the value with
#     highest temporal resolution (daily>monthly>annual)
unify_list <- function(a) {
  merged_list <- list()
  for (item in unique(names(a))){ # iterate over all unique keys
    indices <- which(names(a) == item)
    values <- sapply(a[indices],function(x) x[["resolution"]]) # get the values for the current key
    outval <- a[indices[1]]
    outval[[item]]$resolution <- highest_temporal_res(values)
    merged_list <- c(merged_list, outval)
  }
  return(merged_list)
}

# among list of input values a, return that with highest temporal resolution (daily>monthly>annual) # nolint
highest_temporal_res <- function(a) {
  if ("daily" %in% a) {
    return("daily")
  } else if ("monthly" %in% a) {
    return("monthly")
  } else {
    return("annual")
  }
}
