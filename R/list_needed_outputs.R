# output timesteps have to be added (monthly / yearly)
# differentiation between required and optional outputs (prec, temp)

#' Returns LPJmL outputs required for given metric
#'
#' Function to return a list of strings with the LPJmL output names required for
#' the computation of the metrics M-COL, M-ECO, the biome classification
#' and/or planetary boundary calculations
#'
#' @param metric string/list of strings, containing name of metric to get
#'        required outputs for. Pick from "meco", "mcol", "biome", "all_pbs",
#'               "pb_n", "pb_w", "pb_b", "pb_lsc", "all" (default)
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

  optional_metrics <- c("meco", "mcol", "biome", "all_pbs", "nitrogen",
                        "bluewater", "greenwater", "biosphere", "lsc", "all")
  notin <- metric[!metric %in% optional_metrics]
  if (length(notin) > 0) {
    stop(paste("Metrics not available:", notin, collapse = ", "))
  }
  # TODO: camelCase to snake_case
  varsGAMMA <- list(grid = "annual",
                    fpc = "annual",
                    fpc_bft = "annual",
                    cftfrac = "annual",
                    firec = "annual",
                    npp = "annual",
                    runoff = "annual",
                    transp = "annual",
                    vegc = "annual",
                    firef = "annual",
                    rh = "annual",
                    harvestc = "annual",
                    evap = "annual",
                    interc = "annual",
                    soilc = "annual",
                    litc = "annual",
                    swc = "annual")
  varsGAMMAnitrogen <- list(vegn = "annual",
                            soilnh4 = "annual",
                            soilno3 = "annual",
                            leaching = "annual",
                            n2o_denit = "annual",
                            n2o_nit = "annual",
                            n2o_denit = "annual",
                            n2_emis = "annual",
                            bnf = "annual",
                            n_volatilization = "annual")
  varsHANPP <- list(grid = "annual",
                    npp = "annual",
                    pft_npp = "annual",
                    pft_harvestc = "annual",
                    pft_rharvestc = "annual",
                    firec = "annual",
                    timber_harvestc = "annual",
                    cftfrac = "annual",
                    fpc = "annual")
  varsBiome <- list(grid = "annual",
                    fpc = "annual",
                    vegc = "annual",
                    pft_lai = "annual",
                    temp = "annual")
  vars_pb_n <- list(grid = "annual",
                    runoff = "annual",
                    leaching = "annual",
                    pet = "annual",
                    prec = "annual")
  vars_pb_w <- list(grid = "annual",
                    discharge = "monthly",
                    irrig = "monthly",
                    drainage = "annual")
  vars_pb_lsc <- varsBiome
  vars_pb_b <- c(varsGAMMA, varsHANPP)
  if (with_nitrogen) vars_pb_b <- c(vars_pb_b, varsGAMMAnitrogen)

  outs <- c()

  if ("all" %in% metric) {
    if (with_nitrogen) {
      outs <- c(outs, varsGAMMA, varsGAMMAnitrogen, varsHANPP, varsBiome,
                           vars_pb_n, vars_pb_w, vars_pb_lsc, vars_pb_b)
    } else {
      outs <- c(outs, varsGAMMA, varsHANPP, varsBiome, vars_pb_w, vars_pb_lsc,
                vars_pb_b)
    }
  }
  if ("meco" %in% metric) {
    if (with_nitrogen) {
      outs <- c(outs, varsGAMMA, varsGAMMAnitrogen)
    } else {
      outs <- c(outs, varsGAMMA)
    }
  }
  if ("biome" %in% metric) {
    outs <- c(outs,varsBiome)
  }
  if ("mcol" %in% metric) {
    outs <- c(outs, varsHANPP)
  }
  if ("all_pbs" %in% metric) {
    if (with_nitrogen) {
      outs <- c(outs, varsBiome, vars_pb_n, vars_pb_w, vars_pb_lsc, vars_pb_b)
    } else {
      outs <- c(outs, varsBiome, vars_pb_w, vars_pb_lsc, vars_pb_b)
    }
  }
  if ("pb_n" %in% metric) {
    outs <- c(outs, vars_pb_n)
  }
  if ("pb_w" %in% metric) {
    outs <- c(outs, vars_pb_w)
  }
  if ("pb_lsc" %in% metric) {
    outs <- c(outs, vars_pb_lsc)
  }
  if ("pb_b" %in% metric) {
    outs <- c(outs, vars_pb_b)
  }
  out <- unify_list(outs)
  return(list(outputs = names(out), timesteps = unname(unlist(out))))
}

# for input list a, all duplicate keys are unified, taking the value with
#     highest temporal resolution (daily>monthly>annual)
unify_list <- function(a) {
  merged_list <- list()
  for (item in unique(names(a))){ # iterate over all unique keys
    values <- a[which(names(a) == item)] # get the values for the current key
    merged_list <- c(merged_list,
                     setNames(highest_temporal_res(values), item))
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