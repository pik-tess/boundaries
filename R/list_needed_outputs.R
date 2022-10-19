# output timesteps have to be added (monthly / yearly)
# differentiation between required and optional outputs (prec, temp)

#' Returns LPJmL outputs required for given metric
#'
#' Function to return a list of strings with the LPJmL output names required for
#' the computation of the metrics M-COL, M-ECO, the biome classification
#' and/or planetary boundary calculations
#'
#' @param metric string: which metric? "meco", "mcol", "biome", "all_pbs",
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
list_needed_outputs <- function(metric="all", with_nitrogen = TRUE) {

  optional_metrics <- c("meco", "mcol", "biome", "all_pbs", "pb_n",
                        "pb_w", "pb_b", "pb_lsc", "all")
  metrics <- match.arg(metric, optional_metrics)

  varsGAMMA <- c("grid", "fpc", "fpc_bft", "cftfrac", "firec", "rh_harvest",
                 "npp", "runoff", "transp", "vegc", "firef", "rh", "harvestc",
                 "evap", "interc", "soilc", "litc", "swc")
  varsGAMMAnitrogen <- c("vegn", "soilnh4", "soilno3", "leaching", "n2o_denit",
                         "n2o_nit", "n2o_denit", "n2_emis", "bnf",
                         "n_volatilization")
  varsHANPP <- c("grid", "mnpp", "pft_npp", "pft_harvest", "pft_rharvest",
                 "firec", "timber_harvest", "cftfrac", "fpc")
  varsBiome <- c("grid", "fpc", "vegc", "pft_lai", "temp")
  vars_pb_n <- c("grid", "runoff", "leaching", "pet", "prec")
  vars_pb_w <- c("grid", "discharge", "irrig", "drainage")
  vars_pb_lsc <- varsBiome
  vars_pb_b <- c()

  if (metric == "all") {
    if (with_nitrogen) {
      return(unique(sort(c(varsGAMMA, varsGAMMAnitrogen, varsHANPP, varsBiome,
                           vars_pb_n, vars_pb_w, vars_pb_lsc, vars_pb_b))))
    } else {
      return(unique(sort(c(varsGAMMA, varsHANPP, varsBiome, vars_pb_w,
                           vars_pb_lsc, vars_pb_b))))
    }
  } else if (metric == "meco") {
    if (with_nitrogen) {
      return(unique(sort(c(varsGAMMA, varsGAMMAnitrogen))))
    } else {
      return(unique(sort(varsGAMMA)))
    }
  } else if (metric == "biome") {
    return(unique(sort(varsBiome)))
  } else if (metric == "mcol") {
    return(unique(sort(varsHANPP)))
  } else if (metric == "all_pbs") {
    if (with_nitrogen) {
      return(unique(sort(varsBiome, vars_pb_n, vars_pb_w, vars_pb_lsc,
                         vars_pb_b)))
    } else {
      return(unique(sort(varsBiome, vars_pb_w, vars_pb_lsc, vars_pb_b)))
    }
  } else if (metric == "pb_n") {
    return(unique(sort(vars_pb_n)))
  } else if (metric == "pb_w") {
    return(unique(sort(vars_pb_w)))
  } else if (metric == "pb_lsc") {
    return(unique(sort(vars_pb_lsc)))
  } else if (metric == "pb_b") {
    return(unique(sort(vars_pb_b)))
  }
}
