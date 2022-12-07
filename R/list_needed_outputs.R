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

  optional_metrics <- c("meco", "mcol", "biome", "all_pbs", "pb_n",
                        "pb_w", "pb_b", "pb_lsc", "all")
  notin <- metric[!metric %in% optional_metrics]
  if (length(notin) > 0) {
    stop(paste("Metrics not available:", notin, collapse = ", "))
  }
  varsGAMMA <- c("grid", "fpc", "fpc_bft", "cftfrac", "firec", "rh_harvest",
                 "npp", "runoff", "transp", "vegc", "firef", "rh", "harvestc",
                 "evap", "interc", "soilc", "litc", "swc")
  varsGAMMAnitrogen <- c("vegn", "soilnh4", "soilno3", "leaching", "n2o_denit",
                         "n2o_nit", "n2o_denit", "n2_emis", "bnf",
                         "n_volatilization")
  varsHANPP <- c("grid", "npp", "pft_npp", "pft_harvestc", "pft_rharvestc",
                 "firec", "timber_harvestc", "cftfrac", "fpc")
  varsBiome <- c("grid", "fpc", "vegc", "pft_lai", "temp")
  vars_pb_n <- c("grid", "runoff", "leaching", "pet", "prec")
  vars_pb_w <- c("grid", "discharge", "irrig", "drainage")
  vars_pb_lsc <- varsBiome
  vars_pb_b <- c()

  outs <- c()
  tsteps <- c() # todo: add timestep functionality - needs to be clever

  if ("all" %in% metric) {
    if (with_nitrogen) {
      outs <- c(outs,varsGAMMA, varsGAMMAnitrogen, varsHANPP, varsBiome,
                           vars_pb_n, vars_pb_w, vars_pb_lsc, vars_pb_b)
    } else {
      outs <- c(outs,varsGAMMA, varsHANPP, varsBiome, vars_pb_w, vars_pb_lsc,
                vars_pb_b)
    }
  }
  if ("meco" %in% metric) {
    if (with_nitrogen) {
      outs <- c(outs,varsGAMMA, varsGAMMAnitrogen)
    } else {
      outs <- c(outs,varsGAMMA)
    }
  }
  if ("biome" %in% metric) {
    outs <- c(outs,varsBiome)
  }
  if ("mcol" %in% metric) {
    outs <- c(outs,varsHANPP)
  }
  if ("all_pbs" %in% metric) {
    if (with_nitrogen) {
      outs <- c(outs,varsBiome, vars_pb_n, vars_pb_w, vars_pb_lsc, vars_pb_b)
    } else {
      outs <- c(outs,varsBiome, vars_pb_w, vars_pb_lsc, vars_pb_b)
    }
  }
  if ("pb_n" %in% metric) {
    outs <- c(outs,vars_pb_n)
  }
  if ("pb_w" %in% metric) {
    outs <- c(outs,vars_pb_w)
  }
  if ("pb_lsc" %in% metric) {
    outs <- c(outs,vars_pb_lsc)
  }
  if ("pb_b" %in% metric ) {
    outs <- c(outs,vars_pb_b)
  }
  return(list(outputs = unique(sort(outs)), timesteps = tsteps))
}
