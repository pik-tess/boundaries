rm(list = ls(all = TRUE))

library(lpjmlkit)
# library(abind)
library(stringr)
library(dplyr)
# library(tidyr)

# TODO:
# 1. improve documentation and formatting. is there a standardised way?
# 2. link to .yaml file
# 3. improve readability (at least two lit values for each variables, cosider ranges, formattable?)
# 4. consistent rounding


# function to benchmark lpjml-based PB calculations to literature values
# currently still using explicit path (in JB folder), for the future use .yaml file
# what should be the name of the main function?
validation_table <- function(dir_output, table_path, timeframe, calculate_pb_status = TRUE) { #nolint
  # load literature table
  ref_table <- read.csv2(ref_table_path)

  # define set of helper functions
  # functions to deal with time and spacial scales of lpjml output
  global_sum_yearly <- function(file, timeframe) {
    output <- read_io(file, subset = list(year = as.character(timeframe))) %>%
      transform(to = "year_month_day") %>%
      as_array(aggregate = list(year = mean)) %>%
      sum(. * terr_area) * 10^-12 # why this factor?
    return(output)
  }

  grid_mean_yearly <- function(file, timeframe) {
    data <- read_io(file, subset = list(year = as.character(timeframe))) %>%
      transform(to = "year_month_day") %>%
      as_array(aggregate = list(year = mean)) %>%
      drop()
    return(data)
  }

  grid_sum_monthly_mean_yearly <- function(file, timeframe) {
    data <- read_io(file, subset = list(year = as.character(timeframe))) %>%
      transform(to = "year_month_day") %>%
      as_array(aggregate = list(
        year = mean,
        month = sum
      )) %>%
      drop()
    return(data)
  }

  # insert the computed variables in the literature table through matching by variable and unit patterns, make sure the structure of original .csv table is mantained
  insert_lpjml_value <- function(ref_table, var_patterns, unit_patterns, calculated_vars) {
    lpjml_values <- rep(NA_real_, nrow(ref_table))
    for (i in seq_along(var_patterns)) {
      matching_indices <- grepl(var_patterns[[i]], ref_table$variable, ignore.case = TRUE) &
        grepl(unit_patterns[[i]], ref_table$unit, ignore.case = TRUE)
      lpjml_values[matching_indices] <- round(eval(parse(text = calculated_vars[[i]]))[1], digits)
    }
    ref_table$lpjml_value <- lpjml_values
    return(ref_table)
  }

  # insert the PB status values (calculated with calc_status) in the benchmarking table
  evaluate_pb_status <- function(ref_table, pb_list, temp2) {
    for (pb in pb_list) { # add calculation of PBs
      ref_table <- ref_table %>%
        mutate(lpjml_value = ifelse(grepl(pb, variable, ignore.case = TRUE) & grepl("PB status", variable, ignore.case = TRUE),
          round(temp2[[pb]][1], digits = 4),
          lpjml_value
        ))
    }
    return(ref_table)
  }

  # end of helper function definition
  # read in terrestrial area
  terr_area <- read_io(paste0(dir_output, "/terr_area.bin.json"))$data %>% drop() # in m2

  cftfrac <- grid_mean_yearly(paste0(dir_output, "/cftfrac.bin.json"), timeframe)

  ncell <- length(terr_area)

  # PB Land-system change

  #-------------- total agricultural area --------------------------------------
  # better based on input than on output (to include areas with failed harvest)?

  area_cft <- apply(cftfrac * terr_area * 10^-10, c(2), sum) # in mio ha

  crop_area <- sum(area_cft[c(1:13, 17:29)])

  pasture_area <- sum(area_cft[c(14, 30)])

  irrig_area <- sum(area_cft[c(17:32)])

  #------------- potential forest extent ---------------------------------------

  fpc <- grid_mean_yearly(paste0(dir_output, "/fpc.bin.json"), timeframe)

  # ids for tree plant funtional types
  tree_ids <- c(2:9)
  fpc_tree_total <- apply(fpc[, tree_ids], c(1), sum)

  is_forest <- fpc_tree_total >= 0.6
  forest_area <- sum(terr_area[is_forest]) * 10^-10
  # TODO: take out savannas in the tropics?

  #------------- deforestation share ---------------------------------------------

  forest <- rep(0, ncell)
  forest[is_forest] <- 1
  deforest_area <- sum(apply(cftfrac, c(1), sum) * terr_area * forest) * 10^-10
  deforest_share <- deforest_area / forest_area * 100

  # PB freshwater change
  #------------- Irrigation withdrawals and consumption ------------------------


  irrig <- grid_sum_monthly_mean_yearly(paste0(dir_output, "/irrig.bin.json"), timeframe)
  conv_loss_evap <- grid_sum_monthly_mean_yearly(paste0(dir_output, "/aconv_loss_evap.bin.json"), timeframe)
  conv_loss_drain <- grid_sum_monthly_mean_yearly(paste0(dir_output, "/aconv_loss_drain.bin.json"), timeframe)
  return_flow_b <- grid_sum_monthly_mean_yearly(paste0(dir_output, "/mreturn_flow_b.bin.json"), timeframe)
  # withdrawals in km3
  wd <- sum((irrig + conv_loss_drain + conv_loss_evap) * terr_area) * 10^-12

  # consumption in km3
  cons <- sum((irrig + conv_loss_evap - return_flow_b) * terr_area) * 10^-12

  # PB Biogechemical flows
  #------------- Nitrogen balance and leaching ---------------------------------

  # N leaching (total) in Tg
  leaching <- global_sum_yearly(
    file = paste0(dir_output, "/mleaching.bin.json"),
    # timeframe = timeframe
    timeframe = as.character(timeframe)
  )

  # N leaching from cropland
  leaching_agr <- global_sum_yearly(
    file = paste0(
      dir_output,
      "/nleaching_agr.bin.json"
    ),
    timeframe = timeframe
  )

  # N use efficiency on cropland
  # inputs
  total_fert <- global_sum_yearly(
    file = paste0(dir_output, "/nfert_agr.bin.json"),
    timeframe
  )
  total_man <- global_sum_yearly(
    file = paste0(dir_output, "/nmanure_agr.bin.json"), # nolint
    timeframe
  )
  total_dep <- global_sum_yearly(
    file = paste0(dir_output, "/ndepo_agr.bin.json"),
    timeframe
  )
  total_bnf <- global_sum_yearly(
    file = paste0(dir_output, "/bnf_agr.bin.json"),
    timeframe
  )
  total_seed <- global_sum_yearly(
    file = paste0(dir_output, "seedn_agr.bin.json"),
    timeframe
  )

  n_inputs <- total_fert + total_man + total_dep + total_bnf + total_seed

  nharvest <- global_sum_yearly(
    file = paste0(
      dir_output,
      "harvestn_agr.bin.json"
    ),
    timeframe
  )

  nue <- nharvest / n_inputs * 100

  # N surplus on cropland
  nsurplus <- n_inputs - nharvest

  # PB biosphere integrity
  #------------- NPP -------------------------------------------
  # NPP LU
  npp_lu <- global_sum_yearly(
    file = paste0(
      dir_output,
      "mnpp.bin.json"
    ),
    timeframe
  ) * 10^-3
  #-------------- total crop production CFT1-12 in freshmatter--------------------

  # assumption: carbon content in dry matter = 0.45
  dm_factor <- 0.45

  # conversion factors from freshmatter to dry matter
  # in gDM gFM-1 (aus Wirsenius (2000): Human use of land and organic materials.
  # Modeling the turnover of biomass in the global food system))
  fm2dm <- c(
    0.88, 0.87, 0.88, 0.88, 0.90, 0.24, 0.35, 0.93, 0.91, 0.94, 0.92,
    0.27, NA, NA
  )
  dm2fm <- 1 / fm2dm

  ### simulated lpjml harvest
  yield <- grid_mean_yearly(paste0(dir_output, "/pft_harvest.pft.bin.json"), timeframe)

  # dry matter harvest g
  harvest_dm <- cftfrac * yield * terr_area * (1 / dm_factor)


  # freshmatter harvest
  crop_fm_prod <- array(0, dim = c(ncell, 12)) # fm production for all CFTs
  for (c in 1:12) {
    crop_fm_prod[, c] <- (harvest_dm[, c] + harvest_dm[, c + 16]) *
      dm2fm[c] / 1000 # kg FM
  }

  crop_sum <- sum(crop_fm_prod) / 10^9



  # lists to match literature character strings with the here computed variables
  var_match <- list("crop", "irrig", "pasture", "forest", "def", "withdr", "cons", "leaching (total)", "leaching from crop", "efficiency", "surplus", "npp", "prod")
  units_match <- list("Mha", "Mha", "Mha", "Mha", "%", "km3", "km3", "TgN", "TgN", "%", "TgN", "PgC", "Gt")
  calc_var <- list(crop_area, irrig_area, pasture_area, forest_area, deforest_share, wd, cons, leaching, leaching_agr, nue, nsurplus, npp_lu, crop_sum)

  # insert lpjml values in table (pb status variables later)
  ref_table <- insert_lpjml_value(ref_table, var_match, units_match, calc_var)

  # check where the match did not work (units of Nitrogen leaching)
  # ref_table %>%
  #   filter(!is.numeric(lpjml_value) | is.na(lpjml_value))

  # include PB status values
  if (calculate_pb_status == TRUE) {
    # define folders as parameters? or through .yaml
    # path_outputs <- "/p/projects/open/Johanna/boundaries/lpjml/configurations/"
    # path_reference <- paste0(path_outputs, "/config_pnv_1500_2017_mg.json")
    # path_scenario <- paste0(path_outputs, "/config_lu_1500_2017_mg.json")


    # calculate PB status
    # timeframe_reference <- as.character(1901:1930)

    # temp2 <- calc_status(boundary = c("lsc" ,"greenwater", "bluewater", "nitrogen"
    #                                   ),
    #   config_scenario = path_scenario,
    #   config_reference = path_reference,
    #   path_baseline = path_reference,
    #   time_span_scenario = timeframe,
    #   time_span_reference = timeframe_reference,
    #   spatial_scale = "global",
    #   with_groundwater_denit = FALSE,
    #   savanna_proxy = list(vegc = 7500),
    #   in_parallel = TRUE,
    #   gridbased = FALSE,
    #   method = list("bluewater" = "porkka2023",
    #                 "greenwater" = "porkka2023",
    #                 "nitrogen" = "schulte_uebbing2022"),
    # )
    # save calc_status output for later analysis
    # save(temp2, file= paste0(my_dir, "/benchmarking/scripts/data/calc_status_data.RData"))

    # load data from old calc_status
    load(file = paste0(my_dir, "/benchmarking/scripts/data/calc_status_data.RData"))

    # list of pb varaibles to insert in the table
    pb_list <- list("lsc", "greenwater", "bluewater", "nitrogen", "biosphere")

    # insert values from calc_status in the table
    ref_table <- evaluate_pb_status(ref_table, pb_list, temp2)
  } # option to calculate pb_status ends here



  # add column for normalized error values
  ref_table <- ref_table %>% mutate(norm_error = signif((lpjml_value - value) / value, digits = 2))

  # check missing values
  # ref_table %>%
  # filter(!is.numeric(lpjml_value) | is.na(lpjml_value))
  return(ref_table)
}

# test
my_dir <- getwd()
dir_output <- "/p/projects/open/Johanna/boundaries/lpjml/output/lu_1901_2017_mg/"
ref_table_path <- "/p/projects/open/Caterina/benchmarking/scripts/data/reference_list.csv"
timeframe <- as.character(c(2010:2017))
benchmark_pb_table <- validation_table(dir_output, table_path, timeframe, calculate_pb_status = TRUE)
