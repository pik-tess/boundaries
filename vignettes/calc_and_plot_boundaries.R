## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo = TRUE, eval = FALSE, highlight = TRUE------------------------------
#  global_status <- calc_status(
#    boundary = c("bluewater", "greenwater", "lsc", "nitrogen", "biosphere"),
#    config_scenario = config_path_scenario,
#    config_reference = config_path_reference,
#    spatial_scale = "global",
#    path_baseline = output_path_baseline,
#    time_span_scenario = as.character(1901:2017),
#    time_span_reference = as.character(1500:1699),
#    time_span_baseline = as.character(1901:2017),
#    time_series_avg = 1, # computation of timeseries, no moving average
#    approach = list("bluewater" = "porkka2024",
#                    "greenwater" = "porkka2024",
#                    "nitrogen" = "schulte_uebbing2022")
#  )

## ----echo = TRUE, eval = FALSE, highlight = TRUE------------------------------
#  plot_status(global_status, filename = "./global_status.png")
#  status_legend(filename = "./legend.png")

## ----echo = TRUE, eval = FALSE, highlight = TRUE------------------------------
#  validation_table <- validate_simulation(
#    config_scenario = config_path_scenario,
#    config_reference = config_path_reference,
#    path_baseline = output_path_baseline,
#    time_span_scenario = as.character(2005:2014),
#    time_span_reference = as.character(1500:1699),
#    filename = "./validation_table.csv",
#    approach = list("bluewater" = "porkka2024",
#                    "greenwater" = "porkka2024",
#                    "nitrogen" = "schulte_uebbing2022")
#  )

## ----setup, echo = TRUE, eval = FALSE, highlight = TRUE-----------------------
#  library(boundaries) # nolint:undesirable_function_linter
#  
#  # set paths to lpjml configuration files and outputs
#  path_configurations <- "./lpjml/configurations/"
#  config_path_scenario <- paste0(path_configurations, "config_scenario.json")
#  config_path_reference <- paste0(path_configurations, "config_reference.json")
#  path_outputs <- "./lpjml/output/"
#  path_outputs_pnv <- paste0(path_outputs, "pnv/")
#  
#  plotpath <- paste0("./R/plots/")
#  

## ----echo = TRUE, eval = FALSE, highlight = TRUE------------------------------
#  ### calc status at the global level
#  
#  # define analysis and reference period
#  time_span_scenario <- as.character(1901:2017)
#  time_span_reference <- as.character(1500:1699)
#  
#  # define the number of years to average over, for a moving average; if set to
#  # 1, no moving average is calculated
#  nyear_window <- 5
#  
#  global_status <- calc_status(
#    # define boundaries to calculate - it can also be only one boundary
#    boundary = c("bluewater", "greenwater", "lsc", "nitrogen", "biosphere"),
#    config_scenario = config_path_scenario,
#    config_reference = config_path_reference,
#    spatial_scale = "global",
#    time_span_scenario = time_span_scenario,
#    time_span_reference = time_span_reference,
#    # moving average over nyear_window years
#    time_series_avg = nyear_window,
#    approach = list("bluewater" = "porkka2024",
#                    "greenwater" = "porkka2024",
#                    "nitrogen" = "schulte_uebbing2022"),
#    # boundary specific parameters, see individual boundary functions for details
#    time_span_baseline = time_span_scenario, # for biosphere integrity
#    path_baseline = paste0(path_outputs, "/pnv/"), # for biosphere integrity
#    savanna_proxy = list(vegc = 7500) # for forest biome definition in the land system change boundary
#  )
#  
#  ### plot timeseries
#  
#  # There are three options for plotting the global status timeseries:
#  # a) one timeseries panel for each boundary control variable, in a Cartesian
#  #    coordinate system
#  # b) all timeseries plotted in one panel, in a Cartesian coordinate system
#  # c) one timeseries panel for each boundary control variable, in a polar
#  #    coordinate system (called "stylized")
#  
#  # a) one timeseries panel for each boundary control variable
#  plot_status(
#    x = global_status,
#    stylized = FALSE,
#    filename = "./global_timeseries_panels.png"
#  )
#  
#  # b) all timeseries plotted in one panel
#  plot_status(
#    x = global_status,
#    stylized = FALSE,
#    all_in_one = TRUE,
#    filename = "./global_timeseries_all_in_one.png"
#  )
#  
#  # c) one timeseries panel for each boundary control variable, in a polar
#  #    coordinate system
#  plot_status(
#    x = global_status,
#    stylized = TRUE,
#    filename = "./global_timeseries_stylized.png"
#  )
#  

## ----echo = TRUE, eval = FALSE, highlight = TRUE------------------------------
#  # In this example, only the bluewater boundary is calculated, for two different
#  # approaches. The results are then jointly plotted in a Cartesian coordinate
#  # system.
#  
#  # define analysis and reference period
#  time_span_scenario <- as.character(1901:2017)
#  time_span_reference <- as.character(1500:1699)
#  
#  # Aprroach A: Following Porkka et al. (2024) (https://doi.org/10.1038/s44221-024-00208-7)
#  # Referring to global area with discharge deviations outside the pre-industrial
#  # range
#  
#  bluewater_status_porkka <- calc_status(
#    # define boundaries to calculate - it can also be only one boundary
#    boundary = c("bluewater"),
#    config_scenario = config_path_scenario,
#    config_reference = config_path_reference,
#    spatial_scale = "global",
#    time_span_scenario = time_span_scenario,
#    time_span_reference = time_span_reference,
#    # no moving average, by setting the number of years to average over to 1:
#    time_series_avg = 1,
#    approach = list("bluewater" = "porkka2024")
#  )
#  
#  # Aprroach B: Following Rockström et al. (2009) (https://doi.org/10.1038/461472a)
#  # Referring to global bluewater consumption, but adapting the defaullt boundary
#  # and high risk values
#  
#  bluewater_status_rockstroem <- calc_status(
#    # define boundaries to calculate - it can also be only one boundary
#    boundary = c("bluewater"),
#    config_scenario = config_path_scenario,
#    config_reference = config_path_reference,
#    spatial_scale = "global",
#    time_span_scenario = time_span_scenario,
#    time_span_reference = time_span_reference,
#    # no moving average, by setting the number of years to average over to 1:
#    time_series_avg = 1,
#    approach = list("bluewater" = "rockstroem2009"),
#    # change default planetary boundary and high risk values, following Gerten
#    # et al. 2013 (https://doi.org/10.1016/j.cosust.2013.11.001),
#    # all referring to km3/yr as defined in `metric_files.yml`
#    thresholds = list("bluewater" = list(holocene = 0,
#                                         pb = 2800,
#                                         high_risk = 4000))
#  )
#  
#  # plot both results in one timeseries panel
#  # merge both results into one list
#  # TODO test if this is working!
#  bluewater_status <- list("bluewater" = bluewater_status_porkka$bluewater,
#                           "bluewater" = bluewater_status_rockstroem$bluewater)
#  
#  plot_status(
#    x = bluewater_status,
#    stylized = FALSE,
#    all_in_one = TRUE,
#    filename = "./global_bluewater_status_porkka_vs_rockstroem.png"
#  )
#  

## ----echo = TRUE, eval = FALSE, highlight = TRUE------------------------------
#  
#  
#  # define analysis and reference period
#  time_span_scenario <- as.character(2008:2017)
#  time_span_reference <- as.character(1500:1699)
#  # if the land use effect is to be isolated, the time span for the reference
#  # period can be set to the same as the scenario period:
#  # time_span_reference <- as.character(2008:2017) #nolint
#  # This way, climate change induced changes are excluded.
#  
#  # calc status at the gridded level
#  gridded_status <- calc_status(
#    boundary = c("lsc", "nitrogen", "bluewater", "biosphere"),
#    config_scenario = config_path_scenario,
#    config_reference = config_path_reference,
#    time_span_scenario = time_span_scenario,
#    time_span_reference = time_span_reference,
#    spatial_scale = "grid",
#    # set time_series_avg to NULL, to not calculate a timeseries, but to average
#    # over the entire scenario time span
#    time_series_avg = NULL,
#    # boundary specific parameters, see individual boundary functions for details
#    time_span_baseline = time_span_scenario, # for biosphere integrity
#    path_baseline = paste0(path_outputs, "/pnv/"), # for biosphere integrity
#    savanna_proxy = list(vegc = 7500) # for forest biome definition in the land system change boundary
#  )
#  
#  # plot status at the gridded level, for details see `?status_maps`
#  # There are two options for plotting a map with the gridded status(es):
#  # a) plot the control variable status of each boundary (e.g. deforestation share
#  #   for the land system change boundary)
#  # b) plot the risk level of each boundary, based on a normalized color scale
#  #    and the boundary and high risk values
#  
#  # a) control variable status
#  plot_status(
#    x = gridded_status,
#    filename = "./gridded_status_control_variable.png",
#    grid_path = paste0(path_outputs_pnv, "grid.bin.json"),
#    risk_level = FALSE
#  )
#  
#  # b) risk level
#  plot_status(
#    x = gridded_status,
#    filename = "./gridded_status_risk_level.png",
#    grid_path = paste0(path_outputs_pnv, "grid.bin.json"),
#    risk_level = TRUE
#  )
#  

## ----echo = TRUE, eval = FALSE, highlight = TRUE------------------------------
#  
#  # define analysis and reference period
#  time_span_scenario <- as.character(2008:2017)
#  time_span_reference <- as.character(1500:1699)
#  
#  # calc status at the regional level
#  regional_status <- calc_status(
#    boundary = c("lsc", "bluewater", "greenwater", "biosphere"),
#    config_scenario = config_path_scenario,
#    config_reference = config_path_scenario,
#    time_span_scenario = time_span_scenario,
#    time_span_reference = time_span_reference,
#    spatial_scale = "regional",
#    approach = list("bluewater" = "porkka2024",
#                    "greenwater" = "porkka2024"),
#    path_baseline = paste0(path_outputs, "/pnv/"), # for biosphere integrity
#    time_span_baseline = time_span_scenario, # for biosphere integrity
#    savanna_proxy = list(vegc = 7500), # for forest biome definition in the land system change boundary
#  )
#  
#  # plot status at the regional level
#  
#  # As for the gridded status, there are two options for plotting a map with the
#  # regional status(es): directly the control variable status or the risk level.
#  
#  # risk level status
#  plot_status(
#    x = regional_status,
#    filename = "./regional_status_risk_level.png",
#    grid_path = paste0(path_outputs_pnv, "grid.bin.json"),
#    risk_level = TRUE
#  )

