#' Validate simulated global PB-relevant variables against literature

#' Calculate a table with global modelled vs literature values for key variables
#' relevant to planetary boundaries
#'
#' @param config_scenario character string. File path to the LPjmL configuration
#' file (json) of the scenario run. The configuration file contains the
#' information about the LPJmL run, e.g. the output directory
#'
#' @param config_reference character string. See config_scenario. For the
#' reference run
#'
#' @param time_span_scenario time span to be used for the scenario run and
#' parallel PNV run, defined as a character string,
#' e.g. `as.character(1982:2011)`
#'
#' @param time_span_reference time span to be used for the reference run, defined
#' as an integer vector, e.g. `1901:1930`. Can differ in offset and length from
#' `time_span_scenario`! If `NULL` value of `time_span_scenario` is used
#'
#' @param path_baseline character string for path to outputs
#' for the baseline run, file names are taken from files scenario
#'
#' @param table_path character string for file path to save the output,
#' (.csv file)
#'
#' @param ... arguments to be passed to [calc_status]
#'
#' @return  table with comparison between lpjml values and literature ranges
#'
#' @examples
#' \dontrun{
#' validation <- validation_table(
#'   config_scenario = "./my_path/config_scenario.json",
#'   config_reference = "./my_path/config_reference.json",
#'   time_span_scenario = as.character(2010:2017),
#'   time_span_reference = as.character(1500:1699),
#'   path_baseline = "./my_path/outputs_baseline/",
#'   table_path = "./my_path/table.csv"
#' )
#' }
#'
#' @md
#' @export
#'

validation_table <- function(
    config_scenario,
    config_reference,
    time_span_scenario,
    time_span_reference,
    path_baseline,
    table_path,
    ...) {

  # use future package for asynchronous parallelization

  rlang::local_options(
    future.globals.maxSize = 3000 * 1024^2
  )
  if (.Platform$OS.type == "windows") {
    future_plan <- future::plan("multisession")
  } else {
    future_plan <- future::plan("multicore")
  }
  on.exit(future::plan(future_plan)) # nolint:undesirable_function_linter


  # please R CMD check for use of dplyr syntax
  lpjml_value <- value <- variable <- range.lower <- range.upper <- unit <- NULL
  X.lower. <- X.upper. <- year <- time_range <- value_range <- Author.s. <- NULL # nolint
  assessment_time <- year.of.publication <- boundary <- literature.range <- NULL # nolint

  config_scen <- lpjmlkit::read_config(config_scenario)
  config_ref <- lpjmlkit::read_config(config_reference)

  # are pft specific outputs grid or stand scaled?
  gridbased <- config_scen$grid_scaled

  if (!all(time_span_scenario %in% get_sim_time(config_scen))) {
    stop("Time span not available in scenario run.")
  }
  if (!all(time_span_reference %in% get_sim_time(config_ref))) {
    stop("Time span not available in reference run.")
  }

  # List required output files for each boundary
  output_files <- list_outputs(
    "benchmark",
    approach = list(benchmark = "benchmark"),
    spatial_scale = "global",
    only_first_filename = FALSE
  )

  # Get filenames for scenario and reference
  files_scenario <- get_filenames(
    config = config_scen,
    output_files = output_files
  )
  files_reference <- get_filenames(
    config = config_ref,
    output_files = output_files
  )

  # read in terrestrial area (in m2)
  terr_area <- lpjmlkit::read_io(files_scenario$terr_area)$data %>% drop()
  ncell <- length(terr_area)

  ####### PB Land-system change ################################################

  # read in cftfrac based on lpjml input to include areas with failed harvest
  cftfrac %<-% aggregate_lpjml_output(
    files_scenario$landuse,
    time_span_scenario,
    aggregate = list(year = mean)
  )

  #-------------- total agricultural area --------------------------------------

  # calculate global area for each band, conversion from m2 to mio ha
  area_cft <- apply(cftfrac * terr_area * 10^-10, "band", sum)

  # global cropland area
  indices_crops <- grep("grass|tree", names(area_cft),
    ignore.case = TRUE, invert = TRUE
  )
  crop_area <- sum(area_cft[indices_crops])

  # global pasture area
  indices_pastures <- grep("grassland", names(area_cft), ignore.case = TRUE)
  pasture_area <- sum(area_cft[indices_pastures])

  # global irrigated area
  indices_irrig <- grep("irrig", names(area_cft), ignore.case = TRUE)
  irrig_area <- sum(area_cft[indices_irrig])

  #------------- potential forest extent ---------------------------------------
  biome_classes <- classify_biomes(
    config_reference = config_reference,
    time_span_reference = time_span_reference
  )

  # read in biome mapping
  biome_mapping <- read_biome_mapping(
    system.file("extdata", "biomes.csv", package = "boundaries")
  )

  is_forest <- rep(FALSE, ncell)
  is_forest[
    biome_classes$biome_id %in% dplyr::filter(
      biome_mapping, `category` == "forest"
    )$id
  ] <- TRUE  # ids for tree plant funtional types

  forest_area <- sum(terr_area[is_forest]) * 10^-10

  ####### PB freshwater change #################################################

  #------------- Irrigation withdrawals and consumption ------------------------

  irrig %<-% aggregate_lpjml_output(
    files_scenario$irrig,
    time_span_scenario,
    aggregate = list(
      year = mean,
      month = sum
    )
  )

  conv_loss_evap %<-% aggregate_lpjml_output(
    files_scenario$conv_loss_evap,
    time_span_scenario,
    aggregate = list(
      year = mean,
      month = sum
    )
  )

  conv_loss_drain %<-% aggregate_lpjml_output(
    files_scenario$conv_loss_drain,
    time_span_scenario,
    aggregate = list(
      year = mean,
      month = sum
    )
  )

  return_flow_b %<-% aggregate_lpjml_output(
    files_scenario$return_flow_b,
    time_span_scenario,
    aggregate = list(
      year = mean,
      month = sum
    )
  )
  # withdrawals in km3
  wd <- global_sum((irrig + conv_loss_drain + conv_loss_evap), area = terr_area)

  # consumption in km3
  cons <- global_sum((irrig + conv_loss_evap - return_flow_b), area = terr_area)

  ####### PB Nitrogen ##########################################################
  #------------- Nitrogen balance and leaching ---------------------------------

  # # N leaching (total) in Tg
  leaching %<-% aggregate_lpjml_output(
    files_scenario$leaching,
    time_span_scenario,
    aggregate = list(
      year = mean
    )
  ) %>% global_sum(area = terr_area)

  # N leaching from cropland
  leaching_agr %<-% aggregate_lpjml_output(
    files_scenario$nleaching_agr,
    time_span_scenario,
    aggregate = list(
      year = mean
    )
  ) %>% global_sum(area = terr_area)

  # N use efficiency on cropland
  # inputs
  total_fert %<-% aggregate_lpjml_output(
    files_scenario$nfert_agr,
    time_span_scenario,
    aggregate = list(
      year = mean
    )
  ) %>% global_sum(area = terr_area)

  total_man %<-% aggregate_lpjml_output(
    files_scenario$nmanure_agr,
    time_span_scenario,
    aggregate = list(
      year = mean
    )
  ) %>% global_sum(area = terr_area)

  # For total_dep
  total_dep %<-% aggregate_lpjml_output(files_scenario$ndepo_agr,
    time_span_scenario,
    aggregate = list(
      year = mean
    )
  ) %>% global_sum(area = terr_area)


  # For total_bnf
  total_bnf %<-% aggregate_lpjml_output(files_scenario$bnf_agr,
    time_span_scenario,
    aggregate = list(
      year = mean
    )
  ) %>% global_sum(area = terr_area)

  total_seed %<-% aggregate_lpjml_output(files_scenario$seedn_agr,
    time_span_scenario,
    aggregate = list(
      year = mean
    )
  ) %>% global_sum(area = terr_area)


  n_inputs <- total_fert + total_man + total_dep + total_bnf + total_seed


  nharvest %<-% aggregate_lpjml_output(
    files_scenario$pft_harvestn,
    time_span_scenario,
    aggregate = list(
      year = mean
    )
  )
  if (gridbased == FALSE) {
    # based on lpjml output
    cftfrac %<-% aggregate_lpjml_output(
      files_scenario$cftfrac,
        time_span_scenario,
        aggregate = list(year = mean)
    )
    nharvest <- nharvest * cftfrac
  }

  # only refer to cropland
  indices_crops <- grep("grass|tree", names(nharvest),
    ignore.case = TRUE, invert = TRUE
  )
  nharvest_crop <- nharvest[, indices_crops]

  nharvest <- nharvest_crop %>% global_sum(area = terr_area)

  nue <- nharvest / n_inputs * 100

  # N surplus on cropland
  nsurplus <- n_inputs - nharvest

  ####### PB biosphere integrity ###############################################
  #------------- NPP -------------------------------------------

  npp_lu %<-% aggregate_lpjml_output(files_scenario$npp,
    time_span_scenario,
    aggregate = list(
      year = mean
    )
  ) %>% global_sum(area = terr_area) * 10^-3 # in PgC


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

  #------------- simulated lpjml harvest ---------------------------------------

  yield %<-% aggregate_lpjml_output(
    files_scenario$pft_harvestc,
    time_span_scenario,
    aggregate = list(year = mean)
  )

  # dry matter harvest g
  if (gridbased == FALSE) {
    # based on lpjml output
    harvest_dm <- cftfrac * yield * terr_area * (1 / dm_factor)
  } else {
    harvest_dm <- yield * terr_area * (1 / dm_factor)
  }

  # freshmatter harvest
  crop_fm_prod <- array(0, dim = c(ncell, 12)) # fm production for all CFTs

  # include "rainfed" and exclude "grass", "tree", and "other" in crop indices
  crop_ids_explicit_rainfed <- grep("^(?!.*(grass|tree|other)).*rainfed.*$",
    names(area_cft),
    ignore.case = TRUE, perl = TRUE
  )
  # #include "irig" and exclude "grass", "tree", and "other" in crop indices
  crop_ids_explicit_irrig <- grep("^(?!.*(grass|tree|other)).*irrig.*$",
    names(area_cft),
    ignore.case = TRUE, perl = TRUE
  )

  for (cc in seq_along(crop_ids_explicit_rainfed)) {
    crop_fm_prod[, cc] <- (
      harvest_dm[, crop_ids_explicit_rainfed[cc]] +
        harvest_dm[, crop_ids_explicit_irrig[cc]]
    ) * dm2fm[c] / 1000 # kg FM
  }

  crop_sum <- sum(crop_fm_prod) / 10^9 # in Mt, 1 Mt = 10^9 kg

  ##### combine simulated values with literature ###############################

  # lists to match literature character strings with the here computed variables
  var_match <- list(
    "crop", "irrig", "pasture", "forest", "withdr",
    "cons", "leaching", "leaching from crop",
    "efficiency", "surplus", "npp", "prod"
  )
  units_match <- list(
    "Mha", "Mha", "Mha", "Mha", "km3", "km3", "Tg N yr-1",
    "Tg N yr-1", "%", "Tg N yr-1", "PgC", "Mt"
  )
  calc_var <- list(
    crop_area, irrig_area, pasture_area, forest_area, wd, cons,
    leaching, leaching_agr, nue, nsurplus, npp_lu, crop_sum
  )

  # load literature table
  ref_table <- system.file("extdata", "global_validation_data.csv",
    package = "boundaries"
  ) %>%
    utils::read.csv2()

  # insert lpjml values in table (pb status variables later)
  ref_table <- insert_lpjml_value(ref_table, var_match, units_match, calc_var)


  # calculate PB status
  pb_status <- calc_status(
    boundary = c("lsc", "greenwater", "bluewater", "biosphere"),
    config_scenario = config_scenario,
    config_reference = config_reference,
    time_span_scenario = time_span_scenario,
    time_span_reference = time_span_reference,
    spatial_scale = "global",
    gridbased = gridbased,
    path_baseline = path_baseline,
    ...
  )

  # list of pb variables to insert in the table
  pb_list <- list(
    "lsc", "greenwater", "bluewater", "nitrogen",
    "biosphere"
  )

  # insert values from calc_status in the table
  ref_table <- evaluate_pb_status(ref_table, pb_list, pb_status)

  # add column for normalized error values
  ref_table <- ref_table %>%
    dplyr::mutate(norm_error = signif((lpjml_value - as.numeric(value)) /
                                        as.numeric(value), digits = 2))

  # create summary table
  summary_tbl <- ref_table %>%
    dplyr::group_by(variable) %>%
    dplyr::mutate(
      value_range = paste(range.lower, "-", range.upper),
      time_range = ifelse(!is.na(X.lower.),
        ifelse(X.lower. == "various", "various",
          paste(X.lower., "-", X.upper.)
        ), ""
      )
    ) %>%
    dplyr::mutate(
      assessment_time = ifelse(
        !is.na(year) & year != "",
        as.character(year),
        as.character(time_range)
      ),
      literature.range = paste(
        ifelse(
          !is.na(value) & value != "",
          as.character(value),
          as.character(value_range)
        ),
        paste(
          "[", as.character(assessment_time), "] (",
          Author.s., year.of.publication, ")"
        )
      )
    ) %>%
    dplyr::select(boundary, variable, lpjml_value, literature.range, unit) %>%
    dplyr::mutate(
      literature.range = paste0(literature.range, collapse = "; ")
    ) %>%
    dplyr::distinct()


  utils::write.csv(summary_tbl, file = table_path, row.names = FALSE)

  return(ref_table)
}


# define set of helper functions
# functions to deal with time and spacial scales of lpjml output

aggregate_lpjml_output <- function(
    file,
    timespan,
    aggregate = list()) {
  # file <- lpjmlkit::read_io(
  variable <- lpjmlkit::read_io(
    file,
    subset = list(year = timespan)
  ) %>%
    lpjmlkit::transform(to = c("year_month_day")) %>%
    lpjmlkit::as_array(aggregate = aggregate) # %>%
  return(variable)
}


# conversion from g/m2 to Tg
global_sum <- function(output, area) {
  sum_result <- sum(output * area) * 10^-12
  return(sum_result)
}

# insert computed variables in literature table through matching by variable and unit patterns,
# make sure structure of original .csv table is mantained

insert_lpjml_value <- function(
    ref_table, var_patterns, unit_patterns,
    calculated_vars, set_nan_if_empty = TRUE) {
  lpjml_values <- rep(NA_real_, nrow(ref_table))
  for (i in seq_along(var_patterns)) {
    matching_indices <- grepl(var_patterns[[i]], ref_table$variable, ignore.case = TRUE) &
      grepl(unit_patterns[[i]], ref_table$unit, ignore.case = TRUE)
    lpjml_values[matching_indices] <- ifelse(
      set_nan_if_empty && (calculated_vars[[i]] == 0 || is.nan(calculated_vars[[i]])),
      NaN, # set to NaN if set_nan_if_empty is TRUE and the calculated value is empty or NaN
      round(eval(parse(text = calculated_vars[[i]]))[1], digits = 2)
    )
  }
  ref_table$lpjml_value <- lpjml_values
  return(ref_table)
}


# insert the PB status values (calculated with calc_status) in the benchmarking table
evaluate_pb_status <- function(ref_table, pb_list, temp2) {
  # Please R CMD check for use of dplyr syntax
  lpjml_value <- NULL
  for (pb in pb_list) { # add calculation of PBs
    ref_table <- ref_table %>%
      dplyr::mutate(lpjml_value = ifelse(
        grepl(pb, ref_table$variable, ignore.case = TRUE) &
          grepl("PB status", ref_table$variable, ignore.case = TRUE),
        round(temp2[[pb]][1], digits = 2),
        lpjml_value
      ))
  }
  return(ref_table)
}
