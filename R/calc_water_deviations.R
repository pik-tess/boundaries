#' Calculate water status based on deviations of a monthly scenario variable
#'  from a corresponding monthly reference variable
#'
#' Calculate deviations (<q5 / >q95) for a monhtly variable in a scenario LPJmL
#' run as compared to a reference LPJmL run, wither referring to global area
#' share with deviations (spatial_resolution: global), or to number of years
#' with deviations (spatial resolution: cell). From this, calculate a global or
#' gridded PB status
#'
#' @param file_scenario character string with path to monthly variable of the
#'         scenario LPJmL run.
#'
#' @param file_reference character string with path to monthly variable of the
#'         reference LPJmL run.
#'
#' @param grid_path character string with path to a grid file
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#'        as a character string, e.g. `as.character(1982:2011)` (default)
#'
#' @param time_span_reference time span to be used for the reference run, defined
#'        defined as a character string (e.g. `as.character(1901:1930)`).
#'        Can differ in offset and length from `time_span_scenario`!
#'        If `NULL` value of `time_span_scenario` is used
#'
#' @param method method (character string) to be used , currently available
#'        method is `c("wang-erlandsson2022")` based on
#'        [Wang-Erlandsson et al. 2022](https://doi.org/10.1038/s43017-022-00287-8)
#'        (referring only to the driest/wettest month of each year) or
#'        `porkka_2023` based on
#'        [Porkka et al. 2023](https://eartharxiv.org/repository/view/3438/)
#'        (referring to each month of a year)
#' @param spatial_resolution character string indicating spatial resolution
#'        either "grid" for calculation of number of years with transgression
#'        (for wang-erlandsson2022: dim(ncell, nyears);
#'         for porkka_2023: dim(ncell, nyears, months)) or
#'        "global" for calculation of the share (%) of total global area with
#'        deviations (either one value per year (wang-erlandsson2022) or one
#'        value per year and month (porkka_2023))
#'
#' @param avg_nyear_args list of arguments to be passed to
#'        \link[pbfunctions]{average_nyear_window} (see for more info).
#'        To be used for time series analysis
#'
#' @examples
#' \dontrun{
#'  calc_deviations(file_scenario, file_reference, grid_path,
#'                 time_span_reference, spatial_resolution)
#' }
#'
#' @md
#' @export



calc_water_status <- function(file_scenario,
                              file_reference,
                              grid_path,
                              time_span_scenario = as.character(1982:2011),
                              time_span_reference = NULL,
                              method = "wang-erlandsson2022",
                              avg_nyear_args = list(),
                              spatial_resolution
                                   ) {

  # -------------------------------------------------------------------------- #
  # calc deviations, for global resolution referring to area with deviations;
  # for grid resolution referring to number of years/months with deviations
  deviations <- calc_deviations(
    file_scenario = file_scenario,
    file_reference = file_reference,
    grid_path = grid_path,
    time_span_scenario = time_span_scenario,
    time_span_reference =  time_span_reference,
    method = method,
    avg_nyear_args = avg_nyear_args,
    spatial_resolution = spatial_resolution
  )

  # -------------------------------------------------------------------------- #
  if (spatial_resolution == "grid") {
    #prop.test to test for significance in departure increases

    p_dry <- p_wet <- array(NA, 67420) #TODO make flexible
    if (method == "wang-erlandsson2022")
      trials <- c(length(time_span_reference), length(time_span_scenario))
    else if (method == "porkka_2023") {
      trials <- c(length(time_span_reference) * 12,
                   length(time_span_scenario) * 12)
    }
    for (i in seq_len(67420)) { #TODO make flexible
      test_dry <- prop.test(x = c(deviations$reference$dry[i],
                                  deviations$scenario$dry[i]),
                              n = trials,
                              alternative = "less") #TODO verify!
      p_dry[i] <- test_dry$p.value

      test_wet <- prop.test(x = c(deviations$reference$wet[i],
                                  deviations$scenario$wet[i]),
                              n = trials,
                              alternative = "less")
      p_wet[i] <- test_wet$p.value
    }
    #TODO understand where NA values come from
    p_wet[is.na(p_wet)] <- 1
    p_dry[is.na(p_dry)] <- 1
    #TODO translation into PB status only prelimary. Thresholds should be
    # definable in function parameters
    pb_status <- array(0,
                       dim = dim(p_dry))
    pb_status[p_wet < 0.05 | p_dry < 0.05] <- 2
    pb_status[p_wet < 0.05 & p_dry < 0.05] <- 3
    pb_status[p_wet >= 0.05 & p_dry >= 0.05] <- 1

  # -------------------------------------------------------------------------- #
  } else if (spatial_resolution == "global") {
    #calculate q95 for area with deviations in the reference
    q95_area <- quantile(deviations$reference$wet_or_dry,
                         probs = 0.95, na.rm = T)
    q50_area <- quantile(deviations$reference$wet_or_dry,
                         probs = 0.5, na.rm = T)
    control_variable <- mean(deviations$scenario$wet_or_dry)

    attr(control_variable,"thresholds") <- c(holocene = q50_area, pb = q95_area) # h: 0, pb: 1

  }
  return(control_variable)
}

# todo: finish
as_risk_level <- function(control_variable, type = "continuous"){
  thresholds <- attr(control_variable,"thresholds")
  if (type == "continuous"){
    # pb value normalized to 0-1, 1 is the
    # threshold between safe and increasing risk, >1 is transgressed
    pb_status$continuous_normalized <- (pb_status$control_variable -
                                          thresholds[["hol"]]) / (thresholds[["pb"]] - thresholds[["hol"]])
    pb_status$risk_level <- ifelse(pb_status$control_variable > q95_area,
                                   2, 1) # only safe and transgressed
  }else{ # type == "discrete"

  }
  return(pb_status)
}

calc_deviations <- function(file_scenario,
                            file_reference,
                            grid_path,
                            time_span_scenario = as.character(1982:2011),
                            time_span_reference = NULL,
                            method = "porkka_2023",
                            avg_nyear_args = list(),
                            spatial_resolution
                                   ) {

# reference
  var_reference <- lpjmlkit::read_io(
      file_reference, subset = list(year = time_span_reference)
      ) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      lpjmlkit::as_array(aggregate = list(band = sum)) %>%
      suppressWarnings()

  # scenario
  var_scenario <- lpjmlkit::read_io(
      file_scenario, subset = list(year = time_span_scenario)
      ) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      lpjmlkit::as_array(aggregate = list(band = sum)) %>%
      suppressWarnings()

  # -------------------------------------------------------------------------- #
  #calculate the 5% and 95% quantiles of the baseline period
  quants <- calc_water_baseline(var_reference,
                                method = method)

  # calculate number of months/years with dry & wet departures (grid resolution)
  # or area with dry/wet departures (global resolution)

  ref_depart <- calc_water_depart(var_reference, quants,
                                    spatial_resolution = spatial_resolution,
                                    method = method)
  scen_depart <- calc_water_depart(var_scenario, quants,
                                     spatial_resolution = spatial_resolution,
                                     method = method)

  return(list(reference = ref_depart, scenario = scen_depart))
}



# quantile functions
q5  <- function(x) quantile(x, probs = 0.05, na.rm = T)
q95 <- function(x) quantile(x, probs = 0.95, na.rm = T)

# calculate the baseline quantiles
calc_water_baseline <- function(file_reference, method) {
  if (method == "wang-erlandsson2022") {
    dry_base_yr <- apply(file_reference, c(1, 3), min)
    wet_base_yr <- apply(file_reference, c(1, 3), max)
  } else if (method == "porkka_2023") {
    dry_base_yr <- wet_base_yr <- file_reference
  }
  # calc 5% and 95% percentile for each cell
  # for each year over all months over baseline period
  dim_remain <- names(dim(dry_base_yr))[names(dim(dry_base_yr)) != "year"]
  q5_base  <- apply(dry_base_yr, dim_remain, q5)
  q95_base <- apply(wet_base_yr, dim_remain, q95)
  quants <- list(q5 = q5_base, q95 = q95_base)
  return(quants)
}

# calculate GW dry & wet departures and return mean annual area of departure
# (global resolution) or number of years/months with wet/dry departures
# (grid resolution)

calc_water_depart <- function(file_scenario, quants, spatial_resolution,
                              method) {

  q5_base <- rep(quants[["q5"]], dim(file_scenario)["year"])
  q95_base <- rep(quants[["q95"]], dim(file_scenario)["year"])

  if (method == "wang-erlandsson2022") {
    # driest/ wettest month per gridcell for each year -> ignores which month
    dry <- apply(file_scenario, c("cell", "year"), min)
    wet <- apply(file_scenario, c("cell", "year"), max)
  } else if (method == "porkka_2023") {
    dry <- wet <- file_scenario
  }

  # identify cells with dry/wet departures
  dry[dry >= q5_base] <- NA
  dry[dry >= 0] <- 1
  wet[wet <= q95_base] <- NA
  wet[wet >= 0] <- 1
  wet[is.na(wet)] <- 0
  dry[is.na(dry)] <- 0

  result <- list()
  if (spatial_resolution == "grid") {
    #dim_remain <- names(dim(dry))[names(dim(dry)) != "year"]
    result$dry <- apply(dry, "cell", sum)
    result$wet <- apply(wet, "cell", sum)

  } else if (spatial_resolution == "global") {
    grid <- lpjmlkit::read_io(
      grid_path,
      silent = TRUE
      )
    # TODO this should rather be the percentage of ice-free land surface!
    cell_area <- lpjmlkit::calc_cellarea(grid)
    dim_remain <- names(dim(dry))[names(dim(dry)) != "cell"]
    result$dry <- apply((dry * cell_area) / sum(cell_area) * 100, dim_remain,
                           sum)
    result$wet <- apply((wet * cell_area) / sum(cell_area) * 100, dim_remain,
                           sum)
    wet_or_dry <- dry + wet
    wet_or_dry[wet_or_dry > 1] <- 1
    result$wet_or_dry <- apply((wet_or_dry * cell_area) / sum(cell_area) * 100,
                                dim_remain, sum)
  }
  return(result)
}
