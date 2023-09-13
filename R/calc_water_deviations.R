#' Calculate water status based on deviations of a monthly scenario variable
#'  from a corresponding monthly reference variable
#'
#' Calculate deviations (<q5 / >q95) for a monhtly variable in a scenario LPJmL
#' run as compared to a reference LPJmL run, either referring to global area
#' share with deviations (spatial_resolution: global), or to number of months or
#' years with deviations (spatial resolution: cell). From this, calculate a
#' global or gridded PB status
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
#' @param time_span_reference time span to be used for the reference run,
#'        defined as a character string (e.g. `as.character(1901:1930)`).
#'        Can differ in offset and length from `time_span_scenario`!
#'        If `NULL` value of `time_span_scenario` is used
#'
#' @param method method (character string) to be used , currently available
#'        method is `c("wang-erlandsson2022")` based on
#'        [Wang-Erlandsson et al. 2022](https://doi.org/10.1038/s43017-022-00287-8)
#'        (referring only to the driest/wettest month of each year) or
#'        `porkka2023` based on
#'        [Porkka et al. 2023](https://eartharxiv.org/repository/view/3438/)
#'        (referring to each month of a year; default)
#'
#' @param thresholds list with thresholds to be used to
#'        define the safe, increasing risk and high risk zone,
#'        e.g. c(holocene = 0.5, pb = 0.95, highrisk = 0.99).
#'        For spatial resolution = "grid", this refers to the p value
#'        (significance level of increases in deviations) with the default:
#'        c(holocene = 1, pb = 0.05, highrisk = 0.01).
#'        For spatial resolution = "global", this refers to the quantiles of
#'        the global area with deviations in the reference period. The dafault
#'        for global resolution is: c(holocene = 0.5, pb = 0.95,
#'        highrisk = 0.99).
#'        If set to NULL, the respective default is taken (see above; matching
#'        the spatial_resolution)
#' 
#' @param spatial_resolution character string indicating spatial resolution
#'        either "grid" for calculation of number of years with transgression
#'        (for wang-erlandsson2022: dim(ncell, nyears);
#'         for porkka2023: dim(ncell, nyears, months)) or
#'        "global" for calculation of the share (%) of total global area with
#'        deviations (either one value per year (wang-erlandsson2022) or one
#'        value per year and month (porkka2023))
#'
#' @param avg_nyear_args list of arguments to be passed to
#'        \link[pbfunctions]{average_nyear_window} (see for more info).
#'        To be used for time series analysis
#'
#' @examples
#' \dontrun{
#'  calc_water_status(file_scenario, file_reference, grid_path,
#'                 time_span_reference, spatial_resolution)
#' }
#'
#' @md
#' @export
calc_water_status <- function(file_scenario,
                              file_reference,
                              grid_path,
                              time_span_scenario = as.character(1982:2011),
                              time_span_reference,
                              method = "porkka2023",
                              thresholds = NULL,
                              avg_nyear_args = list(),
                              spatial_resolution = "grid") {

  # read in reference and scenario output
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
  # calculate the 5% and 95% quantiles of the baseline period for each cell
  # either for each month of the year or only for driest/wettest month
  # depending on the defined method
  quants <- calc_baseline(var_reference,
                                method = method)

  # -------------------------------------------------------------------------- #
  # calculate number of months/years with dry & wet departures (grid resolution)
  # or area with dry/wet departures (global resolution) for each cell
  ref_depart <- calc_departures(var_reference, grid_path, quants,
                                    spatial_resolution = spatial_resolution,
                                    method = method)
  scen_depart <- calc_departures(var_scenario, grid_path, quants,
                                     spatial_resolution = spatial_resolution,
                                     method = method)


  # -------------------------------------------------------------------------- #
  if (spatial_resolution == "grid") {

    # define the number of trials (= number of years or months within the
    # reference and scenario timeframe)
    if (method == "wang-erlandsson2022") {
      trials <- c(length(time_span_reference), length(time_span_scenario))
      # TODO adapt for timeseries analysis
    } else if (method == "porkka2023") {
      # calculate number of months within the time_span
      trials <- c(length(time_span_reference) * 12,
                   length(time_span_scenario) * 12)
    }
    # test for significance in departure increases,
    # based on vectorized prop.test

    # calculate the proportions
    ref_prop <- ref_depart$wet_or_dry / trials[1]
    scen_prop <- scen_depart$wet_or_dry / trials[2]

    # Calculate the standard errors
    reference_se <- sqrt(ref_prop * (1 - ref_prop) / trials[1])
    new_se <- sqrt(scen_prop * (1 - scen_prop) / trials[2])

    # Calculate the z-scores,  taking into account the variability of the
    # proportions and the sample sizes
    z_scores <- (
      scen_prop - ref_prop - 0.5 * (1 / trials[2] - 1 / trials[1])) /
      sqrt(reference_se^2 + new_se^2 + (1 / trials[1] + 1 / trials[2]) * 0.25
    )
    # Calculate the p-value
    control_variable <- pnorm(z_scores, lower.tail = TRUE)

    # # old, non-vectorized calculation based on the prop.test
    # # test for significance in departure increases based on prop.test
    # p_wet_or_dry <- array(NA, 67420) #TODO make flexible
    # # TODO test if possible to hand over table/matrix for all cells instead
    # # of single cells
    # for (i in seq_len(67420)) { #TODO make flexible
    #   test_wet_or_dry <- prop.test(x = c(ref_depart$wet_or_dry[i],
    #                               scen_depart$wet_or_dry[i]),
    #                               n = trials,
    #                               alternative = "less") #TODO verify!
    #   p_wet_or_dry[i] <- test_wet_or_dry$p.value
    # }
    # # remove NA values in polar regions
    # p_wet_or_dry[is.na(p_wet_or_dry)] <- 1
    # # 1 - x to be compatible with other boundaries (holocene < pb < highrisk):
    # control_variable <- 1 - p_wet_or_dry

    # thresholds for translation into pb status, set to default if thresholds
    # are not explicitely defined
    if (is.null(thresholds)) {
      thresholds <- c(holocene = 1,
                        pb = 0.05,
                        highrisk = 0.01)
    }

    #TODO translation into PB status only prelimary.
    attr(control_variable, "thresholds") <-  1 - thresholds

  # -------------------------------------------------------------------------- #
  } else if (spatial_resolution == "global") {


    # calculate areas corresponding to the quantiles defined in thresholds
    area_high_risk <- quantile(ref_depart$wet_or_dry,
                         probs = thresholds[["highrisk"]], na.rm = TRUE)
    area_pb <- quantile(ref_depart$wet_or_dry,
                         probs = thresholds[["pb"]], na.rm = TRUE)
    area_holocene <- quantile(ref_depart$wet_or_dry,
                         probs = thresholds[["holocene"]], na.rm = TRUE)

    # transform to array for compatability with average_nyear_window
    # TODO better implement option in average_nyear_window to work with one
    # dimensional arrays (hard coded cell = 1)
    # alternative: apply the moving average funtion already to var_scenario?
    scen_depart_array <- array(scen_depart$wet_or_dry,
                               dim = c(cell = 1,
                                       year = length(scen_depart$wet_or_dry)
                                       ),
                               dimnames = list(cell = 1,
                                               year = names(scen_depart$wet_or_dry)
                                               )) #dummy cell dimension

    control_variable <- do.call(average_nyear_window,
                                      append(list(x = scen_depart_array),
                                             avg_nyear_args))[1, ]

    attr(control_variable, "thresholds") <- list(holocene = area_holocene,
                                              pb = area_pb,
                                              highrisk = area_high_risk)
    attr(control_variable, "control variable") <- "area with wet/dry departures (%)"

  }
  return(control_variable)
}



# quantile functions
q5  <- function(x) quantile(x, probs = 0.05, na.rm = TRUE)
q95 <- function(x) quantile(x, probs = 0.95, na.rm = TRUE)
#TODO make q5/95 flexibel = parameters?

# calculate the baseline quantiles
calc_baseline <- function(file_reference, method) {
  if (method == "wang-erlandsson2022") {
    dry_base_yr <- apply(file_reference, c(1, 3), min)
    wet_base_yr <- apply(file_reference, c(1, 3), max)
  } else if (method == "porkka2023") {
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

calc_departures <- function(data, grid_path, quants,
                              spatial_resolution, method) {

  q5_base <- rep(quants[["q5"]], dim(data)["year"])
  q95_base <- rep(quants[["q95"]], dim(data)["year"])

  if (method == "wang-erlandsson2022") {
    # driest/ wettest month per gridcell for each year -> ignores which month
    dry <- apply(data, c("cell", "year"), min)
    wet <- apply(data, c("cell", "year"), max)
  } else if (method == "porkka2023") {
    dry <- wet <- dry_or_wet <- data
  }

  # identify cells with dry/wet departures
  dry[dry >= q5_base] <- NA
  dry[dry >= 0] <- 1
  wet[wet <= q95_base] <- NA
  wet[wet >= 0] <- 1
  wet[is.na(wet)] <- 0
  dry[is.na(dry)] <- 0
  dry_or_wet <- ifelse((dry == 1 | wet == 1), 1, 0)

  result <- list()
  if (spatial_resolution == "grid") {
    #dim_remain <- names(dim(dry))[names(dim(dry)) != "year"]
    result$dry <- apply(dry, "cell", sum)
    result$wet <- apply(wet, "cell", sum)
    result$wet_or_dry <- apply(dry_or_wet, "cell", sum)
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
