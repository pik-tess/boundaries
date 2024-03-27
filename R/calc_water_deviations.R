#' Calculate water status based on deviations of a monthly scenario variable
#'  from a corresponding monthly reference variable
#'
#' Calculate deviations (<q5 / >q95) for a monhtly variable in a scenario LPJmL
#' run as compared to a reference LPJmL run, either referring to global area
#' share with deviations (spatial_scale: global), or to number of months or
#' years with deviations (spatial resolution: cell). From this, calculate a
#' global or gridded PB status
#'
#' @param files_scenario list with variable names and corresponding file paths
#' (character string) of the scenario LPJmL run. All needed files are
#' provided in XXX. E.g.: list(leaching = "/temp/leaching.bin.json")
#'
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the reference LPJmL run. All needed files are
#' provided in XXX. E.g.: list(leaching = "/temp/leaching.bin.json"). If not
#' needed for the applied approach, set to NULL.
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#' as character string, e.g. `as.character(1982:2011)` (default)
#'
#' @param time_span_reference time span to be used for the reference run,
#' defined as a character string (e.g. `as.character(1901:1930)`).
#' Can differ in offset and length from `time_span_scenario`!
#' If `NULL` value of `time_span_scenario` is used
#'
#' @param approach approach (character string) to be used , currently available
#' approach is `c("wang-erlandsson2022")` based on
#' [Wang-Erlandsson et al. 2022](https://doi.org/10.1038/s43017-022-00287-8)
#' (referring only to the driest/wettest month of each year) or
#' `porkka2023` based on
#' [Porkka et al. 2023](https://eartharxiv.org/repository/view/3438/)
#' (referring to each month of a year; default)
#'
#' @param thresholds list with thresholds to be used to
#' define the safe, increasing risk and high risk zone,
#' For spatial_scale = "global" and "subglobal", this refers to the quantiles of
#' the global/basin area with deviations in the reference period. The default
#' is: c(holocene = 0.5, pb = 0.95, highrisk = NULL).
#' If set to NULL, the  default is taken from metric_files.yml
#' For highrisk, the value is currently hard-coded to 0.5
#' (following Richardson et al. 2023)
#'
#' @param spatial_scale character string indicating spatial scale;
#' "global" or "subglobal" for calculation of the share (%) of total
#' global/basin area with deviations (either one value per year
#' (wang-erlandsson2022) or one value per year and month (porkka2023)); "grid"
#' not yet defined
#'
#' @param time_aggregation_args list of arguments to be passed to
#' [`aggregate_time`] (see for more info).
#' To be used for time series analysis
#'
#' @param config_args list of arguments to be passed on from the model
#' configuration.
#'
#' @param variable character string with the name of the variable to be used
#' for the calculation of the water deviations. Default is "rootmoist"
#'
#' @examples
#' \dontrun{
#'  calc_water_deviations(file_scenario, file_reference, terr_area_path,
#'                 time_span_reference, spatial_scale)
#' }
#'
#' @md
#' @export
calc_water_deviations <- function(files_scenario,
                                  files_reference,
                                  spatial_scale = "subglobal",
                                  time_span_scenario = NULL,
                                  time_span_reference,
                                  approach = "porkka2023",
                                  thresholds = NULL,
                                  time_aggregation_args = list(),
                                  config_args = list(),
                                  variable = "rootmoist") {

  # read in reference and scenario output
  # please R CMD check for use of future operator
  var_reference <- var_scenario <- NULL
  # reference
  var_reference %<-% read_io_format(
    files_reference[[variable]],
    time_span_reference,
    aggregate = list(band = sum),
    spatial_subset = config_args$spatial_subset
  )

  # scenario
  var_scenario %<-% read_io_format(
    files_scenario[[variable]],
    time_span_scenario,
    aggregate = list(band = sum),
    spatial_subset = config_args$spatial_subset
  )

  # -------------------------------------------------------------------------- #
  # calculate the 5% and 95% quantiles of the baseline period for each cell
  # either for each month of the year or only for driest/wettest month
  # depending on the defined approach
  quants <- calc_baseline(var_reference,
                          approach = approach)

  # -------------------------------------------------------------------------- #
  # calculate the area with dry/wet departures (subglobal and global resolution)

  if (spatial_scale == "subglobal") {
    # TODO: needs commenting
    # calculate for each basin and year: area with wet/dry departures

    cellinfo <- indexing_drainage(drainage_file = files_scenario$drainage)
    if (!is.null(config_args$spatial_subset)) {
      cellinfo <- lpjmlkit::asub(cellinfo, cell = config_args$spatial_subset)
    }
    endcell <- lpjmlkit::asub(cellinfo, band = "endcell")
  } else {
    endcell <- NULL
  }

  # please R CMD check for use of future operator
  ref_depart <- scen_depart <- NULL
  ref_depart %<-% calc_departures(
    data = var_reference,
    files_path = files_reference,
    quants = quants,
    spatial_scale = spatial_scale,
    approach = approach,
    spatial_subset = config_args$spatial_subset,
    endcell = endcell
  )

  scen_depart %<-% calc_departures(
    data = var_scenario,
    files_path = files_scenario,
    quants = quants,
    spatial_scale = spatial_scale,
    approach = approach,
    spatial_subset = config_args$spatial_subset,
    endcell = endcell
  )

  # please R CMD check for use of future operator
  area_high_risk <- area_pb <- area_holocene <- NULL
  # -------------------------------------------------------------------------- #
  if (spatial_scale == "global") {

    # calculate areas corresponding to the quantiles defined in thresholds
    if (is.null(thresholds[["highrisk"]])) {
      area_high_risk <- 50 # hard-coded based on Richardson et al. 2023 (SI)
    } else {
      area_high_risk %<-% stats::quantile(
        ref_depart$wet_or_dry,
        probs = thresholds[["highrisk"]],
        na.rm = TRUE
      )
    }

    area_pb %<-% stats::quantile(
      ref_depart$wet_or_dry,
      probs = thresholds[["pb"]],
      na.rm = TRUE
    )
    area_holocene %<-% stats::quantile(
      ref_depart$wet_or_dry,
      probs = thresholds[["holocene"]],
      na.rm = TRUE
    )

    if (approach == "porkka2023") {
      scen_depart$wet_or_dry <- apply(
        scen_depart$wet_or_dry,
        "year",
        mean,
        na.rm = TRUE
      )
    }
    control_variable <- do.call(
      aggregate_time,
      append(
        list(x = scen_depart$wet_or_dry),
        time_aggregation_args
      )
    )

    attr(control_variable, "thresholds") <- list(
      holocene = area_holocene,
      pb = area_pb,
      highrisk = area_high_risk
    )

  } else if (spatial_scale == "subglobal") {

    # calculate for each basin: thresholds for translation into pb status
    area_high_risk %<-% apply(
      ref_depart$wet_or_dry,
      "basin",
      function(x) {
        y <- stats::quantile(x, probs = thresholds[["highrisk"]], na.rm = TRUE)
        y
      }
    )

    area_pb %<-% apply(
      ref_depart$wet_or_dry,
      "basin",
      function(x) {
        y <- stats::quantile(x, probs = thresholds[["pb"]], na.rm = TRUE)
        y
      }
    )

    area_holocene %<-% apply(
      ref_depart$wet_or_dry,
      "basin",
      function(x) {
        y <- stats::quantile(x, probs = thresholds[["holocene"]], na.rm = TRUE)
        y
      }
    )

    # transform to array with cell dim as in lpjml
    # for vectors:
    basin_to_cell <- function(data, endcell) {
      ncells <- length(endcell)
      basins <- unique(endcell)
      out <- array(NA, dim = c(cell = ncells),
                   dimnames = list(cell = 1:ncells))
      for (i in seq_along(basins)){
        basin <- basins[i]
        out[which(endcell == basin)] <- data[i]
      }
      return(out)
    }

    area_high_risk %<-% basin_to_cell(area_high_risk, endcell)
    area_pb %<-% basin_to_cell(area_pb, endcell)
    area_holocene %<-% basin_to_cell(area_holocene, endcell)

    # for two dimensional array
    scen_departures <- array(
      NA,
      dim = c(
        cell = length(endcell),
        year = dim(scen_depart$wet_or_dry)[["year"]]
      ),
      dimnames = list(
        cell = seq_along(endcell),
        year = dimnames(scen_depart$wet_or_dry)[["year"]]
      )
    )

    basins <- unique(endcell)
    for (i in seq_along(basins)) {
      y <- which(endcell == basins[i])
      scen_departures[y, ] <- scen_depart$wet_or_dry[i, ]
    }

    control_variable <- do.call(
      aggregate_time,
      append(
        list(x = scen_departures),
        time_aggregation_args
      )
    ) %>%
      drop()

    attr(control_variable, "thresholds") <- list(
      holocene = area_holocene,
      pb = area_pb,
      highrisk = area_high_risk
    )
  }
  attr(control_variable, "control_variable") <-
    "area with wet/dry departures (%)"
  attr(control_variable, "spatial_scale") <- spatial_scale

  class(control_variable) <- c("control_variable")
  return(control_variable)
}




# calculate GW dry & wet departures and return mean annual area of departure
#   (global or basic scale)
calc_departures <- function(
  data,
  files_path,
  quants,
  spatial_scale,
  approach,
  spatial_subset,
  endcell
) {

  q5_base <- rep(quants[["q5"]], dim(data)["year"])
  q95_base <- rep(quants[["q95"]], dim(data)["year"])

  if (approach == "wang-erlandsson2022") {
    # driest/ wettest month per gridcell for each year -> ignores which month
    dry <- apply(data, c("cell", "year"), min)
    wet <- apply(data, c("cell", "year"), max)
  } else if (approach == "porkka2023") {
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

  # calc terrestrial area if spatial resolution is subglobal or global
  # TODO this should be better the percentage of ice-free land surface!
  if (spatial_scale == "global" | spatial_scale == "subglobal") { # nolint
    terr_area <- lpjmlkit::read_io(
      files_path$terr_area,
      silent = TRUE
    ) %>%
      lpjmlkit::transform("year_month_day") %>%
      lpjmlkit::as_array()
  }

  control_variable <- list()
  if (spatial_scale == "subglobal") {
    # TODO: needs commenting
    # calculate for each basin and year: area with wet/dry departures

    nbasins <- length(unique(endcell))
    control_variable$wet_or_dry <- array(
      NA,
      dim = c(basin = nbasins, dim(dry)["year"]),
      dimnames = list(basin = unique(endcell),
                      year = dimnames(dry)[["year"]])
    )

    for (b in unique(endcell)) { # go through all basins
      control_variable$wet_or_dry[which(unique(endcell) == b), ] <- (
        apply(
          (lpjmlkit::asub(dry_or_wet, cell = which(endcell == b)) *
             lpjmlkit::asub(terr_area, cell = which(endcell == b))) /
            sum(lpjmlkit::asub(terr_area, cell = which(endcell == b))) * 100,
          c("year"), sum, na.rm = TRUE
        )
      )
    }

  } else if (spatial_scale == "global") {
    # TODO: needs commenting
    terr_area <- drop(terr_area)

    dim_remain <- names(dim(dry))[names(dim(dry)) != "cell"]
    control_variable$dry <- apply(
      (dry * terr_area) / sum(terr_area) * 100,
      dim_remain,
      sum
    )
    control_variable$wet <- apply(
      (wet * terr_area) / sum(terr_area) * 100,
      dim_remain,
      sum
    )
    wet_or_dry <- dry + wet
    wet_or_dry[wet_or_dry > 1] <- 1
    control_variable$wet_or_dry <- apply(
      (wet_or_dry * terr_area) / sum(terr_area) * 100,
      dim_remain,
      sum
    )
  }

  class(control_variable) <- c("control_variable")
  return(control_variable)
}


# quantile functions
q5  <- function(x) stats::quantile(x, probs = 0.05, na.rm = TRUE)
q95 <- function(x) stats::quantile(x, probs = 0.95, na.rm = TRUE)


# calculate the baseline quantiles
calc_baseline <- function(file_reference, approach) {

  if (approach == "wang-erlandsson2022") {
    dry_base_yr <- apply(file_reference, c("cell", "year"), min)
    wet_base_yr <- apply(file_reference, c("cell", "year"), max)

  } else if (approach == "porkka2023") {
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
