#' Calculate water status based on deviations of a monthly scenario variable
#'  from a corresponding monthly reference variable
#'
#' Calculate deviations (<q5 / >q95) for a monhtly variable in a scenario LPJmL
#' run as compared to a reference LPJmL run, either referring to global area
#' share with deviations (spatial_scale: global), or to number of months or
#' years with deviations (spatial resolution: cell). From this, calculate a
#' global or gridded PB status
#'
#' @param file_scenario character string with path to monthly variable of the
#' scenario LPJmL run.
#'
#' @param file_reference character string with path to monthly variable of the
#' reference LPJmL run.
#'
#' @param terr_area_path character string with path to a grid file
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#' as a character string, e.g. `as.character(1982:2011)` (default)
#'
#' @param time_span_reference time span to be used for the reference run,
#' defined as a character string (e.g. `as.character(1901:1930)`).
#' Can differ in offset and length from `time_span_scenario`!
#' If `NULL` value of `time_span_scenario` is used
#'
#' @param method method (character string) to be used , currently available
#' method is `c("wang-erlandsson2022")` based on
#' [Wang-Erlandsson et al. 2022](https://doi.org/10.1038/s43017-022-00287-8)
#' (referring only to the driest/wettest month of each year) or
#' `porkka2023` based on
#' [Porkka et al. 2023](https://eartharxiv.org/repository/view/3438/)
#' (referring to each month of a year; default)
#'
#' @param thresholds list with thresholds to be used to
#' define the safe, increasing risk and high risk zone,
#' For spatial scale = "global" and "subglobal", this refers to the quantiles of
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
#' @param avg_nyear_args list of arguments to be passed to
#' \link[boundaries]{average_nyear_window} (see for more info).
#' To be used for time series analysis
#'
#' @examples
#' \dontrun{
#'  calc_water_status(file_scenario, file_reference, terr_area_path,
#'                 time_span_reference, spatial_scale)
#' }
#'
#' @md
#' @export
calc_water_status <- function(file_scenario,
                              file_reference,
                              terr_area_path,
                              drainage_path,
                              time_span_scenario = NULL,
                              time_span_reference,
                              method = "porkka2023",
                              thresholds = NULL,
                              avg_nyear_args = list(),
                              spatial_scale = "subglobal") {

  # read in reference and scenario output
  # reference
  var_reference %<-% read_io_format(
    file_reference,
    time_span_reference,
    aggregate = list(band = sum)
  )

  # scenario
  var_scenario %<-% read_io_format(
    file_scenario,
    time_span_scenario,
    aggregate = list(band = sum)
  )
 

  # -------------------------------------------------------------------------- #
  # calculate the 5% and 95% quantiles of the baseline period for each cell
  # either for each month of the year or only for driest/wettest month
  # depending on the defined method
  quants <- calc_baseline(var_reference,
                          method = method)

  # -------------------------------------------------------------------------- #
  # calculate the area with dry/wet departures (subglobal and global resolution)

  ref_depart %<-% calc_departures(var_reference, terr_area_path, drainage_path,
                                quants, spatial_scale = spatial_scale,
                                method = method)
  scen_depart %<-% calc_departures(var_scenario, terr_area_path, drainage_path,
                                 quants, spatial_scale = spatial_scale,
                                 method = method)

  # -------------------------------------------------------------------------- #
  if (spatial_scale == "global") {

    # calculate areas corresponding to the quantiles defined in thresholds
    if (is.null(thresholds[["highrisk"]])) {
      area_high_risk <- 50 # hard-coded based on Richardson et al. 2023 (SI)
    } else {
      area_high_risk %<-% quantile(ref_depart$wet_or_dry,
                               probs = thresholds[["highrisk"]], na.rm = TRUE)
    }
    
    area_pb %<-% quantile(ref_depart$wet_or_dry,
                        probs = thresholds[["pb"]], na.rm = TRUE)
    area_holocene %<-% quantile(ref_depart$wet_or_dry,
                              probs = thresholds[["holocene"]], na.rm = TRUE)

    if (method == "porkka2023") {
      scen_depart$wet_or_dry <- apply(scen_depart$wet_or_dry, "year",
                                      mean, na.rm = TRUE)
    }
    control_variable <- do.call(average_nyear_window,
                                append(list(x = scen_depart$wet_or_dry),
                                       avg_nyear_args))

    attr(control_variable, "thresholds") <- list(holocene = area_holocene,
                                                 pb = area_pb,
                                                 highrisk = area_high_risk)
    attr(control_variable, "control variable") <- (
      "area with wet/dry departures (%)"
    )

  } else if (spatial_scale == "subglobal") {
    # calculate for each basin: thresholds for translation into pb status
    area_high_risk %<-% apply(ref_depart$wet_or_dry, "basin",
                              function(x) {
                                y <- quantile(x, probs = thresholds[["highrisk"]], # nolint:line_length_linter
                                              na.rm = TRUE)
                              })
    area_pb %<-% apply(ref_depart$wet_or_dry, "basin",
                       function(x) {
                         y <- quantile(x, probs = thresholds[["pb"]],
                                       na.rm = TRUE)
                       })
    area_holocene %<-% apply(ref_depart$wet_or_dry, "basin",
                             function(x) {
                               y <- quantile(x, probs = thresholds[["holocene"]], # nolint:line_length_linter
                                             na.rm = TRUE)
                             })

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
    scen_departures <- array(NA,
                             dim = c(cell = length(endcell),
                                     year = dim(scen_depart$wet_or_dry)[["year"]]), #nolint
                             dimnames = list(cell = 1:length(endcell),
                                             year = dimnames(scen_depart$wet_or_dry)[["year"]])) #nolint

    basins <- unique(endcell)
    for (i in 1:length(basins)) {
      y <- which(endcell == basins[i])
      scen_departures[y, ] <- scen_depart$wet_or_dry[i, ]
    }
    control_variable <- do.call(average_nyear_window,
                                append(list(x = scen_departures),
                                       avg_nyear_args)) %>%
      drop()
    attr(control_variable, "thresholds") <- list(holocene = area_holocene,
                                                 pb = area_pb,
                                                 highrisk = area_high_risk)
    attr(control_variable, "control variable") <-
      "area with wet/dry departures (%)"

  }
  return(control_variable)
}



# quantile functions
q5  <- function(x) quantile(x, probs = 0.05, na.rm = TRUE)
q95 <- function(x) quantile(x, probs = 0.95, na.rm = TRUE)

# calculate the baseline quantiles
calc_baseline <- function(file_reference, method) {
  if (method == "wang-erlandsson2022") {
    dry_base_yr <- apply(file_reference, c("cell", "year"), min)
    wet_base_yr <- apply(file_reference, c("cell", "year"), max)
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
# (global or basic scale)

calc_departures <- function(data, terr_area_path, drainage_path, quants,
                            spatial_scale, method) {

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

  # calc terrestrial area if spatial resolution is subglobal or global
  # TODO this should be better the percentage of ice-free land surface!
  if (spatial_scale == "global" | spatial_scale == "subglobal") {
    terr_area <- lpjmlkit::read_io(
      terr_area_path,
      silent = TRUE
    )$data %>% drop()
  }

  result <- list()
  if (spatial_scale == "subglobal") {
    # calculate for each basin and year: area with wet/dry departures
    cellinfo <- indexing_drainage(drainage_file = drainage_path)
    endcell <<- lpjmlkit::asub(cellinfo, band = "endcell")
    nbasins <- length(unique(endcell))
    result$wet_or_dry <- array(NA, dim = c(basin = nbasins, dim(dry)["year"]),
                               dimnames = list(basin = unique(endcell),
                                               year = dimnames(dry)[["year"]]))
    for (b in unique(endcell)) { # go through all basins
      print(b)
      if (length(which(endcell == b) == TRUE) == 1) {
        result$wet_or_dry[which(unique(endcell) == b), ] <-
          (dry_or_wet[endcell == b, ] * terr_area[endcell == b]) /
          sum(terr_area[endcell == b]) * 100
      } else {
        result$wet_or_dry[which(unique(endcell) == b), ] <-
          apply((dry_or_wet[endcell == b, ] * terr_area[endcell == b]) /
                  sum(terr_area[endcell == b]) * 100,
                c("year"), sum, na.rm = TRUE)
      }
    }

  } else if (spatial_scale == "global") {
    dim_remain <- names(dim(dry))[names(dim(dry)) != "cell"]
    result$dry <- apply((dry * terr_area) / sum(terr_area) * 100, dim_remain,
                        sum)
    result$wet <- apply((wet * terr_area) / sum(terr_area) * 100, dim_remain,
                        sum)
    wet_or_dry <- dry + wet
    wet_or_dry[wet_or_dry > 1] <- 1
    result$wet_or_dry <- apply((wet_or_dry * terr_area) / sum(terr_area) * 100,
                               dim_remain, sum)
  }
  return(result)
}
