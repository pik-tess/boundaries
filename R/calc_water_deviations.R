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
#' `porkka2024` based on
#' [Porkka et al. 2023](https://eartharxiv.org/repository/view/3438/)
#' (referring to each month of a year; default)
#'
#' @param thresholds list with thresholds to be used to
#' define the safe, increasing risk and high risk zone,
#' For spatial_scale = "global" and "subglobal", this refers to the quantiles of
#' the global/basin area with deviations in the reference period. The default
#' is: c(holocene = 50, pb = 95, highrisk = NULL).
#' If set to NULL, the  default is taken from metric_files.yml
#' For highrisk, the value is currently hard-coded to 0.5
#' (following Richardson et al. 2023)
#'
#' @param spatial_scale character string indicating spatial scale;
#' "global" or "subglobal" for calculation of the share (%) of total
#' global/basin area with deviations (either one value per year
#' (wang-erlandsson2022) or one value per year and month (porkka2024)); "grid"
#' not yet defined
#'
#' @param time_series_avg integer. Number of years to be used for the moving
#' average calculation. If `NULL`, all years are averaged for calculation,
#' for `1` the whole time span is used to calculate a time
#' series.
#'
#' @param config_args list of arguments to be passed on from the model
#' configuration.
#'
#' @param variable character string with the name of the variable to be used
#' for the calculation of the water deviations. Default is "rootmoist"
#'
#' @md
calc_water_deviations <- function(files_scenario,
                                  files_reference,
                                  spatial_scale = "subglobal",
                                  time_span_scenario = NULL,
                                  time_span_reference,
                                  approach = "porkka2024",
                                  thresholds = NULL,
                                  time_series_avg = NULL,
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

  # get basin information for subglobal resolution
  if (spatial_scale == "subglobal") {
    cellinfo <- indexing_drainage(drainage_file = files_scenario$drainage)
    if (!is.null(config_args$spatial_subset)) {
      cellinfo <- lpjmlkit::asub(cellinfo, cell = config_args$spatial_subset)
    }
    endcell <- lpjmlkit::asub(cellinfo, band = "endcell")
  } else {
    endcell <- NULL
  }

  # calculate ice free area
  icefree_area <- calc_icefree_area(
    files_path = files_scenario,
    time_span = time_span_scenario,
    spatial_subset = config_args$spatial_subset
  )

  # calculate the area with wet/dry departues in the reference and scenario
  # please R CMD check for use of future operator
  ref_depart <- scen_depart <- NULL
  ref_depart %<-% calc_departures(
    data = var_reference,
    icefree_area = icefree_area,
    quants = quants,
    spatial_scale = spatial_scale,
    approach = approach,
    endcell = endcell
  )

  scen_depart %<-% calc_departures(
    data = var_scenario,
    icefree_area = icefree_area,
    quants = quants,
    spatial_scale = spatial_scale,
    approach = approach,
    endcell = endcell
  )

  # -------------------------------------------------------------------------- #
  # calculate holocene, pb and highrisk area thresholds
  # please R CMD check for use of future operator
  area_high_risk <- area_pb <- area_holocene <- NULL

  if (spatial_scale == "global") {

    if (approach == "porkka2024") {
      # from monthly to yearly mean area with deviations
      scen_depart$wet_or_dry <- apply(
        scen_depart$wet_or_dry,
        "year",
        mean,
        na.rm = TRUE
      )
      ref_depart$wet_or_dry <- apply(
        ref_depart$wet_or_dry,
        "year",
        mean,
        na.rm = TRUE
      )
    }

    # calculate areas corresponding to the quantiles defined in thresholds
    if (is.null(thresholds[["highrisk"]])) {
      area_high_risk <- 50 # hard-coded based on Richardson et al. 2023 (SI)
    } else {
      area_high_risk %<-% stats::quantile(
        ref_depart$wet_or_dry,
        probs = thresholds[["highrisk"]] / 100,
        na.rm = TRUE
      )
    }

    area_pb %<-% stats::quantile(
      ref_depart$wet_or_dry,
      probs = thresholds[["pb"]] / 100,
      na.rm = TRUE
    )
    area_holocene %<-% stats::quantile(
      ref_depart$wet_or_dry,
      probs = thresholds[["holocene"]] / 100,
      na.rm = TRUE
    )

    control_variable <- aggregate_time(
      x = scen_depart$wet_or_dry,
      time_series_avg = time_series_avg
    )

    control_variable <- set_attributes(
      control_variable,
      thresholds = list(
        holocene = area_holocene,
        pb = area_pb,
        highrisk = area_high_risk
      )
    )

  } else if (spatial_scale == "subglobal") {

    if (approach == "porkka2024") {
      # calculate mean yearly area with transgression for each basin
      dim_remain <- dim(scen_depart$wet_or_dry)[names(dim(scen_depart$wet_or_dry)) != "month"]
      scen_depart$wet_or_dry <- apply(
        scen_depart$wet_or_dry,
        names(dim_remain),
        mean,
        na.rm = TRUE
      )
      dim_remain <- dim(ref_depart$wet_or_dry)[names(dim(ref_depart$wet_or_dry)) != "month"]
      ref_depart$wet_or_dry <- apply(
        ref_depart$wet_or_dry,
        names(dim_remain),
        mean,
        na.rm = TRUE
      )
    }

    # calculate for each basin: thresholds for translation into pb status
    if (is.null(thresholds[["highrisk"]])) {
      # hard-coded based on Richardson et al. 2023 (SI)
      area_high_risk <- rep(50, length(unique(endcell)))
      names(area_high_risk) <- unique(endcell)
    } else {
      area_high_risk %<-% apply(
        ref_depart$wet_or_dry,
        "basin",
        function(x) {
          y <- stats::quantile(x, probs = thresholds[["highrisk"]] / 100, na.rm = TRUE)
          y
        }
      )
    }

    area_pb <- apply(
      ref_depart$wet_or_dry,
      "basin",
      function(x) {
        y <- stats::quantile(x, probs = thresholds[["pb"]] / 100, na.rm = TRUE)
        y
      }
    )

    area_holocene <- apply(
      ref_depart$wet_or_dry,
      "basin",
      function(x) {
        y <- stats::quantile(x, probs = thresholds[["holocene"]] / 100, na.rm = TRUE)
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

    # fill global array with basin values (i.e. assign the basin value to each
    # cell belonging to the respective basin
    scen_departures <- array(
      NA,
      dim = c(
        year = dim(scen_depart$wet_or_dry)[["year"]],
        cell = length(endcell)
      ),
      dimnames = list(
        year = dimnames(scen_depart$wet_or_dry)[["year"]],
        cell = seq_along(endcell)
      )
    )

    basins <- unique(endcell)
    for (i in seq_along(basins)) {
      y <- which(endcell == basins[i])
      scen_departures[, y] <- scen_depart$wet_or_dry[i, ]
    }

    control_variable <- aggregate_time(
      x = scen_departures,
      time_series_avg = time_series_avg
    )

    # create array with thresholds
    threshold_attr <- array(
      NA,
      dim = c(
        year = dim(control_variable)[["year"]],
        cell = length(endcell),
        thresholds = 3
      ),
      dimnames = list(
        year = dimnames(control_variable)[["year"]],
        cell = seq_along(endcell),
        thresholds = names(thresholds)
      )
    )
    threshold_attr[, , "pb"] <- area_pb
    threshold_attr[, , "highrisk"] <- area_high_risk
    threshold_attr[, , "holocene"] <- area_holocene

    control_variable <- set_attributes(
      control_variable,
      thresholds = threshold_attr
    )

    # set ice areas to NA
    control_variable[, is.na(icefree_area[, 1, 1])] <- NA

  }

  return(control_variable)
}




# calculate dry & wet departures and return mean annual/monthly area of
# departure # (global or basin scale)
calc_departures <- function(
  data,
  icefree_area,
  quants,
  spatial_scale,
  approach,
  endcell
) {

  q5_base <- rep(quants[["q5"]], dim(data)["year"])
  q95_base <- rep(quants[["q95"]], dim(data)["year"])

  if (approach == "wang-erlandsson2022") {
    # driest/ wettest month per gridcell for each year -> ignores which month
    dry <- apply(data, c("cell", "year"), min)
    wet <- apply(data, c("cell", "year"), max)
  } else if (approach == "porkka2024") {
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

  control_variable <- list()
  if (spatial_scale == "subglobal") {
    # calculate for each basin and year/month: area with wet/dry departures
    dim_remain <- dim(dry)[names(dim(dry)) != "cell"]
    dimnames_remain <- dimnames(dry)[names(dim(dry)) != "cell"]

    nbasins <- length(unique(endcell))
    control_variable$wet_or_dry <- array(
      NA,
      dim = c(basin = nbasins, dim_remain),
      dimnames = c(list(basin = unique(endcell)), dimnames_remain)
    )
    if (approach == "porkka2024") {
      # calculate mean per basin
      for (b in unique(endcell)) { # go through all basins
        control_variable$wet_or_dry[which(unique(endcell) == b), , ] <- (
          apply(
            (lpjmlkit::asub(dry_or_wet, cell = which(endcell == b), drop = FALSE) *
               lpjmlkit::asub(icefree_area, cell = which(endcell == b))) /
              sum(lpjmlkit::asub(icefree_area, cell = which(endcell == b)),
                  na.rm = TRUE) * 100,
            names(dim_remain), sum, na.rm = TRUE
          )
        )
      }
    } else if (approach == "wang-erlandsson2022") {
      for (b in unique(endcell)) { # go through all basins
        control_variable$wet_or_dry[which(unique(endcell) == b), ] <-
          apply(lpjmlkit::asub(dry_or_wet, cell = which(endcell == b), drop = FALSE) *
                  lpjmlkit::asub(icefree_area, cell = which(endcell == b)) /
                  sum(lpjmlkit::asub(icefree_area, cell = which(endcell == b)),
                      na.rm = TRUE) * 100, names(dim_remain), sum, na.rm = TRUE)
      }
    }

  } else if (spatial_scale == "global") {
    # calculate for year/month: global area with wet/dry departures

    terr_area <- drop(icefree_area)

    dim_remain <- names(dim(dry))[names(dim(dry)) != "cell"]
    control_variable$dry <- apply(
      (dry * terr_area) / sum(terr_area, na.rm = TRUE) * 100,
      dim_remain,
      sum,
      na.rm = TRUE
    )
    control_variable$wet <- apply(
      (wet * terr_area) / sum(terr_area, na.rm = TRUE) * 100,
      dim_remain,
      sum,
      na.rm = TRUE
    )
    wet_or_dry <- dry + wet
    wet_or_dry[wet_or_dry > 1] <- 1
    control_variable$wet_or_dry <- apply(
      (wet_or_dry * terr_area) / sum(terr_area, na.rm = TRUE) * 100,
      dim_remain,
      sum,
      na.rm = TRUE
    )
  }

  return(control_variable)
}

# calculate the baseline quantiles
calc_baseline <- function(file_reference, approach) {

  if (approach == "wang-erlandsson2022") {
    # vectorized calculation of the 5% and 95% quantiles for each cell
    #    only available for matrices (2 dimensions)
    q5_base <- matrixStats::rowQuantiles(
      apply(file_reference, c("cell", "year"), min),
      probs = c(0.05)
    )
    q95_base <- matrixStats::rowQuantiles(
      apply(file_reference, c("cell", "year"), max),
      probs = c(0.95)
    )
    quants <- list(q5 = q5_base, q95 = q95_base)

  } else if (approach == "porkka2024") {
    quants_array <- apply(
      file_reference,
      c("cell", "month"),
      stats::quantile,
      probs = c(0.05, 0.95)
    )
    quants <- list(
      q5 = abind::asub(quants_array, idx = 1, dim = 1),
      q95 = abind::asub(quants_array, idx = 2, dim = 1)
    )
  }

  return(quants)
}


# calculate ice free area, return array with area per cell, cells with
# ice and rock are set to NA
calc_icefree_area <- function(files_path, time_span, spatial_subset) {
  terr_area <- lpjmlkit::read_io(
    files_path$terr_area,
    silent = TRUE
  ) %>%
    lpjmlkit::transform("year_month_day") %>%
    lpjmlkit::as_array()

  # set ice and rock cells to 0 to only refer to ice-free land surface
  # read in fpc
  fpc <- read_io_format(
    files_path$fpc,
    time_span,
    spatial_subset = spatial_subset
  )

  fpc_total <- apply(
    lpjmlkit::asub(fpc, band = -1, drop = FALSE),
    c("cell"),
    sum,
    na.rm = TRUE
  )

  # Temperature
  temp <- {
    lpjmlkit::read_io(
      files_path$temp,
      subset = list(year = time_span),
      silent = TRUE
    ) %>%
      conditional_subset(spatial_subset) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      lpjmlkit::as_array(aggregate =
                           list(month = mean, day = mean, band = mean,
                                year = mean)) %>%
      suppressWarnings()
  }

  # Rocks and Ice
  is_rocks_and_ice <- {
    fpc_total == 0 &
      temp < 0 # defined in classify_biomes
  }
  terr_area[is_rocks_and_ice, , ] <- NA
  return(terr_area)
}


# show the full downstream drainage route for cell ind until final drainage
#   (ocean or inland sink)
show_route <- function(ind, routing_table) {
  if (routing_table[ind] < 1) {
    # can be 0 or -8 -> endcell or nacell
    return(ind)
  } else {
    return(c(ind,
             show_route(routing_table[ind], routing_table)))
  }
}


# calculate the routing information based on the drainage file
indexing_drainage <- function(drainage_file) {

  drainage <- tryCatch(
    lpjmlkit::read_io(drainage_file)$data %>%
      suppressWarnings() %>%
      drop(),
    error = function(e) {
      lpjmlkit::read_io(
        drainage_file,
        datatype = 2
      )$data %>%
        suppressWarnings() %>%
        drop()
    }
  )
  # -------------------------------------------------------------------------- #

  # add 1 since in C indexing starts at 0 but in R at 1
  routing <- drainage[, 1] + 1
  ncell <- length(routing)
  cellinfo <- array(0, dim = c(ncell, 3))
  cellinfo[, 3] <- routing
  for (cell in 1:ncell) {
    route <- show_route(cell, routing)
    # (over)writing all cells from route
    cellinfo[route, 1] <- seq(length(route), 1, -1) # rank (former cellindex)
    cellinfo[cell, 2] <- route[length(route)] # endcell
  }
  dimnames(cellinfo) <- list(cell = 0:(ncell - 1),
                             band = c("rank", "endcell", "routing"))
  return(cellinfo)
}
