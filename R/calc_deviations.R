
calc_deviations <- function(file_scenario,
                            file_reference,
                            grid_path,
                            time_span_scenario = as.character(1982:2011),
                            time_span_reference = NULL,
                            method = "porkka_2023",
                            avg_nyear_args = list(),
                            spatial_resolution
                                   ) {

  # verify available methods
  method <- match.arg(method, c("wang-erlandsson2022", "porkka_2023"))
   # verify available spatial resolution
  spatial_resolution <- match.arg(spatial_resolution, c("grid", "global"))
  # TODO not yet compatible with avg_nyear_args

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

  # grid for cellarea, only needed for global resolution #TODO
  grid <- lpjmlkit::read_io(
      grid_path,
      silent = TRUE
      )
  cell_area <- lpjmlkit::calc_cellarea(grid)

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

  result <- list()
  if (spatial_resolution == "grid") {
    n_depart <- function(x) {
      length(which(x == 1))
    }
    result$dry <- apply(dry, "cell", n_depart)
    result$wet <- apply(wet, "cell", n_depart)

  } else if (spatial_resolution == "global") {
    dim_remain <- names(dim(dry))[names(dim(dry)) != "cell"]
    result$dry <- mean(
                    apply((dry * cell_area) / sum(cell_area) * 100, dim_remain,
                           sum, na.rm = TRUE)
                )
    result$wet <- mean(
                    apply((wet * cell_area) / sum(cell_area) * 100, dim_remain,
                           sum, na.rm = TRUE)
                )
  }
  return(result)
}
