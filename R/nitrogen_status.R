#' Status calculation of the nitrogen boundary
#'
#' Planetary Boundary status calculation of the the nitrogen boundary based on
#' a scenario LPJmL run and if `approach == "braun2022_minusref"` a reference
#' LPJmL run.
#'
#' @param files_scenario list with variable names and corresponding file paths
#' (character string) of the scenario LPJmL run. Handled automatically via
#' [`calc_status()`].
#'
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the files_reference LPJmL run. Handled automatically via
#' [`calc_status()`].
#'
#' @param spatial_scale character. Spatial resolution, available options
#'        are `"global"` and `"grid"`
#'
#' @param time_span_scenario time span to use output from the scenario run,
#' e.g. `1982:2011`.
#'
#' @param time_span_reference time span use output from the reference run,
#' e.g. `1901:1930`.
#'
#' @param approach (character string) to be used , currently available
#' approach is `"braun2022"` based on unpublished suggestion by Johanna Braun.
#' Second approach option is `"braun2022_minusref"` to subtract reference run
#' output
#'
#' @param time_series_avg integer. Number of years to be used for the moving
#' average calculation. If `NULL`, all years are averaged for one status
#' calculation, for `1` the whole time span is used to calculate a status time
#' series.
#'
#' @param config_args list of arguments to be passed on from the model
#' configuration.
#'
#' @param thresholds list with highrisk and pb threshold for N concentration
#' (mg N/l) in runoff to surface water
#' Default: highrisk = 5, pb = 2
#' (based on Schulte-Uebbing et al. 2022,
#' https://doi.org/10.1038/s41586-022-05158-2:
#' "we used a threshold for N concentration in run-off to surface water. This
#' threshold was set to 5.0 mgN/l, based on the assumption that on average 50%
#' of N entering surface water is removed through retention and sedimentation"))
#'
#' @param cut_arid double. Exclude boundary calculations below the defined
#' threshold for aridity (annual precipitation / annual potential
#' evapotranspiration); Default: 0.2
#'
#' @param cut_runoff double. Exclude boundary calculations below the defined
#' runoff threshold; Default: 0 mm per year (no treshold)
#'
#' @param with_groundwater_denit logical. Include global assumptions made on
#' groundwater denitrification losses. Defaults to TRUE ( = simulated leaching
#' is multiplied with 0.71 based on simulated denitrification losses in ground
#' water from Bouwman et al 2013)
#'
#'@return Object of class `control_variable` with the boundary status of the
#' nitrogen boundary.
#'
#' @examples
#' \dontrun{
#' boundary_status <- calc_status(
#'   boundary = "nitrogen",
#'   config_scenario = "path/to/config_scenario.json",
#'   config_reference = "path/to/config_reference.json",
#'   spatial_scale = "global",
#'   time_span_scenario = 1901:2019,
#'   time_span_reference = 1901:1930
#' )
#' }
#'
#' @md
#' @importFrom future %<-%
#' @export
nitrogen_status <- function(
  files_scenario,
  files_reference,
  spatial_scale = "grid",
  time_span_scenario = 1982:2011,
  time_span_reference = time_span_scenario,
  approach = "braun2022",
  time_series_avg = NULL,
  config_args = list(),
  thresholds = NULL,
  cut_arid = 0.2,
  cut_runoff = 0,
  with_groundwater_denit = TRUE
) {

  if (spatial_scale == "grid") {
    # verify available methods
    approach <- match.arg(
      approach,
      c("braun2022", "braun2022_minusref")
    )

    # sub function to be used for scenario and reference run
    # (braun2022_minusref)
    calc_nitrogen_leach <- function(
      path_data,
      time_span,
      with_groundwater_denit,
      time_series_avg,
      time_repeat = NULL
    ) {


      # please R CMD check for use of future operator
      runoff <- NULL
      # read runoff ---------------------------------------------------------- #
      runoff %<-% read_io_format(
        file = path_data$runoff,
        timespan = time_span,
        aggregate = list(month = sum),
        spatial_subset = config_args$spatial_subset
      )

      # please R CMD check for use of future operator
      leaching <- NULL
      # read leaching -------------------------------------------------------- #
      leaching %<-% read_io_format(
        file = path_data$leaching,
        timespan = time_span,
        aggregate = list(month = sum),
        spatial_subset = config_args$spatial_subset
      )

      # ---------------------------------------------------------------------- #
      # please R CMD check for use of future operator
      avg_runoff <- avg_leaching <- NULL
      # average runoff
      avg_runoff %<-% aggregate_time(
        x = runoff,
        time_series_avg = time_series_avg,
        time_repeat = time_repeat
      )

      # average leaching
      avg_leaching %<-% aggregate_time(
        x = leaching,
        time_series_avg = time_series_avg,
        time_repeat = time_repeat
      )

      if (with_groundwater_denit) {
        # temporary solution using shares of global losses after
        # Bouwman et al 2013
        # (figure 4) https://doi.org/10.1098/rstb.2013.0112
        # net_flow = net soil N flow + net groundwater N flow (inflow - outflow)
        net_flow <- 65 + 15 - 6
        # gross_flow = gross soil & groundwater N flow + gross surface N flow
        gross_flow <- 93 + 11
        loss_factor <- net_flow / gross_flow
      } else {
        loss_factor <- 1
      }

      n_conc <- ifelse(
        avg_runoff > 0,
        (avg_leaching * 1e3 * loss_factor) /
          (avg_runoff), 0
      )

      return(n_conc)
    }

    # apply defined approach
    switch(approach,
      braun2022 = {
        n_conc <- calc_nitrogen_leach(
          path_data = files_scenario,
          time_span = time_span_scenario,
          with_groundwater_denit = with_groundwater_denit,
          time_series_avg = time_series_avg
        )
      },
      braun2022_minusref = {

        # calculate leaching concentration and loss rate for scenario output
        n_conc_scenario <- calc_nitrogen_leach(
          path_data = files_scenario,
          time_span = time_span_scenario,
          with_groundwater_denit = with_groundwater_denit,
          time_series_avg = time_series_avg
        )

        # calculate leaching concentration and loss rate for reference output
        n_conc_reference <- calc_nitrogen_leach(
          path_data = files_reference,
          time_span = time_span_reference,
          with_groundwater_denit = with_groundwater_denit,
          time_series_avg = NULL,
          time_repeat = length(time_span_scenario)
        )

        # subtract scenario leaching concentration and loss rate from reference
        n_conc <- n_conc_scenario - n_conc_reference
        n_conc[n_conc < 0] <- 0
      }
    )

    # please R CMD check for use of future operator
    pet <- NULL
    # read potential evapotranspiration -------------------------------------- #
    pet %<-% read_io_format(
      file = files_scenario$pet,
      timespan = time_span_scenario,
      aggregate = list(month = sum, band = sum),
      spatial_subset = config_args$spatial_subset
    )

    # please R CMD check for use of future operator
    prec <- NULL
    # read precipitation ----------------------------------------------------- #
    prec %<-% read_io_format(
      file = files_scenario$prec,
      timespan = time_span_scenario,
      aggregate = list(month = sum, band = sum, day = sum),
      spatial_subset = config_args$spatial_subset
    )
    # ------------------------------------------------------------------------ #
    # please R CMD check for use of future operator
    avg_pet <- NULL
    # average pet
    avg_pet <- aggregate_time(
      x = pet,
      time_series_avg = time_series_avg
    )
    # please R CMD check for use of future operator
    avg_prec <- NULL
    # average precipitation
    avg_prec <- aggregate_time(
      x = prec,
      time_series_avg = time_series_avg
    )

    # calculate global aridity index (AI) as an indicator for a level under
    #   which the calculation of leaching just cannot show realistic behavior,
    #   see also on the AI: https://doi.org/10.6084/m9.figshare.7504448.v4%C2%A0
    #   & (first descr.) https://wedocs.unep.org/xmlui/handle/20.500.11822/30300
    #   on nitrogen processes in arid areas:
    #   https://www.jstor.org/stable/45128683
    #     -> indicates boundary to "arid" as thresholds (=< 0.2)
    #   on "arid threshold" (indirectly): https://doi.org/10.1038/ncomms5799
    #     -> threshold for behaviour change in nitrogen cycling (=< 0.32)
    global_aridity_index <- avg_prec / avg_pet + 1e-9

    # ------------------------------------------------------------------------
    # please R CMD check for use of future operator
    runoff <- NULL
    # read runoff
    runoff %<-% read_io_format(
      file = files_scenario$runoff,
      timespan = time_span_scenario,
      aggregate = list(month = sum),
      spatial_subset = config_args$spatial_subset
    )

    # please R CMD check for use of future operator
    runoff_annual <- NULL
    # average runoff
    runoff_annual %<-% aggregate_time(
      x = runoff,
      time_series_avg = time_series_avg
    )

    # to display arid cells (leaching behaviour threshold) in other color (grey)
    cells_arid <- array(FALSE,
                        dim = dim(n_conc),
                        dimnames = dimnames(n_conc))
    cells_arid[which(global_aridity_index <= cut_arid)] <- TRUE

    # to display cells with low runoff in other color (grey):
    cells_low_runoff <- array(FALSE,
                              dim = dim(n_conc),
                              dimnames = dimnames(n_conc))
    cells_low_runoff[which(runoff_annual <= cut_runoff)] <- TRUE

    # add thresholds as attributes
    if ("band" %in% names(dimnames(n_conc))) {
      control_variable <- abind::adrop(
        n_conc,
        drop = which(names(dimnames(n_conc)) == "band")
      )
    } else {
      control_variable <- n_conc
    }

    # non applicable cells
    control_variable[cells_arid] <- NA
    control_variable[cells_low_runoff] <- NA

  } else if (spatial_scale == "global") {
    # verify available methods
    approach <- match.arg(approach, c("schulte_uebbing2022"))
    # thresholds from rockström et al. 2023
    # https://doi.org/10.1038/s41586-023-06083-8

    # please R CMD check for use of future operator
    napplied_mg <- NULL
    # read in fertilizer and manure input on managed land
    napplied_mg %<-% read_io_format(
      file = files_scenario$napplied_mg,
      time_span_scenario,
      aggregate = list(month = sum, band = sum),
      spatial_subset = config_args$spatial_subset
    )

    # please R CMD check for use of future operator
    bnf <- NULL
    # read in biological nitrogen fixation on managed land
    bnf %<-% read_io_format(
      file = files_scenario$bnf_mg,
      time_span_scenario,
      aggregate = list(month = sum, band = sum),
      spatial_subset = config_args$spatial_subset
    )

    # please R CMD check for use of future operator
    dep <- NULL
    # read in nitrogen deposition on managed land
    dep %<-% read_io_format(
      file = files_scenario$ndepo_mg,
      time_span_scenario,
      aggregate = list(month = sum, band = sum),
      spatial_subset = config_args$spatial_subset
    )

    # please R CMD check for use of future operator
    flux_estabn <- NULL
    # read in establishemnt input on managed land
    flux_estabn %<-% read_io_format(
      file = files_scenario$flux_estabn_mg,
      time_span_scenario,
      aggregate = list(month = sum, band = sum),
      spatial_subset = config_args$spatial_subset
    )

    # please R CMD check for use of future operator
    harvest <- NULL
    # read in N removal on managed land
    harvest %<-% read_io_format(
      file = files_scenario$harvestn,
      time_span_scenario,
      aggregate = list(month = sum, band = sum),
      spatial_subset = config_args$spatial_subset
    )

    # calc terrestrial area
    terr_area <- lpjmlkit::read_io(
      files_scenario$terr_area
    ) %>%
      conditional_subset(config_args$spatial_subset) %>%
      lpjmlkit::as_array()

    terr_area <- terr_area[, , 1]
    # calc n surplus
    nsurplus <- (napplied_mg + bnf + dep + flux_estabn - harvest) *
      terr_area * 10^-12

    # N surplus on cropland (n inputs minus n harvest)
    # conversion to TgN/year
    #nsurplus <- (total_fert + total_man + total_dep + total_bnf + total_seed -
    #             harvest) * cellarea * 10^-12

    # average over time
    avg_nsurplus <- aggregate_time(
      x = nsurplus,
      time_series_avg = time_series_avg
    )

    # aggregate to global value
    dim_remain <- names(dim(avg_nsurplus))[names(dim(avg_nsurplus)) != "cell"]
    control_variable <- apply(avg_nsurplus, dim_remain, sum, na.rm = TRUE)

  }

  control_variable <- set_attributes(
    control_variable,
    approach,
    "nitrogen",
    spatial_scale,
    thresholds
  )

  return(control_variable)
}

# read file function --------------------------------------------------------- #
read_io_format <- function(
  file,
  timespan,
  aggregate = list(),
  spatial_subset = NULL
) {
  file <- lpjmlkit::read_io(
    file, subset = list(year = timespan)
  ) %>%
    conditional_subset(spatial_subset) %>%
    lpjmlkit::transform(to = c("year_month_day")) %>%
    lpjmlkit::as_array(aggregate = aggregate) %>%
    suppressWarnings()
  return(file)
}
