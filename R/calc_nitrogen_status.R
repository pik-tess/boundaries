#' Calculate the planetary boundary status for the nitrogen boundary
#'
#' Calculate the PB status for the nitrogen boundary based on a scenario LPJmL
#' run and if `method == "braun2022_minusref"` a reference LPJmL run.
#'
#' @param files_scenario list with variable names and corresponding file paths
#' (character string) of the scenario LPJmL run. All needed files are
#' provided in XXX. E.g.: list(leaching = "/temp/leaching.bin.json")
#'
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the reference LPJmL run. All needed files are
#' provided in XXX. E.g.: list(leaching = "/temp/leaching.bin.json"). If not
#' needed for the applied method, set to NULL.
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#' as character string, e.g. `as.character(1982:2011)` (default)
#'
#' @param time_span_reference time span to be used if scenario run used, defined
#' as character string, e.g. `as.character(1901:1930)`. Can differ in offset and
#' length from `time_span_scenario`! If `NULL` value of `time_span_scenario` is
#' used
#'
#' @param method (character string) to be used , currently available
#' method is `"braun2022"` based on unpublished suggestion by Johanna Braun.
#' Second method option is `"braun2022_minusref"` to subtract reference run
#' output
#'
#' @param spatial_scale character. Spatial resolution, available options
#'        are `"global"` and `"grid"`
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
#' @param thresholds list with highrisk and pb threshold for N concentration
#' (mg N/l) in runoff to surface water
#' Default: highrisk = 5, pb = 2
#' (based on Schulte-Uebbing et al. 2022,
#' https://doi.org/10.1038/s41586-022-05158-2:
#' "we used a threshold for N concentration in run-off to surface water. This
#' threshold was set to 5.0 mgN/l, based on the assumption that on average 50%
#' of N entering surface water is removed through retention and sedimentation"))
#'
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#' time series analysis
#'
#' @examples
#' \dontrun{
#'  calc_nitrogen_status(files_scenario, files_reference)
#' }
#'
#' @md
#' @importFrom future %<-%
#' @export
calc_nitrogen_status <- function(files_scenario,
                                 files_reference,
                                 time_span_scenario = as.character(1982:2011),
                                 time_span_reference = NULL,
                                 method = "braun2022",
                                 thresholds = NULL,
                                 spatial_scale = "grid",
                                 cut_arid = 0.2,
                                 cut_runoff = 0,
                                 with_groundwater_denit = TRUE,
                                 avg_nyear_args = list()
                                 ) {

  if (spatial_scale == "grid") {
    # verify available methods
    method <- match.arg(method, c("braun2022",
                                  "braun2022_minusref"))

    # sub function to be used for scenario and reference run
    # (braun2022_minusref)
    calc_nitrogen_leach <- function(path_data,
                                    time_span,
                                    with_groundwater_denit,
                                    avg_nyear_args
                                    ) {

      # read runoff ---------------------------------------------------------- #
      runoff %<-% read_io_format(
        file = path_data$runoff,
        timespan = time_span,
        aggregate = list(month = sum)
      )

      # read leaching -------------------------------------------------------- #
      leaching %<-% read_io_format(
        file = path_data$leaching,
        timespan = time_span,
        aggregate = list(month = sum)
      )

      # ---------------------------------------------------------------------- #
      # average runoff
      avg_runoff %<-% do.call(average_nyear_window,
                              append(list(x = runoff),
                                     avg_nyear_args))

      # average leaching
      avg_leaching %<-% do.call(average_nyear_window,
                                append(list(x = leaching),
                                       avg_nyear_args))

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

    # apply defined method
    switch(method,
      braun2022 = {
        n_conc <- calc_nitrogen_leach(
          path_data = files_scenario,
          time_span = time_span_scenario,
          with_groundwater_denit = with_groundwater_denit,
          avg_nyear_args = avg_nyear_args)
      },
      braun2022_minusref = {

        # calculate leaching concentration and loss rate for scenario output
        n_conc_scenario <- calc_nitrogen_leach(
          path_data = files_scenario,
          time_span = time_span_scenario,
          with_groundwater_denit = with_groundwater_denit,
          avg_nyear_args = avg_nyear_args)


        if (length(time_span_reference) < length(time_span_scenario)) {
          avg_nyear_args["nyear_reference"] <- length(time_span_scenario)
        }

        # calculate leaching concentration and loss rate for reference output
        n_conc_reference <- calc_nitrogen_leach(
          path_data = files_reference,
          time_span = time_span_reference,
          with_groundwater_denit = with_groundwater_denit,
          avg_nyear_args = avg_nyear_args)

        # subtract scenario leaching concentration and loss rate from reference
        n_conc <- n_conc_scenario - n_conc_reference
        n_conc[n_conc < 0] <- 0
      }
    )

    # read potential evapotranspiration -------------------------------------- #
    pet %<-% read_io_format(
      file = files_scenario$pet,
      timespan = time_span_scenario,
      aggregate = list(month = sum)
    )

    # read precipitation ----------------------------------------------------- #
    prec %<-% read_io_format(
      file = files_scenario$prec,
      timespan = time_span_scenario,
      aggregate = list(month = sum)
    )
    # ------------------------------------------------------------------------ #
    # average pet
    avg_pet <- do.call(average_nyear_window,
                       append(list(x = pet),
                              avg_nyear_args))

    # average precipitation
    avg_prec <- do.call(average_nyear_window,
                        append(list(x = prec),
                               avg_nyear_args))

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
    # read runoff
    runoff %<-% read_io_format(
      file = files_scenario$runoff,
      timespan = time_span_scenario,
      aggregate = list(month = sum)
    )

    # average runoff
    runoff_annual %<-% do.call(
      average_nyear_window,
      append(
        list(x = runoff),
        avg_nyear_args
      )
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
    control_variable <- n_conc
    attr(control_variable, "thresholds") <- thresholds
    # non applicable cells
    control_variable[cells_arid] <- NA
    control_variable[cells_low_runoff] <- NA

  } else if (spatial_scale == "global") {
    # verify available methods
    method <- match.arg(method, c("schulte_uebbing2022"))
    # thresholds from rockström et al. 2023
    # https://doi.org/10.1038/s41586-023-06083-8

    # TODO better use cftoutput?
    cftfrac %<-% read_io_format(
      file = files_scenario$cftfrac,
      time_span_scenario,
      aggregate = list(month = sum)
    )
    cftfrac_combined <- apply(cftfrac, c("cell", "year"), sum)

    # read in fertilizer for agricultural land
    fert_agr %<-% read_io_format(
      file = files_scenario$nfert_agr,
      time_span_scenario,
      aggregate = list(month = sum)
    )

    # read in fertilizer for managed grassland
    fert_mgrass %<-% read_io_format(
      file = files_scenario$nfert_mgrass,
      time_span_scenario,
      aggregate = list(month = sum)
    )

    # read in manure for agricultural land
    manure_agr %<-% read_io_format(
      file = files_scenario$nmanure_agr,
      time_span_scenario,
      aggregate = list(month = sum)
    )
    # read in manure for managed grassland
    manure_mgrass %<-% read_io_format(
      file = files_scenario$nmanure_mgrass,
      time_span_scenario
    )

    # read in biological nitrogen fixation
    bnf %<-% read_io_format(
      file = files_scenario$pft_bnf,
      time_span_scenario,
      aggregate = list(month = sum)
    )

    # read in nitrogen deposition
    dep %<-% read_io_format(
      file = files_scenario$ndepos,
      time_span_scenario,
      aggregate = list(month = sum)
    )

    # find band names which start with rainfed or irrigated to exclude
    # natural PFT bands
    bandnames_bnf <- dimnames(bnf)[["band"]][
      grepl("rainfed|irrigated", dimnames(bnf)[["band"]])
    ]
    bnf <- bnf[, , bandnames_bnf, drop = FALSE]

    harvest %<-% read_io_format(
      file = files_scenario$harvest,
      time_span_scenario,
      aggregate = list(month = sum)
    ) # TODO better with asub!

    seed <- read_io_format(
      file = files_scenario$seed,
      time_span_scenario,
      aggregate = list(month = sum)
    )

    flux_estabn_mgrass <- read_io_format(
      file = files_scenario$flux_estabn_mgrass,
      time_span_scenario,
      aggregate = list(month = sum)
    )

    # calc terrestrial area
    # TODO better asub?
    terr_area <- lpjmlkit::read_io(files_scenario$terr_area) %>% as_array()
    terr_area <- terr_area[, , 1]

    # calc n surplus
    nsurplus <- (apply(bnf * cftfrac * terr_area, c("cell", "year"), sum) +
                    apply((seed + flux_estabn_mgrass + fert_agr + manure_agr +
                     fert_mgrass + manure_mgrass - harvest), c("cell", "year"), sum) * terr_area +
                     apply(dep, c("cell", "year"), sum) * cftfrac_combined * terr_area) * 10^-12

    # N surplus on cropland (n inputs minus n harvest)
    # conversion to TgN/year
    #nsurplus <- (total_fert + total_man + total_dep + total_bnf + total_seed -
    #             harvest) * cellarea * 10^-12

    # average over time
    avg_nsurplus <- do.call(average_nyear_window,
                            append(list(x = nsurplus),
                                   avg_nyear_args))

    # aggregate to global value
    dim_remain <- names(dim(avg_nsurplus))[names(dim(avg_nsurplus)) != "cell"]
    control_variable <- apply(avg_nsurplus, dim_remain, sum, na.rm = TRUE)

    attr(control_variable, "thresholds") <- thresholds
    attr(control_variable, "control variable") <- (
      "nitrogen surplus on agricultural land"
    )
  }
  return(control_variable)
}


# read file function --------------------------------------------------------- #
read_io_format <- function(file, timespan, aggregate=list()) {
  file <- lpjmlkit::read_io(
    file, subset = list(year = timespan)
  ) %>%
    lpjmlkit::transform(to = c("year_month_day")) %>%
    lpjmlkit::as_array(aggregate = aggregate) %>%
    suppressWarnings()
  return(file)
}
