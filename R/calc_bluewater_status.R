#' Calculate the planetary boundary status for the bluewater boundary
#'
#' Calculate the PB status for the bluewater (former freshwater) boundary based
#' on a scenario LPJmL run and a reference LPJmL run.
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
#' @param time_span_reference time span to be used for the scenario run, defined
#' as character string, e.g. `as.character(1901:1930)`. Can differ in offset and
#' length from `time_span_scenario`! If `NULL` value of `time_span_scenario` is
#' used
#'
#' @param method method (character string) to be used , currently available
#' method is `c("gerten2020")` based on
#' [Gerten et al. 2020](https://doi.org/10.1038/s41893-019-0465-1).
#'
#' @param cut_min double. Exclude boundary calculations for Q < cut_min and
#' dismiss EFR transgresssions if < cut_min for "gerten2020" method,
#' Default: 0.0864 hm3/day (=1 m³/s)
#'
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#' time series analysis
#'
#' @param irrmask_basin logical, if true: all cells in river basins without
#' irrigation will be masked (= no boundary transgression)
#'
#' @param spatial_resolution character string indicating spatial resolution
#'        either "grid" for calculation of number of years with transgression
#'        (for wang-erlandsson2022: dim(ncell, nyears);
#'         for porkka_2023: dim(ncell, nyears, months)) or
#'        "global" for calculation of the share (%) of total global area with
#'        deviations (either one value per year (wang-erlandsson2022) or one
#'        value per year and month (porkka_2023)) - note: not applied for
#'        the method "gerten2020" (only at the grid cell level)
#'
#'@return todo: describe returned object
#'
#' @examples
#' \dontrun{
#'  calc_bluewater_status(files_scenario, files_reference)
#' }
#'
#' @md
#' @export
calc_bluewater_status <- function(files_scenario,
                                  files_reference,
                                  time_span_scenario = as.character(1982:2011),
                                  time_span_reference = time_span_scenario,
                                  method = "gerten2020",
                                  cut_min = 0.0864,
                                  avg_nyear_args = list(),
                                  irrmask_basin = FALSE,
                                  spatial_resolution) {
  # verify available methods
  method <- match.arg(method, c("gerten2020",
                                "wang-erlandsson2022",
                                "porkka_2023"))

  # apply defined method
  if (method == "gerten2020") {
    if (spatial_resolution == "global") {
      stop(paste0("Global resolution not yet defined for method ", method, ". ",
                  "For aggregation to global level use wang-erlandsson2020 or ",
                  "porkka_2023"))
    }
    # reference discharge ---------------------------------------------------- #
    discharge_reference <- lpjmlkit::read_io(
      files_reference$discharge, subset = list(year = time_span_reference)
      ) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      lpjmlkit::as_array(aggregate = list(band = sum)) %>%
      suppressWarnings()
    #TODO add warning, if discharge output is not monthly but yearly?
    #TODO not yet working for header files (month dimension is calles band)
    # scenario discharge ----------------------------------------------------- #

    discharge_scenario <- lpjmlkit::read_io(
      files_scenario$discharge, subset = list(year = time_span_scenario)
      ) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      lpjmlkit::as_array(aggregate = list(band = sum)) %>%
      suppressWarnings()
    # ------------------------------------------------------------------------ #

    if (length(time_span_reference) < length(time_span_scenario)) {
      nyear_ref <- length(time_span_scenario)
    } else {
      nyear_ref <- NULL
    }

    # average discharge reference
    avg_discharge_reference <- do.call(average_nyear_window,
                                       append(list(x = discharge_reference,
                                                   nyear_reference = nyear_ref),
                                              avg_nyear_args))

    # average discharge scenario
    avg_discharge_scenario <- do.call(average_nyear_window,
                                      append(list(x = discharge_scenario),
                                             avg_nyear_args))

    # calc efrs for vmf_min and vmf_max
    efr_uncertain <- calc_efrs(discharge_reference,
                                   "vmf_min",
                                   avg_nyear_args)
    efr_safe <- calc_efrs(discharge_reference,
                              "vmf_max",
                              avg_nyear_args)
    # calculation of EFR transgressions = EFR deficits in LU run
    efr_deficit <- efr_safe - avg_discharge_scenario
    # dismiss small EFR deficits #TODO check relevance
    efr_deficit[efr_deficit < cut_min] <- 0

    # calculation of uncertainty zone
    uncertainty_zone <- efr_safe - efr_uncertain

    # calculate boundary status based on transgression to uncertainty ratio
    # (as in Steffen 2015; degree to which EFRs are undermined)

    status_frac_monthly <- ifelse(uncertainty_zone > 0,
                                    efr_deficit / uncertainty_zone,
                                    0)

    third_dim <- names(dim(status_frac_monthly))[
      !names(dim(status_frac_monthly)) %in% c("cell", "month")
    ] %>% {
      if (rlang::is_empty(.)) NULL else .
    }

    # to average the ratio only over months which are not "safe"
    status_frac_monthly[status_frac_monthly <= 0.05] <- NA
    status_frac <- apply(
      status_frac_monthly,
      names(dim(status_frac_monthly))[
        names(dim(status_frac_monthly)) %in% c("cell", third_dim)
      ],
      mean,
      na.rm = TRUE)
    # set cells with NA (all months safe) to 0
    status_frac[is.na(status_frac)] <- 0

    # check if vector was returned (loss if dimnames) -> reconvert to array
    if (is.null(dim(status_frac))) {
      status_frac <- array(
        status_frac,
        dim = c(cell = dim(status_frac_monthly)[["cell"]], 1),
        dimnames = list(cell = dimnames(status_frac_monthly)[["cell"]], 1)
      )
    }

    # ommit boundary status calculation if PNV discharge is < cut_min
    # (marginal discharge)
    cells_marginal_discharge <- array(FALSE,
                                     dim = dim(status_frac),
                                     dimnames = dimnames(status_frac))
    cells_marginal_discharge[
      which(
        apply(
          avg_discharge_scenario, c("cell", third_dim), mean
        ) < cut_min
      )
    ] <- TRUE
    status_frac[cells_marginal_discharge] <- NA

    # ommit boundary status calculation in basins without irrigation?
    if (irrmask_basin) {
      # calc irrigation mask to exclude non irrigated basins
      irrmask_basin <- calc_irrigation_mask(files_scenario,
                                            time_span = time_span_scenario,
                                            avg_nyear_args = avg_nyear_args)
      status_frac[irrmask_basin == 0] <- NA
    }

    #   if ratio is above >5%: within uncertainty range (yellow)
    #   if ratio is above >75% transgression (red)
    # TODO define in parameter?
    # define PB thresholds as attributes
    attr(status_frac, "thresholds") <- c(holocene = 0, pb = 0.05,
                                         highrisk = 0.75)

    pb_status <- status_frac

  } else if (method %in% c("wang-erlandsson2022", "porkka_2023")) {
    #TODO also account for cut_min?
    pb_status <- calc_water_status(
     file_scenario = files_scenario$discharge,
     file_reference = files_reference$discharge,
     grid_path = files_reference$grid,
     time_span_scenario = time_span_scenario,
     time_span_reference =  time_span_reference,
     method = method,
     avg_nyear_args = avg_nyear_args,
     spatial_resolution = spatial_resolution
   )
  }

  return(pb_status)
}
