#' Calculate the planetary boundary status for the bluewater boundary
#'
#' Calculate the PB status for the bluewater (former freshwater) boundary based
#' on a scenario LPJmL run and a reference LPJmL run.
#'
#' @param path_scenario output directory (character string) of the scenario
#' LPJmL run where binary files (soon with metafiles) are written
#'
#' @param path_reference output directory (character string) of the reference
#' LPJmL run where binary files (soon with metafiles) are written
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#' as an integer vector, e.g. `1982:2011` (default)
#'
#' @param time_span_reference time span to be used for the scenario run, defined
#' as an integer vector, e.g. `1901:1930`. Can differ in offset and length from
#' `time_span_scenario`! If `NULL` value of `time_span_scenario` is used
#'
#' @param method method (character string) to be used , currently available
#' method is `c("gerten2020")` based on
#' [Gerten et al. 2020](https://doi.org/10.1038/s41893-019-0465-1).
#'
#' @param temporal_resolution character. Temporal resolution, available options
#' are `"annual"` (default) and `"monthly"`.
#'
#' @param cut_min double. Exclude boundary calculations for Q < cut_min
#'
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#' time series analysis.
#'
#' @examples
#' \dontrun{
#'  calc_bluewater(path_scenario, path_reference)
#' }
#'
#' @md
#' @export
calc_bluewater_status <- function(path_scenario,
                                  path_reference,
                                  time_span_scenario = c(1982, 2011),
                                  time_span_reference = NULL,
                                  method = "gerten2020",
                                  temporal_resolution = "annual",
                                  # Q < smaller than 1m³/s
                                  cut_min = 0.0864,
                                  avg_nyear_args = list(),
                                  # to be replaced by lpjmlKit::read_output
                                  start_year = 1901) {
  # verify available methods
  method <- match.arg(method, c("gerten2020",
                                "steffen2015"))
  # verify available temporal resolution
  temporal_resolution <- match.arg(temporal_resolution, c("annual",
                                                          "monthly"))

  # check time_spans of scenario and reference runs
  if (is.null(time_span_reference)) {
    time_span_reference <- time_span_scenario
    nyear_ref <- NULL
  } else {
    if (diff(time_span_reference) > diff(time_span_scenario)) {
      stop(paste0("time_span_reference is longer than time_span_scenario.",
                  "Define a time_span_reference that is shorter than",
                  "time_span_scenario"))
    } else if (diff(time_span_reference) < diff(time_span_scenario)) {
      nyear_ref <- length(time_span_scenario[1]:time_span_scenario[2])
    } else {
      nyear_ref <- NULL
    }
  }

  # reference discharge
  # TO BE REPLACED BY lpjmlKit::read_output ---------------------------------- #
  #   hardcoded values to be internally replaced
  nstep <- 12
  nbands <- 1
  ncell <- 67420
  size <- 4
  bl_file <- file(paste(path_reference, "discharge.bin", sep = "/"), "rb")
  seek(bl_file,
       where = (time_span_reference[1] - start_year) *
               nstep * nbands * ncell * size,
       origin = "start")
  discharge_reference <- readBin(bl_file,
                                double(),
                                n = (ncell * nstep * nbands *
                                     (time_span_reference[2] -
                                      time_span_reference[1] + 1)),
                                size = size)
  close(bl_file)
  dim(discharge_reference) <- c(cells = ncell,
                                months = nstep,
                                years = (time_span_reference[2] -
                                         time_span_reference[1] + 1))

  dimnames(discharge_reference) <- list(cells = seq_len(ncell),
                                        months = seq_len(nstep),
                                        years = seq(time_span_scenario[1],
                                                    time_span_scenario[2]))
  # -------------------------------------------------------------------------- #

  # average discharge reference
  avg_discharge_reference <- do.call(average_nyear_window,
                                     append(list(x = discharge_reference,
                                                 nyear_reference = nyear_ref),
                                            avg_nyear_args))

  # scenario discharge
  # TO BE REPLACED BY lpjmlKit::read_output ---------------------------------- #
  #   hardcoded values to be internally replaced
  s_path <- file(paste(path_scenario, "discharge.bin", sep = "/"), "rb")
  seek(s_path,
       where = (time_span_scenario[1] - start_year) *
               nstep * nbands * ncell * size,
       origin = "start")
  discharge_scenario <- readBin(s_path,
                               double(),
                                n = (ncell * nstep * nbands *
                                     (time_span_scenario[2] -
                                      time_span_scenario[1] + 1)),
                               size = size)
  close(s_path)
  dim(discharge_scenario) <- c(cells = ncell,
                               months = nstep,
                               years = (time_span_scenario[2] -
                                        time_span_scenario[1] + 1))
  dimnames(discharge_scenario) <- list(cells = seq_len(ncell),
                                       months = seq_len(nstep),
                                       years = seq(time_span_scenario[1],
                                                   time_span_scenario[2]))
  # -------------------------------------------------------------------------- #

  # average discharge reference
  avg_discharge_scenario <- do.call(average_nyear_window,
                                    append(list(x = discharge_scenario),
                                           avg_nyear_args))

  # apply defined method
  switch(method,
    # "gerten2020" - Gerten et al. 2020
    gerten2020 = {
      # calc efrs for vmf_min and vmf_max
      efr_uncertain <- calc_efrs(discharge_reference,
                                 "vmf_min",
                                 avg_nyear_args)
      efr_safe <- calc_efrs(discharge_reference,
                            "vmf_max",
                            avg_nyear_args)
      # calculation of EFR transgressions = EFR deficits in LU run
      safe_space <- efr_safe - avg_discharge_scenario

      # calculation of uncertainty zone
      uncertainty_zone <- efr_safe - efr_uncertain

      # dismiss boundary status calculation if PNV discharge is < 0.1 m^3/s
      uncertainty_zone[avg_discharge_reference < cut_min] <- 0
      safe_space[safe_space < cut_min] <- 0
      safe_space[avg_discharge_reference < cut_min] <- 0

      # calculate boundary status based on transgression to uncertainty ratio
      #   as in Steffen 2015: degree to which EFRs are undermined: expressed as
      #   the transgression-to uncertainty ratio
      #   if ratio is above >5%: within uncertainty range (yellow)
      #   if ratio is above >75% transgression (red)
      status_frac_monthly <- ifelse(uncertainty_zone > 0,
                                    safe_space / uncertainty_zone,
                                    0)

      third_dim <- names(dim(status_frac_monthly))[
        !names(dim(status_frac_monthly)) %in% c("cells", "months")
      ] %>%
      ifelse(length(.) == 0, NA, .)

      if (temporal_resolution == "annual") {
        # to average the ratio only over months which are not "safe"
        status_frac_monthly[status_frac_monthly <= 0.05] <- NA
        status_frac <- apply(
          status_frac_monthly,
          names(dim(status_frac_monthly))[
            names(dim(status_frac_monthly)) %in% na.omit(c("cells", third_dim))
          ],
          mean,
          na.rm = TRUE)
      } else {
        status_frac <- status_frac_monthly
      }

      # to display cells with marginal discharge in other color (grey):
      cells_maginal_discharge <- array(0, ncell)
      cells_maginal_discharge[
        which(apply(avg_discharge_scenario, 1, mean) < cut_min)] <- (-1)

      # calc irrigation mask to exclude non irrigated basins
      irrmask_basin <- calc_irrigation_mask(path_scenario,
                                            time_span = time_span_scenario,
                                            avg_nyear_args = avg_nyear_args,
                                            start_year = start_year)
      # init pb_status based on status_frac
      pb_status <- status_frac
      # high risk
      pb_status[status_frac >= 0.75] <- 3
      # increasing risk
      pb_status[status_frac < 0.75 & status_frac >= 0.05] <- 2
      # safe zone
      pb_status[status_frac < 0.05] <- 1
      pb_status[irrmask_basin == 0] <- 1
      pb_status[is.na(status_frac)] <- 1
      # non applicable cells
      pb_status[cells_maginal_discharge < 0] <- 0
    },
    steffen2015 = {
      stop("This method is currently not defined.")
    }
  )
  return(pb_status)
}