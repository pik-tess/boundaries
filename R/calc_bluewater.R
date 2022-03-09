#' Calculate bluewater planetary boundary status
#'
#' Calculate ...
#'
#' @param path_reference
#'
#' @param method
#'
#' @param time_span_scenario
#'
#' @param bin_size
#'
#' @examples
#' \dontrun{
#'  calc_bluewater(path_reference)
#' }
#'
#' @md
#' @export
calc_bluewater <- function(path_scenario,
                           path_reference,
                           time_span_scenario = c(1982, 2011),
                           time_span_reference = NULL,
                           method = "gerten2020",
                           temporal_resolution = "annual",
                           spatial_format = "grid",
                           cut_min = 0.0864, # Q < smaller than 1m³/s
                           avg_nyear_args = list(),
                           # to be replaced internally by lpjmlKit::read_output
                           start_year = 1901,
                           end_year = 2011) {
  # verify available methods
  method <- match.arg(method, c("gerten2020",
                                "steffen2015"))
  # verify available temporal resolution
  temporal_resolution <- match.arg(temporal_resolution, c("annual",
                                                          "monthly"))

  # verify available temporal resolution
  spatial_format <- match.arg(spatial_format, c("grid",
                                                "continent",
                                                "country",
                                                "global"))

  if (is.null(time_span_reference)) {
    time_span_reference <- time_span_scenario
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
  s_path <- file(paste(path_scenario, "discharge.bin", sep = ""), "rb")
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
      efr_uncertain <- calc_efrs(discharge_reference, "vmf_min")
      efr_safe <- calc_efrs(discharge_reference, "vmf_max")

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
      status_frac <- ifelse(uncertainty_zone > 0,
                            safe_space / uncertainty_zone,
                            0)

      if (temporal_resolution == "annual") {
        # to average the ratio only over months which are not "safe"
        status_frac[status_frac <= 0.05] <- NA
        pbw_frac_mean <- apply(
          status_frac,
          names(dim(status_frac))[
            names(dim(status_frac)) %in% c("cells", "years")
          ],
          mean,
          na.rm = TRUE)

      }

    },
    steffen2015 = {
      stop("This method is currently not defined.")
    }
  )
}