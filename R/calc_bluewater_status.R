#' Calculate the planetary boundary status for the bluewater boundary
#'
#' Calculate the PB status for the bluewater (former freshwater) boundary based
#' on a scenario LPJmL run and a reference LPJmL run.
#'
#' @param files_scenario list with variable names and corresponding file paths
#' (character string) of the scenario LPJmL run. All needed files are
#' provided in metric_files.yaml.
#' E.g.: list(leaching = "/temp/leaching.bin.json")
#'
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the reference LPJmL run. All needed files are
#' provided in metric_files.yml.
#' E.g.: list(leaching = "/temp/leaching.bin.json").
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#' as character string
#'
#' @param time_span_reference time span to be used for the scenario run, defined
#' as character string, e.g. `as.character(1901:1930)`. Can differ in offset and
#' length from `time_span_scenario`! If `NULL` value of `time_span_scenario` is
#' used
#'
#' @param method method (character string) to be used , currently available
#' method is `"gerten2020"` based on
#' [Gerten et al. 2020](https://doi.org/10.1038/s41893-019-0465-1)
#' for spatial_scale = "grid" and
#' "wang_erlandsson2022" as well as "porkka2023" for
#' spatial_scale = "global" or "subglobal"
#'
#' @param cut_min double. Exclude boundary calculations for discharge < cut_min
#' and dismiss EFR transgresssions if < cut_min for "gerten2020" method,
#' Default: 0.0864 hm3/day (=1 m3/s)
#'
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#' time series analysis
#'
#' @param spatial_scale character string indicating spatial resolution
#' options: "global", "subglobal", "grid";
#' for "grid" the method "gerten2020" is applicable based on EFR calculations;
#' for "global"/"subglobal" the share (%) of total global/basin area with
#' deviations is calculated
#'
#' @param thresholds named character string with thresholds to be used to
#' define the safe, increasing risk and high risk zone, the method and scale
#' specific default thresholds are defined in metric_files.yml are are applied
#' if thresholds are set to NULL. 
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
                                  time_span_scenario = NULL,
                                  time_span_reference = time_span_scenario,
                                  method = "gerten2020",
                                  cut_min = 0.0864,
                                  avg_nyear_args = list(),
                                  spatial_scale,
                                  thresholds = NULL) {

  # verify available methods and resolution
  method <- match.arg(method, c("gerten2020",
                                "wang-erlandsson2022",
                                "porkka2023"))
  spatial_scale <- match.arg(spatial_scale, c("global", "subglobal", "grid"))

  # apply defined method
  if (method == "gerten2020") {
    if (spatial_scale != "grid") {
      stop("Method \"gerten2020\" is only applicable for spatial_scale = \"grid\"") #nolint
    }
    # reference discharge ---------------------------------------------------- #
    discharge_reference <- lpjmlkit::read_io(
      files_reference$discharge, subset = list(year = time_span_reference),
      silent = TRUE
    ) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      lpjmlkit::as_array(aggregate = list(band = sum))

    # scenario discharge ----------------------------------------------------- #
    discharge_scenario <- lpjmlkit::read_io(
      files_scenario$discharge, subset = list(year = time_span_scenario),
      silent = TRUE
    ) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      lpjmlkit::as_array(aggregate = list(band = sum))

    # ------------------------------------------------------------------------ #

    # TODO understand what this does
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

    pb_status <- apply(status_frac_monthly,
                       names(dim(status_frac_monthly))[
                         names(dim(status_frac_monthly)) %in%
                           c("cell", third_dim)
                       ],
                       mean,
                       na.rm = TRUE)
    # set cells with NA (all months safe) to 0
    pb_status[is.na(pb_status)] <- 0

    # check if vector was returned (loss of dimnames) -> reconvert to array
    #TODO check if necessary
    if (is.null(dim(pb_status))) {
      pb_status <- array(
        pb_status,
        dim = c(cell = dim(status_frac_monthly)[["cell"]], 1),
        dimnames = list(cell = dimnames(status_frac_monthly)[["cell"]], 1)
      )
    }
    # omit boundary status calculation if PNV discharge is < cut_min
    # (marginal discharge)
    cells_marginal_discharge <- array(FALSE,
                                      dim = dim(pb_status),
                                      dimnames = dimnames(pb_status))
    cells_marginal_discharge[
      which(
        apply(
          avg_discharge_reference, c("cell", third_dim), mean
        ) < cut_min
      )
    ] <- TRUE
    pb_status[cells_marginal_discharge] <- NA

    # if ratio is above >5%: within uncertainty range (yellow)
    # if ratio is above >75% transgression (red)
    # define PB thresholds as attributes

    attr(pb_status, "thresholds") <- thresholds

  } else if (method %in% c("wang-erlandsson2022", "porkka2023")) {
    if (spatial_scale == "grid") {
      stop("Method is only applicable for spatial_scale = \"global\" or \"subglobal\"") #nolint
    }
    #TODO also account for cut_min?
    pb_status <- calc_water_status(
      file_scenario = files_scenario$discharge,
      file_reference = files_reference$discharge,
      terr_area_path = files_reference$terr_area,
      drainage_path = files_reference$drainage,
      time_span_scenario = time_span_scenario,
      time_span_reference =  time_span_reference,
      method = method,
      avg_nyear_args = avg_nyear_args,
      spatial_scale = spatial_scale,
      thresholds = thresholds
    )
  }

  return(pb_status)
}
