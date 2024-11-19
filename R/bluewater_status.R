#' Status calculation of the bluewater boundary
#'
#' Planetary Boundary status calculation of the bluewater boundary (as part of
#' the freshwater boundary) based on a scenario LPJmL run and a reference LPJmL
#' run.
#'
#' @param files_scenario list with variable names and corresponding file paths
#' (character string) of the scenario LPJmL run. Handled automatically via
#' [`calc_status()`].
#'
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the files_reference LPJmL run. Handled automatically
#' via [`calc_status()`].
#'
#' @param spatial_scale character string indicating spatial resolution
#' options: "global", "regional", "grid";
#' for "grid" the approach "gerten2020" is applicable based on EFR calculations;
#' for "global"/"regional" the share (%) of total global/basin area with
#' deviations is calculated
#'
#' @param time_span_scenario time span to use output from the scenario run,
#' e.g. `1982:2011`.
#'
#' @param time_span_reference time span use output from the reference run,
#' e.g. `1901:1930`.
#'
#' @param approach approach (character string) to be used , currently available
#' approach is `"gerten2020"` based on
#' [Gerten et al. 2020](https://doi.org/10.1038/s41893-019-0465-1)
#' for spatial_scale = "grid" and
#' "wang_erlandsson2022" as well as "porkka2024" for
#' spatial_scale = "global" or "regional"
#'
#' @param time_series_avg integer. Number of years to be used for the moving
#' average calculation. If `NULL`, all years are averaged for one status
#' calculation, for `1` the whole time span is used to calculate a status time
#' series.
#'
#' @param config_args list of arguments to be passed on from the model
#' configuration.
#'
#' @param thresholds named character string with thresholds to be used to
#' define the safe, increasing risk and high risk zone, the approach and scale
#' specific default thresholds are defined in metric_files.yml are are applied
#' if thresholds are set to NULL.
#'
#' @param cut_min double. Exclude boundary calculations for discharge < cut_min
#' and dismiss EFR transgresssions if < cut_min for "gerten2020" approach,
#' Default: 0.0864 hm3/day (=1 m3/s)
#'
#'@return Object of class `control_variable` with the boundary status of the
#' bluewater boundary.
#'
#' @examples
#' \dontrun{
#' boundary_status <- calc_status(
#'   boundary = "bluewater",
#'   config_scenario = "path/to/config_scenario.json",
#'   config_reference = "path/to/config_reference.json",
#'   spatial_scale = "global",
#'   time_span_scenario = 1901:2019,
#'   time_span_reference = 1901:1930,
#'   approach = "porkka2024"
#' )
#' }
#'
#' @md
#' @export
bluewater_status <- function(
  files_scenario,
  files_reference,
  spatial_scale,
  time_span_scenario = as.character(1982:2011),
  time_span_reference = time_span_scenario,
  approach = "gerten2020",
  time_series_avg = NULL,
  config_args = list(),
  thresholds = NULL,
  cut_min = 0.0864
) {

  # verify available methods and resolution
  approach <- match.arg(
    approach,
    c("gerten2020", "wang-erlandsson2022", "porkka2024")
  )
  spatial_scale <- match.arg(spatial_scale, c("global", "regional", "grid"))

  # apply defined approach
  if (spatial_scale == "grid") {
    if (approach != "gerten2020") {
      stop(
        "Approach \"gerten2020\" is the only available approach for ",
        "spatial_scale = \"grid\"."
      )
    }
    control_variable <- calc_bluewater_efrs(
      files_scenario = files_scenario,
      files_reference = files_reference,
      spatial_scale = spatial_scale,
      time_span_scenario = time_span_scenario,
      time_span_reference =  time_span_reference,
      approach = approach,
      time_series_avg = time_series_avg,
      config_args = config_args,
      thresholds = thresholds,
      cut_min = cut_min
    )

  } else if (spatial_scale %in% c("regional", "global")) {
    if (!approach %in% c("wang-erlandsson2022", "porkka2024")) {
      stop(
        "Approach \"",
        approach,
        "\" is not available for spatial_scale = ",
        spatial_scale,
        "."
      )
    }

    control_variable <- calc_water_deviations(
      files_scenario = files_scenario,
      files_reference = files_reference,
      spatial_scale = spatial_scale,
      time_span_scenario = time_span_scenario,
      time_span_reference =  time_span_reference,
      approach = approach,
      time_series_avg = time_series_avg,
      config_args = config_args,
      thresholds = thresholds,
      variable = "discharge"
    )
  }
  attr(control_variable, "long_name") <- list_long_name("bluewater")

  return(control_variable)
}


calc_bluewater_efrs <- function(
  files_scenario,
  files_reference,
  spatial_scale,
  time_span_scenario = time_span_scenario,
  time_span_reference = NULL,
  approach = "gerten2020",
  time_series_avg = NULL,
  config_args = list(),
  thresholds = NULL,
  cut_min = 0.0864
) {

  if (spatial_scale != "grid") {
    stop("Approach \"gerten2020\" is only applicable for spatial_scale = \"grid\"") #nolint
  }

  # please R CMD check for use of future operator
  discharge_reference <- NULL
  # reference discharge ---------------------------------------------------- #
  discharge_reference %<-% read_io_format(
    files_reference$discharge,
    time_span_reference,
    aggregate = list(band = sum),
    spatial_subset = config_args$spatial_subset
  )

  # please R CMD check for use of future operator
  discharge_scenario <- NULL
  # scenario discharge ----------------------------------------------------- #
  discharge_scenario %<-% read_io_format(
    files_scenario$discharge,
    time_span_scenario,
    aggregate = list(band = sum),
    spatial_subset = config_args$spatial_subset
  )

  # ------------------------------------------------------------------------ #
  # please R CMD check for use of future operator
  avg_discharge_scenario <- NULL
  # average discharge scenario
  avg_discharge_scenario %<-% aggregate_time(
    discharge_scenario,
    time_series_avg = time_series_avg
  )

  # please R CMD check for use of future operator
  avg_discharge_reference <- NULL
  # average discharge reference
  avg_discharge_reference %<-% aggregate_time(
    discharge_reference,
    time_repeat = if (is.null(time_series_avg)) NULL else dim(avg_discharge_scenario)["year"]
  )

  # please R CMD check for use of future operator
  efr_uncertain <- efr_safe <- NULL
  # calc efrs for vmf_min and vmf_max
  efr_uncertain %<-% calc_efrs(
    avg_discharge_reference,
    "vmf_min"
  )

  efr_safe %<-% calc_efrs(
    avg_discharge_reference,
    "vmf_max"
  )
  # calculation of EFR transgressions = EFR deficits in LU run
  efr_deficit <- efr_safe - avg_discharge_scenario
  # dismiss small EFR deficits
  efr_deficit[efr_deficit < cut_min] <- 0
  # calculation of uncertainty zone
  uncertainty_zone <- efr_safe - efr_uncertain

  # calculate boundary status based on transgression to uncertainty ratio
  # (as in Steffen 2015; degree to which EFRs are undermined)

  status_frac_monthly <- ifelse(uncertainty_zone > 0,
                                efr_deficit / uncertainty_zone * 100,
                                0)

  third_dim <- names(dim(status_frac_monthly))[
    !names(dim(status_frac_monthly)) %in% c("cell", "month")
  ] %>% {
    if (rlang::is_empty(.)) NULL else .
  }

  # to average the ratio only over months which are not "safe"
  status_frac_monthly[status_frac_monthly <= thresholds[["pb"]]] <- NA

  control_variable <- apply(status_frac_monthly,
                            names(dim(status_frac_monthly))[
                              names(dim(status_frac_monthly)) %in%
                                c("cell", third_dim)
                            ],
                            mean,
                            na.rm = TRUE)
  # set cells with NA (all months safe) to 0
  control_variable[is.na(control_variable)] <- 0

  # omit boundary status calculation if PNV discharge is < cut_min
  # (marginal discharge)
  cells_marginal_discharge <- array(FALSE,
                                    dim = dim(control_variable),
                                    dimnames = dimnames(control_variable))
  cells_marginal_discharge[
    which(
      apply(
        avg_discharge_reference, c("cell", third_dim), mean
      ) < cut_min
    )
  ] <- TRUE
  control_variable[cells_marginal_discharge] <- NA

  # if ratio is above >5%: within uncertainty range (yellow)
  # if ratio is above >75% transgression (red)
  # define PB thresholds as attributes

  attr(control_variable, "control_variable") <-
    "EFR transgression to uncertainty ratio"
  attr(control_variable, "thresholds") <- thresholds
  attr(control_variable, "spatial_scale") <- spatial_scale
  attr(control_variable, "unit") <- list_unit("bluewater", approach,
                                              spatial_scale)
  attr(control_variable, "long_name") <- list_long_name("bluewater")

  class(control_variable) <- c("control_variable")
  return(control_variable)
}
