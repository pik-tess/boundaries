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
#' @param time_aggregation_args list of arguments to be passed to
#' \link[boundaries]{aggregate_time} (see for more info). To be used for
#' time series analysis
#'
#' @param config_args list of arguments to be passed on from the model
#' configuration.
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
                                  spatial_scale,
                                  time_span_scenario = time_span_scenario,
                                  time_span_reference = NULL,
                                  method = "gerten2020",
                                  time_aggregation_args = list(),
                                  config_args = list(),
                                  thresholds = NULL,
                                  cut_min = 0.0864) {

  # verify available methods and resolution
  method <- match.arg(method, c("gerten2020",
                                "wang-erlandsson2022",
                                "porkka2023"))
  spatial_scale <- match.arg(spatial_scale, c("global", "subglobal", "grid"))

  # apply defined method
  if (spatial_scale == "grid") {
    if (method != "gerten2020") {
      stop(
        "Method \"gerten2020\" is the only available method for ",
        "spatial_scale = \"grid\"."
      )
    }
    control_variable <- calc_bluewater_efrs(
      files_scenario = files_scenario,
      files_reference = files_reference,
      spatial_scale = spatial_scale,
      time_span_scenario = time_span_scenario,
      time_span_reference =  time_span_reference,
      method = method,
      time_aggregation_args = time_aggregation_args,
      config_args = config_args,
      thresholds = thresholds
    )

  } else if (spatial_scale %in% c("subglobal", "global")) {
    if (!method %in% c("wang-erlandsson2022", "porkka2023")) {
      stop(
        "Method \"",
        method,
        "\" is not available for spatial_scale = ",
        spatial_scale,
        "."
      )
    }

    #TODO also account for cut_min?
    control_variable <- calc_water_deviations(
      files_scenario = files_scenario,
      files_reference = files_reference,
      spatial_scale = spatial_scale,
      time_span_scenario = time_span_scenario,
      time_span_reference =  time_span_reference,
      method = method,
      time_aggregation_args = time_aggregation_args,
      config_args = config_args,
      thresholds = thresholds,
      variable = "discharge"
    )
  }

  return(control_variable)
}


calc_bluewater_efrs <- function(
  files_scenario,
  files_reference,
  spatial_scale,
  time_span_scenario = time_span_scenario,
  time_span_reference = NULL,
  method = "gerten2020",
  cut_min = 0.0864,
  time_aggregation_args = list(),
  config_args = list(),
  thresholds = NULL
) {

  if (spatial_scale != "grid") {
    stop("Method \"gerten2020\" is only applicable for spatial_scale = \"grid\"") #nolint
  }

  # reference discharge ---------------------------------------------------- #
  discharge_reference %<-% read_io_format(
    files_reference$discharge,
    time_span_reference,
    aggregate = list(band = sum),
    spatial_subset = config_args$spatial_subset
  )

  # scenario discharge ----------------------------------------------------- #
  discharge_scenario %<-% read_io_format(
    files_scenario$discharge,
    time_span_scenario,
    aggregate = list(band = sum),
    spatial_subset = config_args$spatial_subset
  )

  # ------------------------------------------------------------------------ #

  # TODO understand what this does
  if (length(time_span_reference) < length(time_span_scenario)) {
    nyear_ref <- length(time_span_scenario)
  } else {
    nyear_ref <- NULL
  }

  # average discharge reference
  avg_discharge_reference %<-% do.call(
    aggregate_time,
    append(list(x = discharge_reference,
                nyear_reference = nyear_ref),
           time_aggregation_args)
  )

  # average discharge scenario
  avg_discharge_scenario %<-% do.call(
    aggregate_time,
    append(list(x = discharge_scenario),
           time_aggregation_args)
  )


  # calc efrs for vmf_min and vmf_max
  efr_uncertain %<-% calc_efrs(discharge_reference,
                                "vmf_min",
                                time_aggregation_args)
  efr_safe %<-% calc_efrs(discharge_reference,
                          "vmf_max",
                          time_aggregation_args)

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

  # check if vector was returned (loss of dimnames) -> reconvert to array
  #TODO check if necessary
  if (is.null(dim(control_variable))) {
    control_variable <- array(
      control_variable,
      dim = c(cell = dim(status_frac_monthly)[["cell"]], 1),
      dimnames = list(cell = dimnames(status_frac_monthly)[["cell"]], 1)
    )
  }
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
    "environmental flow requirements"
  attr(control_variable, "thresholds") <- thresholds
  attr(control_variable, "spatial scale") <- spatial_scale

  class(control_variable) <- c("control_variable")
  return(control_variable)
}