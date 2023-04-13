#' Calculate the planetary boundary status for the greenwater boundary
#'
#' Calculate the PB status for the greenwater (former freshwater) boundary based
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
#' as an integer vector, e.g. `1982:2011` (default)
#'
#' @param time_span_reference time span to be used for the scenario run, defined
#' as an integer vector, e.g. `1901:1930`. Can differ in offset and length from
#' `time_span_scenario`! If `NULL` value of `time_span_scenario` is used
#'
#' @param method method (character string) to be used , currently available
#' method is `c("wang-erlandsson2022")` based on
#' [Wang-Erlandsson et al. 2022](https://doi.org/10.1038/s43017-022-00287-8).
#'
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#' time series analysis
#'
#' @examples
#' \dontrun{
#'  calc_greenwater_status(path_scenario, path_reference)
#' }
#'
#' @md
#' @export
calc_greenwater_status_global <- function(files_scenario,
                                   files_reference,
                                   time_span_scenario = as.character(1982:2011),
                                   time_span_reference = NULL,
                                   method = "wang-erlandsson2022",
                                   avg_nyear_args = list()
                                   ) {
  #TODO add second method: not only referring to driest/wettest month, but
  # to the whole year
  
  # verify available methods
  method <- match.arg(method, c("wang-erlandsson2022"))

  # check time_spans of scenario and reference runs
  if (is.null(time_span_reference)) {
    time_span_reference <- time_span_scenario
    nyear_ref <- NULL
  } else {
    if (length(time_span_reference) > length(time_span_scenario)) {
      stop(paste0("time_span_reference is longer than time_span_scenario.",
                  "Define a time_span_reference that is shorter than",
                  "time_span_scenario"))
    } else if (length(time_span_reference) < length(time_span_scenario)) {
      nyear_ref <- length(time_span_scenario)
    } else {
      nyear_ref <- NULL
    }
  }

  # -------------------------------------------------------------------------- #
  # reference rootmoisture
  rootmoist_reference <- lpjmlkit::read_io(
      files_reference$rootmoist, subset = list(year = time_span_reference)
      ) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      lpjmlkit::as_array(aggregate = list(band = sum)) %>%
      suppressWarnings()

  # scenario rootmoisture
  rootmoist_scenario <- lpjmlkit::read_io(
      files_scenario$rootmoist, subset = list(year = time_span_scenario)
      ) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      lpjmlkit::as_array(aggregate = list(band = sum)) %>%
      suppressWarnings()

  # grid for cellarea
  grid <- lpjmlkit::read_io(
      files_scenario$grid,
      silent = TRUE
      )
  cell_area <- lpjmlkit::calc_cellarea(grid)

  # -------------------------------------------------------------------------- #
  #calculate the green water 5% and 95% quantiles of the baseline period
  quants <- calc_water_baseline(rootmoist_reference)

  # calculate GW dry & wet departures and return percentage
  # of annual area of departure
  area_with_departure <- calc_water_depart_area(rootmoist_scenario, quants)
  #TODO translate into global PB status
  return(area_with_departure)
}



# calculate GW dry & wet departures and return mean annual area of departure
calc_water_depart_area <- function(file_scenario, quants) {

  q5_base <- unlist(quants[1])
  q95_base <- unlist(quants[2])

  result <- list()
  # for each hist year calc area of rootmoist wet & dry departures
  # identify cells with dry departures
  # driest month per gridcell for year i
  dry <- apply(file_scenario, c(1, 3), min)
  dry[dry >= q5_base] <- NA
  dry[dry >= 0] <- 1
  result$dry <- mean(
                    apply((dry * cell_area) / sum(cell_area) * 100, 2, sum,
                           na.rm = TRUE)
                )

  #identify cells with wet departures
  #wettest month per gridcell for year i -> ignores which month
  wet       <- apply(file_scenario, c(1, 3), max)
  wet[wet <= q95_base] <- NA
  wet[wet >= 0] <- 1
  result$wet <- mean(
                    apply((wet * cell_area) / sum(cell_area) * 100, 2, sum,
                           na.rm = TRUE)
                )
  return(result)
}