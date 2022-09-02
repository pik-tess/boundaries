#' Calculate the planetary boundary status for the nitrogen boundary
#'
#' Calculate the PB status for the nitrogen boundary based on a scenario LPJmL
#' run and if `method == "braun2022_minusref"` a reference LPJmL run.
#'
#' @param path_scenario output directory (character string) of the scenario
#' LPJmL run where binary files (soon with metafiles) are written
#'
#' @param path_reference output directory (character string) of the reference
#' - if used - `method == "braun2022_minusref`
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#' as an integer vector, e.g. `1982:2011` (default)
#'
#' @param time_span_reference time span to be used if scenario run used, defined
#' as an integer vector, e.g. `1901:1930`. Can differ in offset and length from
#' `time_span_scenario`! If `NULL` value of `time_span_scenario` is used
#'
#' @param method method (character string) to be used , currently available
#' method is `"braun2022"` based on unpublished suggestion by Johanna Braun.
#' Second method option is `"braun2022_minusref"` to subtract reference run
#' output
#'
#' @param temporal_resolution character. Temporal resolution, available options
#' are `"annual"` (default) and `"monthly"`
#'
#' @param cut_min double. Exclude boundary calculations for Q < cut_min
#'
#' @param with_groundwater_storage logical. Include global assumptions based on
#' assumed water storages (not included
#' in LPJmL). Defaults to FALSE
#'
#' @param with_groundwater_denit logical. Include global assumptions made on
#' groundwater denitrification losses. Defaults to TRUE
#'
#' @param prefix_monthly_output character. Provide a prefix if required for
#' monthly LPJmL output files, e.g. `"m"` for `"mdischarge.bin"` instead of
#' `"discharge.bin"` (default is `""`)
#'
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#' time series analysis
#'
#' @examples
#' \dontrun{
#'  calc_nitrogen_status(path_scenario, path_reference)
#' }
#'
#' @md
#' @importFrom future %<-%
#' @export
calc_nitrogen_status <- function(path_scenario,
                                 path_reference,
                                 time_span_scenario = c(1982, 2011),
                                 time_span_reference = NULL,
                                 method = "braun2022",
                                 temporal_resolution = "annual",
                                 # global aridity index < than 0.2 (arid)
                                 cut_min = 0.2,
                                 with_groundwater_storage = FALSE,
                                 with_groundwater_denit = TRUE,
                                 prefix_monthly_output = "",
                                 avg_nyear_args = list(),
                                 # to be replaced by lpjmlKit::read_output
                                 start_year = 1901) {
  # verify available methods
  method <- match.arg(method, c("braun2022",
                                "braun2022_minusref"))
  # verify available temporal resolution
  temporal_resolution <- match.arg(temporal_resolution, c("annual",
                                                          "monthly"))

  if (.Platform$OS.type == "windows") {
    future_plan <- future::plan("multisession")
  } else {
    future_plan <- future::plan("multicore")
  }
  on.exit(future::plan(future_plan))

  # sub function to be used for scenario and reference run (braun2022_minusref)
  calc_nitrogen_leach <- function(path_data,
                                  time_span,
                                  temporal_resolution,
                                  cut_min,
                                  with_groundwater_storage,
                                  with_groundwater_denit,
                                  prefix_monthly_output,
                                  avg_nyear_args,
                                  # to be replaced by lpjmlKit::read_output
                                  start_year) {
    # TO BE REPLACED BY lpjmlKit::read_output -------------------------------- #
    # read grid
    ncell <- 67420
    size <- 2
    grid_file <- file(paste(path_data, "grid.bin", sep = "/"), "rb")
    lpjml_grid <- readBin(grid_file, integer(), n = 2 * ncell, size = size) /
                  100
    close(grid_file)
    dim(lpjml_grid) <- c(coordinate = 2, cell = ncell)
    dimnames(lpjml_grid) <- list(coordinate = c("lon", "lat"),
                                 cell = seq_len(ncell))
    # ------------------------------------------------------------------------ #
    cell_area <-  lpjmlKit::calc_cellarea(
      lpjmlKit::subset_array(lpjml_grid, list(coordinate = "lat"))
    )

    # TO BE REPLACED BY lpjmlKit::read_output -------------------------------- #
    #   hardcoded values to be internally replaced
    # read runoff
    runoff %<-% tmp_read_monthly(
      file_name = paste0(path_data, "/", prefix_monthly_output, "runoff.bin"),
      time_span = time_span_scenario,
      start_year = start_year,
      nstep = 12,
      ncell = 67420,
      nbands = 1,
      size = 4
    )
    # ------------------------------------------------------------------------ #

    # TO BE REPLACED BY lpjmlKit::read_output -------------------------------- #
    #   hardcoded values to be internally replaced
    # read runoff
    leaching %<-% tmp_read_monthly(
      file_name = paste0(path_data, "/", prefix_monthly_output, "leaching.bin"),
      time_span = time_span_scenario,
      start_year = start_year,
      nstep = 12,
      ncell = 67420,
      nbands = 1,
      size = 4
    )
    # ------------------------------------------------------------------------ #

    # average runoff
    avg_runoff %<-% do.call(average_nyear_window,
                          append(list(x = runoff),
                                 avg_nyear_args))

    # average runoff
    avg_leaching %<-% do.call(average_nyear_window,
                            append(list(x = leaching),
                                   avg_nyear_args))

    # temporary solution using shares of global losses after Bouwman et al 2013
    #   (figure 4) https://doi.org/10.1098/rstb.2013.0112

    if (with_groundwater_storage) {
      # net_flow = net soil N flow + net groundwater N flow (inflow - outflow)
      net_flow <- 65 + 15 - 6
    } else {
      # net_flow = net soil N flow
      net_flow <- 65
    }

    if (with_groundwater_denit) {
      # gross_flow = gross soil & groundwater N flow + gross surface N flow
      gross_flow <- 93 + 11
    } else {
      # gross_flow = gross soil N flow + gross groundwater N flow +
      #   gross riparian N flow + gross surface N flow
      gross_flow <- 49 + 5 + 15 + 11
    }

    loss_factor <- net_flow / gross_flow

    status_frac_monthly <- ifelse(avg_runoff > 0,
                           (avg_leaching * cell_area * 1e3 * loss_factor) /
                           (avg_runoff * cell_area), 0)

    # get name of third dimension of array (band, year, ...)
    third_dim <- names(dim(status_frac_monthly))[
      !names(dim(status_frac_monthly)) %in% c("cell", "month")
    ] %>% {
      if (rlang::is_empty(.)) NULL else .
    }
    third_dim <<- third_dim

    if (temporal_resolution == "annual") {
      status_frac <- apply(
        status_frac_monthly,
        names(dim(status_frac_monthly))[
          names(dim(status_frac_monthly)) %in% c("cell", third_dim)
        ],
        mean,
        na.rm = TRUE)

      # check if vector was returned (loss if dimnames) -> reconvert to array
      if (is.null(dim(status_frac))) {
        status_frac <- array(
          status_frac,
          dim = c(cell = dim(status_frac_monthly)[["cell"]], 1),
          dimnames = list(cell = dimnames(status_frac_monthly)[["cell"]], 1)
        )
      }
    } else {
      status_frac <- status_frac_monthly
    }
    return(status_frac)
  }

  # apply defined method
  switch(method,
    braun2022 = {
      status_frac <- calc_nitrogen_leach(
        path_data = path_scenario,
        time_span = time_span_scenario,
        temporal_resolution = temporal_resolution,
        cut_min = cut_min,
        with_groundwater_storage = with_groundwater_storage,
        with_groundwater_denit = with_groundwater_denit,
        prefix_monthly_output = prefix_monthly_output,
        avg_nyear_args = avg_nyear_args,
        # to be replaced by lpjmlKit::read_output
        start_year = start_year)
    },
    braun2022_minusref = {
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
      # calculate leaching concentration and loss rate for scenario output
      status_frac_scenario %<-% calc_nitrogen_leach(
        path_data = path_scenario,
        time_span = time_span_scenario,
        temporal_resolution = temporal_resolution,
        cut_min = cut_min,
        with_groundwater_storage = with_groundwater_storage,
        with_groundwater_denit = with_groundwater_denit,
        prefix_monthly_output = prefix_monthly_output,
        avg_nyear_args = avg_nyear_args,
        # to be replaced by lpjmlKit::read_output
        start_year = start_year)

      if (!is.null(nyear_ref)) {
        avg_nyear_args["nyear_reference"] <- nyear_ref
      }

      # calculate leaching concentration and loss rate for reference output
      status_frac_reference %<-% calc_nitrogen_leach(
        path_data = path_reference,
        time_span = time_span_reference,
        temporal_resolution = temporal_resolution,
        cut_min = cut_min,
        with_groundwater_storage = with_groundwater_storage,
        with_groundwater_denit = with_groundwater_denit,
        prefix_monthly_output = prefix_monthly_output,
        avg_nyear_args = avg_nyear_args,
        # to be replaced by lpjmlKit::read_output
        start_year = start_year)

      # subtract scenario leaching concentration and loss rate from reference
      status_frac <- status_frac_scenario - status_frac_reference
      status_frac[status_frac < 0] <- 0
    }
  )

  # TO BE REPLACED BY lpjmlKit::read_output -------------------------------- #
  #   hardcoded values to be internally replaced
  # read runoff
  pet %<-% tmp_read_monthly(
    file_name = paste0(path_scenario, "/", prefix_monthly_output, "pet.bin"),
    time_span = time_span_scenario,
    start_year = start_year,
    nstep = 12,
    ncell = 67420,
    nbands = 1,
    size = 4
  )
  # ------------------------------------------------------------------------ #

  # TO BE REPLACED BY lpjmlKit::read_output -------------------------------- #
  #   hardcoded values to be internally replaced
  # read runoff
  prec %<-% tmp_read_monthly(
    file_name = paste0(path_scenario, "/", prefix_monthly_output, "prec.bin"),
    time_span = time_span_scenario,
    start_year = start_year,
    nstep = 12,
    ncell = 67420,
    nbands = 1,
    size = 4
  )
  # ------------------------------------------------------------------------ #
  # average runoff
  avg_pet %<-% do.call(average_nyear_window,
                     append(list(x = pet),
                            avg_nyear_args))

  # average runoff
  avg_prec %<-% do.call(average_nyear_window,
                     append(list(x = prec),
                            avg_nyear_args))

  # calculate annual potential evapotranspiration mean
  avg_pet_annual %<-% apply(
    avg_pet,
    names(dim(avg_pet))[
      names(dim(avg_pet)) %in% c("cell", third_dim)
    ],
    mean,
    na.rm = TRUE)

  # calculate annual precipiation mean
  avg_prec_annual %<-% apply(
    avg_prec,
    names(dim(avg_prec))[
      names(dim(avg_prec)) %in% c("cell", third_dim)
    ],
    mean,
    na.rm = TRUE)

  # calculate global aridity index (AI) as an indicator for a level under which
  #   the calculation of leaching just cannot show realistic behavior, see also
  #   on the AI: https://doi.org/10.6084/m9.figshare.7504448.v4%C2%A0
  #     & (first descr.) https://wedocs.unep.org/xmlui/handle/20.500.11822/30300
  #   on nitrogen processes in arid areas: https://www.jstor.org/stable/45128683
  #     -> indicates boundary to "arid" as thresholds (=< 0.2)
  #   on "arid threshold" (indirectly): https://doi.org/10.1038/ncomms5799
  #     -> threshold for behaviour change in nitrogen cycling (=< 0.32)
  global_aridity_index <- avg_prec_annual / avg_pet_annual + 1e-9

  # to display arid cells (leaching behaviour threshold) in other color (grey):
  cells_arid <- array(FALSE,
                      dim = dim(status_frac),
                      dimnames = dimnames(status_frac))
  cells_arid[which(global_aridity_index <= cut_min)] <- TRUE

  # init pb_status based on status_frac
  pb_status <- status_frac
  # high risk
  pb_status[status_frac >= 2.5] <- 3
  # increasing risk
  pb_status[status_frac < 2.5 & status_frac >= 1] <- 2
  # safe zone
  pb_status[status_frac < 1] <- 1
  # non applicable cells
  pb_status[cells_arid] <- 0

  return(pb_status)
}


# TO BE REPLACED BY lpjmlKit::read_output ------------------------------------ #
tmp_read_monthly <- function(file_name,
                             time_span,
                             start_year,
                             nstep,
                             ncell,
                             nbands,
                             size) {
  # scenario runofv
  file_con <- file(file_name,
                   "rb")
  seek(file_con,
       where = (time_span[1] - start_year) *
               nstep * nbands * ncell * size,
       origin = "start")
  lpjml_data <- readBin(file_con,
                        double(),
                         n = (ncell * nstep * nbands *
                              (time_span[2] -
                               time_span[1] + 1)),
                        size = size)
  close(file_con)
  dim(lpjml_data) <- c(cell = ncell,
                       month = nstep,
                       year = (time_span[2] -
                                time_span[1] + 1))
  dimnames(lpjml_data) <- list(cell = seq_len(ncell),
                               month = seq_len(nstep),
                               year = seq(time_span[1],
                                           time_span[2]))
  return(lpjml_data)
}