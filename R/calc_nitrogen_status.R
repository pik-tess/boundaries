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
#' @param pb_thresholds list with upper and lower threshold for N concentration
#' (mg N/l) in runoff to surface water
#' Default: upper = 2.5, lower = 1
#' (based on de Vries et al. 2013, https://doi.org/10.1016/j.cosust.2013.07.004)
#' Alternative: upper = 5, lower = 2
#' (based on Schulte-Uebbing et al. 2022,
#' https://doi.org/10.1038/s41586-022-05158-2:
#' "we used a threshold for N concentration in run-off to surface water. This
#' threshold was set to 5.0 mgN/l, based on the assumption that on average 50%
#' of N entering surface water is removed through retention and sedimentation"))
#'
#' @param prefix_monthly_output character. Provide a prefix if required for
#' monthly LPJmL output files, e.g. `"m"` for `"mdischarge.bin"` instead of
#' `"discharge.bin"` (default is `""`)
#'
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#' time series analysis
#' 
#' @param read_args list of arguments for reading output.
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
                                 input_files = list(prec = "/p/projects/lpjml/input/historical/CRUDATA_TS3_23/gpcc_v7_cruts3_23_precip_1901_2013.clm"), #nolint
                                 method = "braun2022",
                                 cut_arid = 0.2,
                                 cut_runoff = 0,
                                 with_groundwater_denit = TRUE,
                                 pb_thresholds = list(lower = 1, upper = 2.5),
                                 prefix_monthly_output = "",
                                 avg_nyear_args = list(),
                                 # to be replaced by lpjmlKit::read_output
                                 read_args = list(start_year = 1901,
                                                  headersize = 0)
                                 ) {
  # verify available methods
  method <- match.arg(method, c("braun2022",
                                "braun2022_minusref"))

  if (.Platform$OS.type == "windows") {
    future_plan <- future::plan("multisession")
  } else {
    future_plan <- future::plan("multicore")
  }
  on.exit(future::plan(future_plan))

  # sub function to be used for scenario and reference run (braun2022_minusref)
  calc_nitrogen_leach <- function(path_data,
                                  time_span,
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
    cell_area <- lpjmliotools::cellarea #lpjmlKit::calc_cellarea(
      #lpjmlKit::subset_array(lpjml_grid, list(coordinate = "lat"))
    #)

    # TO BE REPLACED BY lpjmlKit::read_output -------------------------------- #
    #   hardcoded values to be internally replaced
    # read runoff
    runoff <- tmp_read_monthly(
      file_name = paste0(path_data, "/", prefix_monthly_output, "runoff.bin"),
      time_span = time_span,
      start_year = read_args$start_year,
      nstep = 12,
      ncell = 67420,
      nbands = 1,
      size = 4,
      headersize = read_args$headersize
    )
    # ------------------------------------------------------------------------ #

    # TO BE REPLACED BY lpjmlKit::read_output -------------------------------- #
    #   hardcoded values to be internally replaced
    # read leaching
    leaching <- tmp_read_monthly(
      file_name = paste0(path_data, "/", prefix_monthly_output, "leaching.bin"),
      time_span = time_span,
      start_year = read_args$start_year,
      nstep = 12,
      ncell = 67420,
      nbands = 1,
      size = 4,
      headersize = read_args$headersize
    )
    # ------------------------------------------------------------------------ #

    # average runoff
    monthly_runoff <- do.call(average_nyear_window,
                          append(list(x = runoff),
                                 avg_nyear_args))

    # average leaching
    monthly_leaching <- do.call(average_nyear_window,
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


    # get name of third dimension of array (band, year, ...)
    third_dim <- names(dim(monthly_runoff))[
      !names(dim(monthly_runoff)) %in% c("cell", "month")
    ] %>% {
      if (rlang::is_empty(.)) NULL else .
    }
    third_dim <<- third_dim

    avg_runoff <- apply(
      monthly_runoff,
      names(dim(monthly_runoff))[
        names(dim(monthly_runoff)) %in% c("cell", third_dim)
      ],
      sum,
      na.rm = TRUE)
    avg_leaching <- apply(
      monthly_leaching,
      names(dim(monthly_leaching))[
        names(dim(monthly_leaching)) %in% c("cell", third_dim)
      ],
      sum,
      na.rm = TRUE)

    # check if vector was returned (loss if dimnames) -> reconvert to array
    if (is.null(dim(avg_runoff))) {
      avg_runoff <- array(
        avg_runoff,
        dim = c(cell = dim(monthly_runoff)[["cell"]], 1),
        dimnames = list(cell = dimnames(monthly_runoff)[["cell"]], 1)
      )
      avg_leaching <- array(
        avg_leaching,
        dim = c(cell = dim(monthly_leaching)[["cell"]], 1),
        dimnames = list(cell = dimnames(monthly_leaching)[["cell"]], 1)
      )
    }

    status_frac <- ifelse(avg_runoff > 0,
                           (avg_leaching * cell_area * 1e3 * loss_factor) /
                           (avg_runoff * cell_area), 0)

    return(status_frac)
  }

  # apply defined method
  switch(method,
    braun2022 = {
      status_frac <- calc_nitrogen_leach(
        path_data = path_scenario,
        time_span = time_span_scenario,
        with_groundwater_denit = with_groundwater_denit,
        prefix_monthly_output = prefix_monthly_output,
        avg_nyear_args = avg_nyear_args,
        # to be replaced by lpjmlKit::read_output
        start_year = read_args$start_year)
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
      status_frac_scenario <- calc_nitrogen_leach(
        path_data = path_scenario,
        time_span = time_span_scenario,
        with_groundwater_denit = with_groundwater_denit,
        prefix_monthly_output = prefix_monthly_output,
        avg_nyear_args = avg_nyear_args,
        # to be replaced by lpjmlKit::read_output
        start_year = read_args$start_year)

      if (!is.null(nyear_ref)) {
        avg_nyear_args["nyear_reference"] <- nyear_ref
      }

      # calculate leaching concentration and loss rate for reference output
      status_frac_reference <- calc_nitrogen_leach(
        path_data = path_reference,
        time_span = time_span_reference,
        with_groundwater_denit = with_groundwater_denit,
        prefix_monthly_output = prefix_monthly_output,
        avg_nyear_args = avg_nyear_args,
        # to be replaced by lpjmlKit::read_output
        start_year = read_args$start_year)

      # subtract scenario leaching concentration and loss rate from reference
      status_frac <- status_frac_scenario - status_frac_reference
      status_frac[status_frac < 0] <- 0
    }
  )

  # TO BE REPLACED BY lpjmlKit::read_output -------------------------------- #
  #   hardcoded values to be internally replaced
  # read potential evapotranspiration
  pet <- tmp_read_monthly(
    file_name = paste0(path_scenario, "/", prefix_monthly_output, "pet.bin"),
    time_span = time_span_scenario,
    start_year = read_args$start_year,
    nstep = 12,
    ncell = 67420,
    nbands = 1,
    size = 4,
    headersize = read_args$headersize
  )
  # ------------------------------------------------------------------------ #

  # TO BE REPLACED BY lpjmlKit::read_output -------------------------------- #
  #   hardcoded values to be internally replaced
  # read precipitation
  if (is.null(input_files$prec)) {
    prec <- tmp_read_monthly(
      file_name = paste0(path_scenario, "/", prefix_monthly_output, "prec.bin"),
      time_span = time_span_scenario,
      start_year = read_args$start_year,
      nstep = 12,
      ncell = 67420,
      nbands = 1,
      size = 4,
      headersize = read_args$headersize
  )
  } else {
    prec <- lpjmliotools::autoReadInput(inFile = input_files$prec,
                                       getyearstart = time_span_scenario[1],
                                       getyearstop = time_span_scenario[2]) %>%
            rename_step2month()
    prec <- aperm(prec, c(2, 1, 3))
    dim_names <- list("cell" = as.character(as.numeric(dimnames(prec)[["cell"]]) + 1),
                            "month" = dimnames(prec)[["band"]],
                            "year" = dimnames(prec)[["year"]])
    prec <- array(prec,
                  dim = setNames(sapply(dim_names, length), names(dim_names)),
                  dimnames = dim_names)

    # TODO the actual problem is that the "month" dimension is called "band"
  }
  # ------------------------------------------------------------------------ #
  # average pet
  avg_pet <- do.call(average_nyear_window,
                     append(list(x = pet),
                            avg_nyear_args))

  # average precipitation
  avg_prec <- do.call(average_nyear_window,
                     append(list(x = prec),
                            avg_nyear_args))

  # calculate annual potential evapotranspiration mean
  avg_pet_annual <- apply(
    avg_pet,
    names(dim(avg_pet))[
      names(dim(avg_pet)) %in% c("cell", third_dim)
    ],
    mean,
    na.rm = TRUE)

  # calculate annual precipiation mean
  avg_prec_annual <- apply(
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

  # ------------------------------------------------------------------------
  # TO BE REPLACED BY lpjmlKit::read_output -------------------------------- #
  # hardcoded values to be internally replaced
  # read runoff
  runoff <- tmp_read_monthly(
      file_name = paste0(path_scenario, "/", prefix_monthly_output, "runoff.bin"),
      time_span = time_span_scenario,
      start_year = read_args$start_year,
      nstep = 12,
      ncell = 67420,
      nbands = 1,
      size = 4,
      headersize = read_args$headersize
  )

  # average runoff
  runoff_monthly <- do.call(average_nyear_window,
                          append(list(x = runoff),
                                 avg_nyear_args))

  # calculate annual runoff
  runoff_annual <- apply(
    runoff_monthly,
    names(dim(runoff_monthly))[
      names(dim(runoff_monthly)) %in% c("cell", third_dim)
    ],
    sum,
    na.rm = TRUE)

  # to display arid cells (leaching behaviour threshold) in other color (grey):
  cells_arid <- array(FALSE,
                      dim = dim(status_frac),
                      dimnames = dimnames(status_frac))
  cells_arid[which(global_aridity_index <= cut_arid)] <- TRUE

  # to display cells with low runoff in other color (grey):
  cells_low_runoff <- array(FALSE,
                            dim = dim(status_frac),
                            dimnames = dimnames(status_frac))
  cells_low_runoff[which(runoff_annual <= cut_runoff)] <- TRUE

  # init pb_status based on status_frac
  pb_status <- status_frac
  # high risk
  pb_status[status_frac >= pb_thresholds$upper] <- 3
  # increasing risk
  pb_status[status_frac < pb_thresholds$upper &
            status_frac >= pb_thresholds$lower] <- 2
  # safe zone
  pb_status[status_frac < pb_thresholds$lower] <- 1
  # non applicable cells
  pb_status[cells_arid] <- 0
  pb_status[cells_low_runoff] <- 0

  return(pb_status)
}


# TO BE REPLACED BY lpjmlKit::read_output ------------------------------------ #
tmp_read_monthly <- function(file_name,
                             time_span,
                             start_year,
                             nstep,
                             ncell,
                             nbands,
                             size,
                             headersize) {
  # scenario runofv
  file_con <- file(file_name,
                   "rb")
  seek(file_con,
       where = headersize + (time_span[1] - start_year) *
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
