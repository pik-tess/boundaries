# calculate irrigation mask for basins based on drainage input and irrigation
#   output to exclude non irrigated basins
#   also works with aggregate_time
#   returns lpjml output with either just cell or + year/window
calc_irrigation_mask <- function(
  files_scenario,
  time_span,
  time_aggregation_args = list(),
  config_args = list()
) {

  # -------------------------------------------------------------------------- #
  irrigation_scenario <- lpjmlkit::read_io(
    files_scenario$irrig, subset = list(year = time_span)
  ) %>%
    lpjmlkit::transform(to = c("year_month_day")) %>%
    lpjmlkit::as_array(aggregate = list(month = sum, band = sum)) %>%
    suppressWarnings()
  # -------------------------------------------------------------------------- #

  # average irrigation
  avg_irrigation_scenario <- do.call(
    aggregate_time,
    append(
      list(x = irrigation_scenario),
      time_aggregation_args
    )
  )

  # ------------------------------------------------------------------------- #
  cellinfo <- indexing_drainage(drainage_file = files_scenario$drainage)
  if (!is.null(config_args$spatial_subset)) {
    cellinfo <- lpjmlkit::asub(cellinfo, cell = config_args$spatial_subset)
  }
  endcell <- lpjmlkit::asub(cellinfo, band = "endcell")

  basin_ids <- sort(unique(endcell))
  irrmask_basin <- array(
    0,
    dim = dim(avg_irrigation_scenario)[names(dim(avg_irrigation_scenario)) != "month"], # nolint
    dimnames = dimnames(avg_irrigation_scenario)[names(dimnames(avg_irrigation_scenario)) != "month"] # nolint
  )

  third_dim <- names(dim(avg_irrigation_scenario))[
    !names(dim(avg_irrigation_scenario)) %in% c("cell", "month")
  ] %>% {
    if (rlang::is_empty(.)) NULL else .
  }

  for (id in seq_len(length(basin_ids))) {
    basincell <- which(endcell == basin_ids[id])
    if (is.null(third_dim) | length(dim(drop(avg_irrigation_scenario))) < 3) { # nolint
      check_gt0 <- sum(
        lpjmlkit::asub(avg_irrigation_scenario, cell = basincell)
      )
    } else {
      check_gt0 <- apply(
        lpjmlkit::asub(avg_irrigation_scenario, cell = basincell),
        third_dim,
        sum
      )
    }
    basin_replace <- lpjmlkit::asub(irrmask_basin, cell = basincell) %>%
      `[<-`(check_gt0 > 0, value = 1)

    lpjmlkit:::asub(x = irrmask_basin, cell = basincell) <- basin_replace
  }
  return(irrmask_basin)
}

# show the full downstream drainage route for cell ind until final drainage
#   (ocean or inland sink)
show_route <- function(ind, routing_table) {
  if (routing_table[ind] < 1) {
    # can be 0 or -8 -> endcell or nacell
    return(ind)
  } else {
    return(c(ind,
             show_route(routing_table[ind], routing_table)))
  }
}


indexing_drainage <- function(drainage_file) {

  drainage <- tryCatch(
    lpjmlkit::read_io(drainage_file)$data %>%
      suppressWarnings() %>%
      drop(),
    error = function(e) {
      lpjmlkit::read_io(
        drainage_file,
        datatype = 2
      )$data %>%
        suppressWarnings() %>%
        drop()
    }
  )
  # -------------------------------------------------------------------------- #

  # add 1 since in C indexing starts at 0 but in R at 1
  routing <- drainage[, 1] + 1
  ncell <- length(routing)
  cellinfo <- array(0, dim = c(ncell, 3))
  cellinfo[, 3] <- routing
  for (cell in 1:ncell) {
    route <- show_route(cell, routing)
    # (over)writing all cells from route
    cellinfo[route, 1] <- seq(length(route), 1, -1) # rank (former cellindex)
    cellinfo[cell, 2] <- route[length(route)] # endcell
  }
  dimnames(cellinfo) <- list(cell = 0:(ncell - 1),
                             band = c("rank", "endcell", "routing"))
  return(cellinfo)
}
