# calculate irrigation mask for basins based on drainage input and irrigation
#   output to exclude non irrigated basins
#   also works with average_nyear_window
#   returns lpjml output with either just cell or + year/window
calc_irrigation_mask <- function(path_output,
                                 time_span,
                                 prefix_monthly_output = "",
                                 avg_nyear_args=list(),
                                 start_year = 1901,
                                 path_input = (
                                  paste0("/p/projects/lpjml/input/historical",
                                         "/input_VERSION2/")
                                 )) {

  # TO BE REPLACED BY lpjmlKit::read_output ---------------------------------- #
  #   hardcoded values to be internally replaced

  nstep <- 12
  nbands <- 1
  ncell <- 67420
  size <- 4

  s_path <- file(paste0(path_output, "/", prefix_monthly_output, "irrig.bin"),
                 "rb")
  seek(s_path,
       where = (time_span[1] - start_year) *
               nstep * nbands * ncell * size,
       origin = "start")
  irrigation_scenario <- readBin(s_path,
                                 double(),
                                  n = (ncell * nstep * nbands *
                                       (time_span[2] - time_span[1] + 1)),
                                 size = size)
  close(s_path)
  dim(irrigation_scenario) <- c(cell = ncell,
                                month = nstep,
                                year = (time_span[2] - time_span[1]) + 1)
  dimnames(irrigation_scenario) <- list(cell = seq_len(ncell),
                                        month = seq_len(nstep),
                                        year = seq(time_span[1],
                                                    time_span[2]))
  # -------------------------------------------------------------------------- #

  # average discharge reference
  avg_irrigation_scenario %<-% do.call(average_nyear_window,
                                       append(list(x = irrigation_scenario),
                                              avg_nyear_args))

  # TO BE REPLACED BY lpjmlKit::read_input ----------------------------------- #
  #   hardcoded values to be internally replaced
  input_data_size <- 4
  header <- suppressWarnings(lpjmlKit::read_header(
    filename = paste(path_input, "drainage.bin", sep = "/")
  ))
  headersize <- lpjmlKit::get_headersize(header)
  input_file <- file(paste(path_input, "drainage.bin", sep = "/"), "rb")
  seek(input_file, where = headersize, origin = "start")
  drainage <- readBin(input_file,
                    integer(),
                    n = lpjmlKit::get_header_item(header, "nyear") *
                        lpjmlKit::get_header_item(header, "ncell") *
                        lpjmlKit::get_header_item(header, "nbands"),
                  size = input_data_size) *
              lpjmlKit::get_header_item(header, "scalar")
  close(input_file)      #remove to save space
  dim(drainage) <- c(band = lpjmlKit::get_header_item(header, "nbands"),
                     cell = lpjmlKit::get_header_item(header, "ncell"))
  dimnames(drainage) <- list(
    band = seq_len(lpjmlKit::get_header_item(header, "nbands")),
    cell = seq_len(lpjmlKit::get_header_item(header, "ncell")))
  # -------------------------------------------------------------------------- #

  # add 1 since in C indexing starts at 0 but in R at 1
  routing <- drainage[1, ] + 1
  endcell <- array(0, dim = ncell)
  cellindex <- endcell
  for (cell in 1:ncell) {
    route <- show_route(cell, routing)
    cellindex[route] <- seq(length(route), 1, -1)
    endcell[cell] <- route[length(route)]
  }

  basin_ids <- sort(unique(endcell))
  irrmask_basin <- array(0,
                         dim = dim(avg_irrigation_scenario)[
                           names(dim(avg_irrigation_scenario)) != "month"
                         ],
                         dimnames = dimnames(avg_irrigation_scenario)[
                           names(dimnames(avg_irrigation_scenario)) != "month"
                         ])

  third_dim <- names(dim(avg_irrigation_scenario))[
    !names(dim(avg_irrigation_scenario)) %in% c("cell", "month")
  ] %>% {
    if (rlang::is_empty(.)) NULL else .
  }

  for (id in seq_len(length(basin_ids))) {
    basincell <- which(endcell == basin_ids[id])
    if (is.null(third_dim) | length(dim(drop(avg_irrigation_scenario))) < 3) {
      check_gt0 <- sum(
        lpjmlKit::subset_array(avg_irrigation_scenario,
                               list("cell" = basincell))
      )
    } else {
      check_gt0 <- apply(
        lpjmlKit::subset_array(avg_irrigation_scenario,
                               list(cell = basincell)),
        third_dim,
        sum
      )
    }
    basin_replace <- lpjmlKit::subset_array(irrmask_basin,
                                            list(cell = basincell)) %>%
      `[<-`(check_gt0 > 0, value = 1)

    irrmask_basin <- lpjmlKit::replace_array(
      x = irrmask_basin,
      subset_list = list(cell = basincell),
      y = basin_replace
    )
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
