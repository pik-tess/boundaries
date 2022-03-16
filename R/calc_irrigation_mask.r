# calculate irrigation mask for basins based on drainage input and irrigation
#   output to exclude non irrigated basins
#   also works with average_nyear_window
#   returns lpjml output with either just cells or + years/windows
calc_irrigation_mask <- function(path_output,
                                 time_span,
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

  s_path <- file(paste(path_output, "irrig.bin", sep = "/"), "rb")
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
  dim(irrigation_scenario) <- c(cells = ncell,
                                months = nstep,
                                years = (time_span[2] - time_span[1]) + 1)
  dimnames(irrigation_scenario) <- list(cells = seq_len(ncell),
                                        months = seq_len(nstep),
                                        years = seq(time_span[1],
                                                    time_span[2]))
  # -------------------------------------------------------------------------- #

  # average discharge reference
  avg_irrigation_scenario <- do.call(average_nyear_window,
                                     append(list(x = irrigation_scenario),
                                            avg_nyear_args))

  third_dim <- names(dim(avg_irrigation_scenario))[
    !names(dim(avg_irrigation_scenario)) %in% c("cells", "months")
  ] %>%
    ifelse(length(.) == 0, NA, .)

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
  dim(drainage) <- c(bands = lpjmlKit::get_header_item(header, "nbands"),
                     cells = lpjmlKit::get_header_item(header, "ncell"))
  dimnames(drainage) <- list(
    bands = seq_len(lpjmlKit::get_header_item(header, "nbands")),
    cells = seq_len(lpjmlKit::get_header_item(header, "ncell")))
  # -------------------------------------------------------------------------- #

  # add 1 since in C indexing starts at 0 but in R at 1
  routing <- drainage[1, ] + 1
  endcell <- array(0, dim = ncell)
  cellindex <-  endcell
  for (cell in 1:ncell) {
    route <- show_route(cell, routing)
    cellindex[route] <- seq(length(route), 1, -1)
    endcell[cell] <- route[length(route)]
  }

  basin_ids <- sort(unique(endcell))
  irrmask_basin <- array(0,
                         dim = dim(avg_irrigation_scenario)[
                           names(dim(avg_irrigation_scenario)) != "months"
                         ],
                         dimnames = dimnames(avg_irrigation_scenario)[
                           names(dimnames(avg_irrigation_scenario)) != "months"
                         ])

  for (id in seq_len(length(basin_ids))) {
    basincell <- which(endcell == basin_ids[id])
    if (is.na(third_dim) | length(dim(drop(avg_irrigation_scenario))) < 3) {
      check_gt0 <- sum(
        asub(avg_irrigation_scenario,
             list("cells" = basincell))
      )
    } else {
      check_gt0 <- apply(
        asub(avg_irrigation_scenario,
             list(cells = basincell)),
        third_dim,
        sum
      )
    }
    basin_replace <- asub(irrmask_basin,
                          list(cells = basincell)) %>%
      `[<-`(check_gt0 > 0, value = 1)

    irrmask_basin <- asub_replace(x = irrmask_basin,
                                  subset_list = list(cells = basincell),
                                  y = basin_replace)
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


# EXPORT TO lpjmlKIT  -------------------------------------------------------- #

# https://stackoverflow.com/questions/47790061/r-replacing-a-sub-array-dynamically
asub_replace <- function(x, subset_list, y) {
  argum <- c(alist(x), subarray_argument(x, subset_list), alist(y))
  do.call("[<-", argum)
}


asub <- function(x, subset_list, drop=TRUE) {
  if (drop) {
    argum <- c(alist(x), subarray_argument(x, subset_list))
  } else {
    argum <- c(alist(x), subarray_argument(x, subset_list), drop = FALSE)
  }
  do.call("[", argum)
}


# https://stackoverflow.com/questions/47790061/r-replacing-a-sub-array-dynamically
subarray_argument <- function(x, subset_list) {
  # first a suitable empty list
  match_x <- which(names(dimnames(x)) %in% names(subset_list))
  match_subset <- na.omit(match(names(dimnames(x)), names(subset_list)))
  subset_list <- mapply(
    function(x, y) {
      if (is.character(x)) {
        return(which(y %in% x))
      } else {
        return(x)
      }
    },
    subset_list[match_subset],
    dimnames(x)[match_x],
    SIMPLIFY = FALSE
  )
  argument <- rep(list(bquote()), length(dim(x)))
  # insert the wanted dimension slices
  argument[match_x] <- subset_list
 return(argument)
}

# ---------------------------------------------------------------------------- #
