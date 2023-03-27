
path_scenario <- (
  "/p/projects/open/Jannes/earth4all/pb_status/runs/output/lu_old/"
)

test <- LPJmLTools::lpjmlinfo$cells_raster
test[LPJmLTools::lpjmlinfo$cellnumbers] <- pb_status# [, 1]
plot(test)




time_span_scenario <- c(1982, 2011)
start_year <- 1901
end_year <- 2011

nstep <- 12
nbands <- 1
ncell <- 67420
size <- 4
# TO BE REPLACED BY lpjmlKit::read_output ---------------------------------- #
#   hardcoded values to be internally replaced
s_path <- file(paste(path_scenario, "irrig.bin", sep = ""), "rb")
seek(s_path,
     where = (time_span_scenario[1] - start_year) *
             nstep * nbands * ncell * size,
     origin = "start")
irrigation_scenario <- readBin(s_path,
                             double(),
                              n = (ncell * nstep * nbands *
                                   (time_span_scenario[2] -
                                    time_span_scenario[1] + 1)),
                             size = size)
close(s_path)
dim(irrigation_scenario) <- c(cells = ncell,
                             months = nstep,
                             years = (time_span_scenario[2] -
                                      time_span_scenario[1] + 1))
dimnames(irrigation_scenario) <- list(cells = seq_len(ncell),
                                     months = seq_len(nstep),
                                     years = seq(time_span_scenario[1],
                                                 time_span_scenario[2]))
# -------------------------------------------------------------------------- #

# average discharge reference
avg_irrigation_scenario <- do.call(average_nyear_window,
                                   append(list(x = irrigation_scenario),
                                          avg_nyear_args))

third_dim <- names(dim(avg_irrigation_scenario))[
  !names(dim(avg_irrigation_scenario)) %in% c("cells", "months")
] %>%
  ifelse(length(.) == 0, NA, .)

drainage <- autoReadInput(
    inFile = "/p/projects/lpjml/input/historical/input_VERSION2/drainage.bin",
    manu = TRUE,
    msize = 4,
    mheadersize = 43
)
# add 1 since in C indexing starts at 0 but in R at 1
routing <- drainage[1, ] + 1
# indices of ocean draining cells (due to adding 1, -1 becomes 0)
endcell_indices <- which(routing == 0)
# cells that were not part of STN, LPJmL treats all values below 0 as outflow
#   cells
na_cells <- which(routing == -8)
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
  if (is.na(third_dim)) {
    check_gt0 <- sum(
      lpjmlkit::asub(avg_irrigation_scenario, cells = basincell)
    )
  } else {
    check_gt0 <- apply(
      lpjmlkit::asub(avg_irrigation_scenario, cells = basincell),
      "",
      sum
    )
  }
  basin_replace <- lpjmlkit::asub(irrmask_basin, cells = basincell) %>%
    `[<-`(check_gt0 > 0, value = 1)

  irrmask_basin <- lpjmlkit::asub(irrmask_basin,
                      which(names(dimnames(irrmask_basin)) == "cells"),
                      basincell,
                      basin_replace)
}


test <- LPJmLTools::lpjmlinfo$cells_raster
test[LPJmLTools::lpjmlinfo$cellnumbers] <- irrmask_basin[, 1]
plot(test)


show_route <- function(ind, routing_table) {
  if (routing_table[ind] < 1) {
    # can be 0 or -8 -> endcell or nacell
    return(ind)
  } else {
    return(c(ind,
             show_route(routing_table[ind], routing_table)))
  }
}


# https://stackoverflow.com/questions/47790061/r-replacing-a-sub-array-dynamically
asub_replace <- function(x, dims, idx, y){
  argum <- c(alist(x), subarray_argument(x, dims, idx), alist(y))
  do.call("[<-", argum)
}


# https://stackoverflow.com/questions/47790061/r-replacing-a-sub-array-dynamically
subarray_argument <- function(x, dims, idx){
  if (class(idx) != "list") idx <- list(idx)
  dim_x <- dim(x)
  dim_length <- length(dim_x)
  stopifnot(all(dims >= 0) & all(dims <= dim_length))
  stopifnot(dim_x[dims] >= lapply(idx, max))
  # first a suitable empty list
  argument <- rep(list(bquote()),dim_length)
  # insert the wanted dimension slices
  argument[dims] <- idx
 return(argument)
}