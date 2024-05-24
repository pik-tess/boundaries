#' Calculate averages (mean) for defined window sizes
#'
#' Define window sizes (time_resolution) to be used to calculate moving averages
#' (mean). If time_resolution is not supplied, the function calculates the mean
#' over all years. If time_repeat is supplied, the function replicates the
#' mean values for the defined amount of years.
#'
#' @param x LPJmL output array with `dim(x)=c(cell, month, year)`
#'
#' @param time_resolution integer. Number of years to be used for the moving
#' average calculation. If `NULL`, all years are averaged for
#' `spatial_scale = c("grid", "subglobal")` or the whole time span is used for
#' `spatial_scale = "global"`.
#'
#' @param time_repeat integer, if supplied (default NULL), it defines a
#' length of years to be replicated. Only if time_resolution is not supplied.
#'
#' @return array with same amount of cells and months as x if time_resolution is
#' supplied. If time_repeat is supplied, the array has the same amount of
#' cells and months as x but the amount of years is multiplied by
#' time_repeat.
#'
#' @md
aggregate_time <- function(x,
                           time_resolution = NULL,
                           time_repeat = NULL) {

  if (!is.array(x)) {
    x <- array(
      x,
      dim = c(cell = 1, year = length(x)),
      dimnames = list(cell = "global", year = names(x))
    )
  }

  if (is.null(time_resolution)) {
    time_resolution <- dim(x)["year"]
  }

  third_dim <- names(dim(x))[
    !names(dim(x)) %in% c("cell", "year")
  ] %>% {
    if (rlang::is_empty(.)) NULL else .
  }

  if (length(third_dim) > 1) {
    stop(paste0("x has to have dimensions \"cell\", \"year\" and can have ",
                "one third dimension (e.g. \"month\", \"year\""))
  }
  # check validity of x dimensions
  if (!all(names(dim(x)) %in% c("cell", third_dim, "year"))) {
    stop(paste0("x has to have dimensions \"cell\"",
                ifelse(is.null(third_dim),
                       "",
                       paste0(", \"", third_dim, "\"")),
                " and \"year\"."))
  }

  # if time_resolution is supplied not all years are used for averaging
  if (time_resolution == 1) {
    y <- x

  } else if (time_resolution > 1 & time_resolution < dim(x)["year"]) { # nolint:vector_logic_linter
    y <- aperm(apply(x,
                     c("cell", third_dim),
                     moving_average_fun,
                     time_resolution),
               c("cell", third_dim, ""))
    dim(y) <- dim(x)
    dimnames(y) <- dimnames(x)

  } else if (time_resolution == dim(x)["year"]) {

    y <- apply(x, c("cell", third_dim), mean)
    y_dim <- dim(x)
    y_dim[["year"]] <- 1

    y_dimnames <- dimnames(x)
    y_dimnames[["year"]] <- round(mean(as.integer(dimnames(x)[["year"]])))


    y <- array(
      y,
      dim = y_dim,
      dimnames = y_dimnames
    )

    if (!is.null(time_repeat)) {
      y_dimnames[["year"]] <- rep(y_dimnames[["year"]], time_repeat)
      y <- abind::abind(
        mget(rep("y", time_repeat)),
        use.dnns = TRUE
      )
      dim(y) <- sapply(y_dimnames, length) # nolint:undersirable_function_linter
      dimnames(y) <- y_dimnames
    }

  } else {
    stop(paste0("Amount of time_resolution (", time_resolution, ") not supported."))
  }

  if (all(dimnames(x)[["cell"]] == "global")) {
    y <- abind::adrop(y, "cell")
  }

  return(y)
}


# moving average function with shrinking window size at edges
moving_average_fun <- function(x, n) {
  len <- length(x)
  half_n <- n %/% 2

  # Compute the left edge with shrinking window size
  left_edge <- sapply( # nolint:undesirable_function_linter
    1:half_n,
    function(i) mean(x[1 : (2 * i - 1)])
  )

  # Compute the main section with fixed window size
  main_section <- sapply( # nolint:undesirable_function_linter
    (half_n + 1) : (len - half_n),
    function(i) mean(x[(i - half_n) : (i + half_n)])
  )

  # Compute the right edge with shrinking window size
  right_edge <- sapply( # nolint:undesirable_function_linter
    1 : half_n,
    function(i) mean(x[(len - (2 * i - 2)) : len])
  )
  right_edge <- rev(right_edge)  # Reverse the right edge to correctly align

  # Combine all parts into the final result
  result <- c(left_edge, main_section, right_edge)

  return(result)
}
