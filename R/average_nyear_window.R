#' Calculate averages (mean) for defined window sizes
#'
#' Define window sizes (nyear_window) to be used to calculate averages (mean)
#' for each window (`dim(x)[3] / nyear_window`). Instead of discrete windows,
#' also moving averages can be computed as well as years inbetween interpolated.
#'
#' @param x LPJmL output array with `dim(x)=c(cell, month, year)`
#'
#' @param nyear_window integer, if supplied it defines the years for each window
#' to be averaged over in `dim(x)[3]`. If `nyear_window == 1` values are used
#' directly (instead of calculating an average). nyear_window has to be smaller
#' than `dim(x)[3]` and `dim(x)[3]` is ideally a multipe of nyear_window.
#' Defaults to `NULL`
#'
#' @param moving_average logical. If `TRUE` moving average is computed. start
#' and end are interpolated using spline interpolation.
#'
#' @param interpolate logical. If `TRUE` and nyear_window is defined (with
#' `moving_average == FALSE` years are interpolated (spline) to return array
#' with same dimensions as `x` (mainly`dim(x)[3]` -> year).
#'
#' @param nyear_reference integer, if supplied (default NULL), it defines a
#' time_span for ideally reference runs to be used as a baseline. E.g.
#' `nyear_reference = 30` to be used for preindustrial climate reference.
#'
#' @return array with same amount of cells and months as x. 3rd dimension is
#' defined by nyear_window, basically `dim(x)[3]/nyear_window` or equal to
#' dim(x)[3] if `moving_average == TRUE` or `interpolate == TRUE`
#'
#' @md
#' @importFrom magrittr %>%
#' @export
average_nyear_window <- function(x,
                                 nyear_window = NULL,
                                 moving_average = FALSE,
                                 interpolate = FALSE,
                                 nyear_reference = NULL) {

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

  # moving average function - spline interpolation to fill NAs at start/end
  moving_average_fun <- function(x, n) {
    stats::filter(x, rep(1 / n, n), sides = 2) %>%
      zoo::na.spline()
  }

  # utility function to interpolate inbetween averaging windows via spline
  interpolate_spline <- function(x, y, nyear_window) {
    rep(NA, dim(y)["year"]) %>%
      `[<-`(seq(round(nyear_window / 2), dim(y)["year"], nyear_window),
                value = x) %>%
      zoo::na.spline()
  }

  # if nyear_window is supplied not all years are used for averaging
  if (!is.null(nyear_window)) {
    if (!is.null(nyear_reference)) {
      orig_x <- x
      x <- lpjmlKit::subset_array(x, list(year = 1:nyear_reference))
    }
    # only valid for nyear_window <  years of x (dim(x)[3])
    if (nyear_window > 1 & nyear_window <= dim(x)["year"]) {
      # check if multiple (can also be left out)
      # if (dim(x)[3] %% nyear_window == 0) {
      if (moving_average) {
        y <- aperm(apply(x,
                         c("cell", third_dim),
                         moving_average_fun,
                         nyear_window),
                   c("cell", third_dim, ""))
      } else {
        # calculate mean for discret windows/bins with size of nyear_window
        y <- array(x,
                   # set correct dimensions (with names)
                   dim = c(dim(x)[c("cell", third_dim)],
                           year = nyear_window,
                           window = dim(x)[["year"]] / nyear_window),
                   # set correct dimensions names with nyear and windows
                   dimnames = append(dimnames(x)[c("cell", third_dim)],
                                     list(year = seq_len(nyear_window),
                                          window = dimnames(x)[["year"]][
                                            seq(round(nyear_window / 2),
                                                dim(x)[["year"]],
                                                nyear_window)
                                          ]))) %>%
          apply(c("cell", third_dim, "window"),  mean)
        if (interpolate) {
          y <- aperm(apply(y,
                           c("cell", third_dim),
                           interpolate_spline,
                           x,
                           nyear_window),
                       c("cell", third_dim, ""))
        }
      }
      # }
    } else if (nyear_window == 1) {
      y <- x
    } else {
      stop(paste0("Amount of nyear_window (", nyear_window, ") not supported."))
    }
    # recycle nyear_reference subset for original x (years)
    if (!is.null(nyear_reference)) {
      # multiple factor
      nmultiple <- round(dim(orig_x)[["year"]] / nyear_reference)
      replace_multiple_id <- nmultiple * dim(y)[[3]]
      # if average window is returned as year dimension
      if (!moving_average & !interpolate) {
        z <- array(NA,
                   dim = c(dim(y)[c("cell", third_dim)],
                           window = replace_multiple_id),
                   dimnames = append(dimnames(y)[c("cell", third_dim)],
                                     list(window = rep(dimnames(y)[[3]],
                                                        nmultiple))))
      # return as original year dimension
      } else {
        # years vector also for non multiples (subset only partly recycled)
        years <- rep(NA, dim(orig_x)[["year"]]) %>%
          `[<-`(, value = dimnames(x)[["year"]]) %>%
          suppressWarnings()
        z <- array(NA,
                   dim = dim(orig_x),
                   dimnames = append(dimnames(y)[c("cell", third_dim)],
                                     list(year = years)))
      }
      # recycle subset y for rest of original (x) years in z
      z[, , seq_len(replace_multiple_id)] <- y
      # check if not multiple - then only partly recylce array
      if ((dim(z)[3] - replace_multiple_id) > 0) {
        z[, , (replace_multiple_id + 1):dim(z)[3]] <- (
          y[, , seq_len(dim(z)[3] - replace_multiple_id)]
        )
      }
      return(z)
    } else {
      # rename dimnames of array
      if (all(dim(y) == dim(x))) {
        dim(y) <- dim(x)
        dimnames(y) <- dimnames(x)
      }
    }
  } else {
    y <- apply(x, c("cell", third_dim), mean)
    if (is.null(dim(y))) {
      y <- array(
        y,
        dim = c(cell = dim(x)[["cell"]], 1),
        dimnames = list(cell = dimnames(x)[["cell"]], 1)
      )
    }
  }
  return(y)
}
