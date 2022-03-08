#' Calculate averages (mean) for defined window sizes
#'
#' Define window sizes (nyear) to be used to calculate averages (mean) for each
#' window (`dim(x)[3] / nyear`). Instead of discrete windows, also moving
#' averages can be computed as well as years inbetween interpolated.
#'
#' @param x LPJmL output array with `dim(x)=c(ncell, months, years)`
#'
#' @param nyear integer, if supplied it defines the years for each window to
#' be averaged over in `dim(x)[3]`. If `nyear == 1` values are used directly
#' (instead of calculating an average). nyear has to be smaller than
#' `dim(x)[3]` and `dim(x)[3]` has to be a multiply of nyear.
#' Defaults to `NULL`
#'
#' @param moving_average logical. If `TRUE` moving average is computed. start
#' and end are interpolated using spline interpolation.
#'
#' @param interpolate logical. If `TRUE` and nyear is defined (with
#' `moving_average == FALSE` years are interpolated (spline) to return array
#' with same dimensions as `x` (mainly`dim(x)[3]` -> years).
#'
#' @return array with same amount of cells and months as x. 3rd dimension is
#' defined by nyear, basically `dim(x)[3]/nyear` or equal to dim(x)[3] if
#' `moving_average == TRUE` or `interpolate == TRUE`
#'
#' @md
#' @export
average_nyear_window <- function(x,
                                 nyear = NULL,
                                 moving_average = FALSE,
                                 interpolate = FALSE) {

  # moving average function - spline interpolation to fill NAs at start/end
  moving_average_fun <- function(x, n) {
    stats::filter(x, rep(1 / n, n), sides = 2) %>%
      zoo::na.spline()
  }

  # utility function to interpolate inbetween averaging windows via spline
  interpolate_spline <- function(x, y, nyear) {
    rep(NA, dim(y)[3]) %>%
      `[<-`(seq(round(nyear / 2), dim(y)[3], nyear), value = x) %>%
      zoo::na.spline()
  }

  # if nyear is supplied not all years are used for averaging
  if (!is.null(nyear)) {
    # only valid for nyear <  years of x (dim(x)[3])
    if (nyear > 1 & nyear < dim(x)[3]) {
      # check if multiple (can also be left out)
      if (dim(x)[3] %% nyear == 0) {
        if (moving_average) {
          y <- aperm(apply(x, c(1, 2), moving_average_fun, nyear), c(2, 3, 1))
        } else {
          # calculate mean for discret windows/bins with size of nyear
          y <- array(x, dim = c(dim(x)[1:2], nyear, dim(x)[3] / nyear)
                       ) %>%
            apply(c(1, 2, 4),  mean)
          if (interpolate) {
            y <- aperm(apply(y, c(1, 2), interpolate_spline, x, nyear),
                         c(2, 3, 1))
          }
        }
      }
    } else if (nyear == 1) {
      y <- x
    } else {
      stop(paste0("Amount of nyear (", nyear, ") not supported."))
    }
  } else {
    y <- apply(x, c(1, 2), mean)
  }
  return(y)
}