#' Calculate environmental flow requirements (EFRs)
#'
#' Calculate environmental flow requirements (EFRs) based on the number of years
#' of `dim(x)[3]` or specify a nyear_avg calculate the EFRs for each bin in
#' `dim(x)[3]`.
#'
#' @param x discharge array with `dim(x)=c(ncells, months, years)`
#'
#' @param method EFR method to be used , available methods are `c("vmf",
#' "vmf_min", "vmf_max", "q90q50")` based on
#' [Pastor et al. 2014](https://doi.org/10.5194/hess-18-5041-2014)
#' as well as `"steffen2015"`, a modified version of vmf by
#' [Steffen et al. 2015](https://doi.org/10.1126/science.1259855)
#'
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#' time series analysis.
#'
#' @return EFRs with same unit as `x` (discharge), with `dim(x)=c(ncells, 12)`
#' or `dim(EFRs)=c(ncells, 12, dim(x)[3] / nyear_avg)` if nyear_avg is defined
#'
#' @examples
#' \dontrun{
#' # basic example
#' efrs1 <- calcEFRs(discharge_30y = discharge, method = "vmf")
#'
#' dim(efrs1)
#' # c(67420, 12)
#'
#' # example for using a 30 year average bin for a 90 year discharge and
#' #  interpolate between 3 windows afterwards to return 90 years (interpolated)
#' efrs2 <- calcEFRs(discharge_90y = discharge,
#'                   method="vmf",
#'                   avg_nyear_args = list(nyear_avg = 30, interpolate = TRUE))
#'
#' dim(efrs2)
#' # c(67420, 12, 90)
#' # if interpolate == FALSE  dim(efrs2) returns c(67420, 12, 3)
#'
#' # example for using a 1 year (no average) bin for a 100 year discharge
#' efrs3 <- calcEFRs(discharge_100y = discharge,
#'                   method = "vmfmin",
#'                   avg_nyear_args = list(nyear_avg = 1))
#'
#' dim(efrs3)
#' # c(67420, 12, 100)
#' }
#' @md
#' @export
calc_efrs <- function(x,
                      method = "vmf",
                      avg_nyear_args = list()) {
  # verify available methods
  method <- match.arg(method, c("vmf",
                                "vmf_min",
                                "vmf_max",
                                "q90q50",
                                "steffen2015"))

  # function to repeat maf mean for dimension length within apply
  maf_fun <- function(x) {
    xm <- mean(x)
    rep(xm, 12)
  }

  # if nyear_avg (years to average) is supplied
  #   calculate mean monthly flow (mmf) and mean annual flow (maf)
  mmf <- do.call(average_nyear_window, append(list(x = x), avg_nyear_args))
  # get dimensions without cells to get back standard order cells, months, years
  dim_select <- names(dim(mmf))[which(names(dim(mmf)) != "cells")]
  # calculate maf based on flexibly calculated mmf (with/out "conserved" years)
  maf <- aperm(apply(mmf, dim_select, maf_fun), c("cells", dim_select))

  # initialize efrs array
  efrs <- array(0, dim = dim(mmf))

  # apply defined method
  switch(method,
    # "vmf" - Pastor et al. 2014
    vmf = {
      # low flow months
      efrs[mmf <= 0.4 * maf] <- 0.6 * mmf[mmf <= 0.4 * maf]
      # intermediate flow months
      efrs[mmf > 0.4 * maf & mmf <= 0.8 * maf] <- 0.45 * (
        mmf[mmf > 0.4 * maf & mmf <= 0.8 * maf]
      )
      # high flow months
      efrs[mmf > 0.8 * maf] <- 0.3 * mmf[mmf > 0.8 * maf]
    },
    # "vmf_min" - Pastor et al. 2014
    vmf_min = {
       # low flow months
      efrs[mmf <= 0.4 * maf] <- 0.45 * mmf[mmf <= 0.4 * maf]
       # intermediate flow months
      efrs[mmf > 0.4 * maf & mmf <= 0.8 * maf] <- 0.3 * (
        mmf[mmf > 0.4 * maf & mmf <= 0.8 * maf]
      )
      # high flow months
      efrs[mmf > 0.8 * maf] <- 0.15 * mmf[mmf > 0.8 * maf]
    },
    # "vmf_max" - Pastor et al. 2014
    vmf_max = {
      # low flow months
      efrs[mmf <= 0.4 * maf] <- 0.75 * mmf[mmf <= 0.4 * maf]
      # intermediate flow months
      efrs[mmf > 0.4 * maf & mmf <= 0.8 * maf] <- 0.6 * (
        mmf[mmf > 0.4 * maf & mmf <= 0.8 * maf]
      )
      # high flow months
      efrs[mmf > 0.8 * maf] <- 0.45 * mmf[mmf > 0.8 * maf]
    },
    # "steffen2015" - Steffen et al. 2015 (adjusted "vmf")
    steffen2015 = {
      # low flow months
      efrs[mmf <= 0.4 * maf] <- 0.75 * mmf[mmf <= 0.4 * maf]
      # intermediate flow months
      efrs[mmf > 0.4 * maf & mmf <= 0.8 * maf] <- 0.7 * (
        mmf[mmf > 0.4 * maf & mmf <= 0.8 * maf]
      )
      # high flow months
      efrs[mmf > 0.8 * maf] <- 0.45 * mmf[mmf > 0.8 * maf]
    },
    # "q90q50" - Pastor et al. 2014
    q90q50 = {
      if (length(dim(mmf)) == 3) {
        stop(paste0("Method \"Q90Q50\" is not supported for avg_nyear_args ",
                    "being defined"))
      }
      quantiles <- apply(x,
                         c(1, 2),
                         quantile,
                         probs = c(0.5, 0.9),
                         na.rm = TRUE)
      q90 <- quantiles[2, , ]
      q50 <- quantiles[1, , ]
      # low flow months
      efrs[mmf <= maf] <- q90[mmf <= maf]
      # high flow months
      efrs[mmf > maf] <- q50[mmf > maf]
    }
  )
  return(efrs)
}
