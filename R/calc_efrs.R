#' Calculate environmental flow requirements (EFRs)
#'
#' Calculate environmental flow requirements (EFRs) based on the number of years
#' of `dim(x)[3]` or specify a bin_size calculate the EFRs for each bin in
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
#' @param bin_size integer, if supplied it defines the years for each bin to be
#' averaged over in `dim(x)[3]`. If `bin_size == 1` values are used directly
#' (instead of calculating an average). bin_size has to be smaller than
#' `dim(x)[3]` and `dim(x)[3]` has to be a multiply of bin_size.
#' Defaults to `NULL`
#'
#' @return EFRs with same unit as `x` (discharge), with `dim(x)=c(ncells, 12)`
#' or `dim(EFRs)=c(ncells, 12, dim(x)[3] / bin_size)` if bin_size is defined
#'
#' @examples
#' \dontrun{
#' # basic example
#' efrs1 <- calcEFRs(discharge_30y = discharge, method = "vmf")
#' }
#' dim(efrs1)
#' # c(67420, 12)
#'
#' # example for using a 30 year average bin for a 90 year discharge
#' efrs2 <- calcEFRs(discharge_90y = discharge, method="vmf", bin_size = 30)
#' }
#' dim(efrs1)
#' # c(67420, 12, 3)
#'
#' # example for using a 1 year (no average) bin for a 100 year discharge
#' efrs3 <- calcEFRs(discharge_100y = discharge, method = "vmf", bin_size = 1)
#' }
#' dim(efrs3)
#' # c(67420, 12, 100)
#'
#' @md
#' @export
calc_efrs <- function(x, method = "vmf", bin_size = NULL) {
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

  # if bin_size (years to average) is supplied
  #   calculate mean monthly flow (mmf) and mean annual flow (maf)
  if (!is.null(bin_size)) {
    if (bin_size > 1 & bin_size < dim(x)[3]) {
      # check if multiple
      if (dim(x)[3] %% bin_size == 0) {
        intervals <- rep(seq(1, dim(x)[3] / bin_size), each = bin_size)
        mmf <- aperm(apply(x, c(1, 2), tapply, intervals, mean), c(2, 3, 1))
      }
    } else if (bin_size == 1) {
      mmf <- x
    } else {
      stop(paste0("Amount of bin_size (", bin_size, ") not supported."))
    }
    maf <- aperm(apply(mmf, c(1, 3), maf_fun), c(2, 1, 3))
  } else {
    mmf <- apply(x, c(1, 2), mean)
    maf <- aperm(apply(mmf, c(1), maf_fun), c(2, 1))
  }

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
      if (!is.null(bin_size)) {
        stop(paste0("Method \"Q90Q50\" is not supported with bin_size ",
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
