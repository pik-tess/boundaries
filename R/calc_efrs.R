#' Calculate environmental flow requirements (EFRs)
#'
#' Calculate environmental flow requirements (EFRs) based on the number of years
#' of `dim(x)[3]` or an avg_years interval that is recycled for `dim(x)[3]`.
#'
#' @param x discharge array with `dim(x)=c(ncells, months, years)`
#'
#' @param method EFR method to be used , available methods are `c("VMF",
#' "VMFmin", "VMFmax", "Q90Q50")` based on
#' [Pastor et al. 2014](https://doi.org/10.5194/hess-18-5041-2014)
#' as well as `"PBpaper"`, a modified version of VMF by
#' [Steffen et al. 2015](https://doi.org/10.1126/science.1259855)
#'
#' @param avg_years integer, if supplied it defines the years to average over,
#' avg_years is recycled for all years `dim(x)[3]`. If `avg_years == 1` values
#' are used directly (instead of calculating an average).
#' avg_years has to smaller than `dim(x)[3]` and `dim(x)[3]` has to be a 
#' multiply of avg_years. Defaults to `NULL`
#'
#' @return EFRs with same unit as `x` (discharge), with `dim(x)=c(ncells, 12)`
#' or `dim(EFRs)=c(ncells, 12, dim(x)[3] / avg_years)` if avg_years is defined
#'
#' @examples
#' \dontrun{
#' # basic example
#' efrs1 <- calcEFRs(discharge_30y = discharge, method = "VMF")
#' }
#' dim(efrs1)
#' # c(67420, 12)
#'
#' # example for using 30 year avg_years interval for a 90 year discharge
#' efrs2 <- calcEFRs(discharge_90y = discharge, method="VMF", avg_years = 30)
#' }
#' dim(efrs1)
#' # c(67420, 12, 3)
#'
#' # example for using a 1/no year avg_years interval for a 100 year discharge
#' efrs3 <- calcEFRs(discharge_100y = discharge, method = "VMF", avg_years = 1)
#' }
#' dim(efrs3)
#' # c(67420, 12, 100)
#'
#' @md
#' @export
calc_efrs <- function(x, method = "VMF", avg_years = NULL) {
  method <- match.arg(method, c("VMFmin",
                                "VMF",
                                "VMFmax",
                                "PBpaper",
                                "Q90Q50"))

  # function to repeat maf mean for dimension length within apply
  maf_fun <- function(x) {
    xm <- mean(x)
    rep(xm, 12)
  }

  # if avg_years (years to average) is supplied
  #   calculate mean monthly flow (mmf) and mean annual flow (maf)
  if (avg_years) {
    if (avg_years > 1 & avg_years < dim(x)[3]) {
      # check if multiple
      if (dim(x)[3] %% avg_years == 0) {
        intervals <- rep(seq(1, dim(x)[3] / avg_years), each = avg_years)
        mmf <- aperm(apply(x, c(1, 2), tapply, intervals, mean), c(2, 3, 1))
      }
    } else if (avg_years == 1) {
      mmf <- x
    } else {
      stop(paste0("Amount of avg_years (", avg_years, ") not supported."))
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
    # "VMF" - Pastor et al. 2014
    VMF = {
      # low flow months
      efrs[mmf <= 0.4 * maf] <- 0.6 * mmf[mmf <= 0.4 * maf]
      # intermediate flow months
      efrs[mmf > 0.4 * maf & mmf <= 0.8 * maf] <- 0.45 * (
        mmf[mmf > 0.4 * maf & mmf <= 0.8 * maf]
      )
      # high flow months
      efrs[mmf > 0.8 * maf] <- 0.3 * mmf[mmf > 0.8 * maf]
    },
    # "VMFmin" - Pastor et al. 2014
    VMFmin = {
       # low flow months
      efrs[mmf <= 0.4 * maf] <- 0.45 * mmf[mmf <= 0.4 * maf]
       # intermediate flow months
      efrs[mmf > 0.4 * maf & mmf <= 0.8 * maf] <- 0.3 * (
        mmf[mmf > 0.4 * maf & mmf <= 0.8 * maf]
      )
      # high flow months
      efrs[mmf > 0.8 * maf] <- 0.15 * mmf[mmf > 0.8 * maf]
    },
    # "VMFmax" - Pastor et al. 2014
    VMFmax = {
      # low flow months
      efrs[mmf <= 0.4 * maf] <- 0.75 * mmf[mmf <= 0.4 * maf]
      # intermediate flow months
      efrs[mmf > 0.4 * maf & mmf <= 0.8 * maf] <- 0.6 * (
        mmf[mmf > 0.4 * maf & mmf <= 0.8 * maf]
      )
      # high flow months
      efrs[mmf > 0.8 * maf] <- 0.45 * mmf[mmf > 0.8 * maf]
    },
    # "PBpaper" - Steffen et al. 2015 (adjusted "VMF")
    PBpaper = {
      # low flow months
      efrs[mmf <= 0.4 * maf] <- 0.75 * mmf[mmf <= 0.4 * maf]
      # intermediate flow months
      efrs[mmf > 0.4 * maf & mmf <= 0.8 * maf] <- 0.7 * (
        mmf[mmf > 0.4 * maf & mmf <= 0.8 * maf]
      )
      # high flow months
      efrs[mmf > 0.8 * maf] <- 0.45 * mmf[mmf > 0.8 * maf]
    },
    # "Q90Q50" - Pastor et al. 2014
    Q90Q50 = {
      if (!is.null(avg_years)) {
        stop(paste0("Method \"Q90Q50\" is not supported with avg_years ",
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
