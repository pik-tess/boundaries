#' Calculate environmental flow requirements (EFRs)
#'
#' Calculate environmental flow requirements (EFRs) based on the number of years
#' of `dim(x)[3]` or specify a nyear_avg calculate the EFRs for each bin in
#' `dim(x)[3]`.
#'
#' @param x discharge array with `dim(x)=c(cell, month, year)`
#'
#' @param approach EFR approach to be used , available methods are `c("vmf",
#' "vmf_min", "vmf_max", "q90q50")` based on
#' [Pastor et al. 2014](https://doi.org/10.5194/hess-18-5041-2014)
#' as well as `"steffen2015"`, a modified version of vmf by
#' [Steffen et al. 2015](https://doi.org/10.1126/science.1259855)
#'
#' @return EFRs with same unit as `x` (discharge), with `dim(x)=c(ncells, 12)`
#' or `dim(EFRs)=c(ncells, 12, dim(x)[3] / nyear_avg)` if nyear_avg is defined
#'
#' @examples
#' \dontrun{
#' # basic example
#' efrs1 <- calcEFRs(discharge_30y = discharge, approach = "vmf")
#'
#' dim(efrs1)
#' # c(67420, 12)
#'
#' # example for using a 30 year average bin for a 90 year discharge and
#' #  interpolate between 3 windows afterwards to return 90 years (interpolated)
#' efrs2 <- calcEFRs(
#'   discharge_90y = discharge,
#'   approach="vmf"
#' )
#'
#' dim(efrs2)
#' # c(67420, 12, 90)
#' # if interpolate == FALSE  dim(efrs2) returns c(67420, 12, 3)
#'
#' # example for using a 1 year (no average) bin for a 100 year discharge
#' efrs3 <- calcEFRs(discharge_100y = discharge,
#'                   approach = "vmfmin")
#'
#' dim(efrs3)
#' # c(67420, 12, 100)
#' }
#' @md
calc_efrs <- function(x,
                      approach = "vmf") {
  # verify available methods
  approach <- match.arg(approach,
    c("vmf", "vmf_min", "vmf_max", "q90q50", "steffen2015")
  )

  # function to repeat maf mean for dimension length within apply
  maf_fun <- function(x) {
    xm <- mean(x)
    rep(xm, 12)
  }


  # get dimensions without cells to get back standard order cell, month, year
  dim_select <- names(dim(x))[
    which(!names(dim(x)) %in% c("cell", "month"))
  ]
  # calculate maf based on flexibly calculated x (with/out "conserved" years)
  maf <- apply(x, c("cell", dim_select), maf_fun) %>%
    # month dimension is "" here
    aperm(c("cell", "", dim_select))

  # initialize efrs array
  efrs <- array(0, dim = dim(x), dimnames = dimnames(x))

  # apply defined approach
  switch(approach,
    # "vmf" - Pastor et al. 2014
    vmf = {
      # low flow months
      efrs[x <= 0.4 * maf] <- 0.6 * x[x <= 0.4 * maf]
      # intermediate flow months
      efrs[x > 0.4 * maf & x <= 0.8 * maf] <- 0.45 * (
        x[x > 0.4 * maf & x <= 0.8 * maf]
      )
      # high flow months
      efrs[x > 0.8 * maf] <- 0.3 * x[x > 0.8 * maf]
    },
    # "vmf_min" - Pastor et al. 2014
    vmf_min = {
      # low flow months
      efrs[x <= 0.4 * maf] <- 0.45 * x[x <= 0.4 * maf]
      # intermediate flow months
      efrs[x > 0.4 * maf & x <= 0.8 * maf] <- 0.3 * (
        x[x > 0.4 * maf & x <= 0.8 * maf]
      )
      # high flow months
      efrs[x > 0.8 * maf] <- 0.15 * x[x > 0.8 * maf]
    },
    # "vmf_max" - Pastor et al. 2014
    vmf_max = {
      # low flow months
      efrs[x <= 0.4 * maf] <- 0.75 * x[x <= 0.4 * maf]
      # intermediate flow months
      efrs[x > 0.4 * maf & x <= 0.8 * maf] <- 0.6 * (
        x[x > 0.4 * maf & x <= 0.8 * maf]
      )
      # high flow months
      efrs[x > 0.8 * maf] <- 0.45 * x[x > 0.8 * maf]
    },
    # "steffen2015" - Steffen et al. 2015 (adjusted "vmf")
    steffen2015 = {
      # low flow months
      efrs[x <= 0.4 * maf] <- 0.75 * x[x <= 0.4 * maf]
      # intermediate flow months
      efrs[x > 0.4 * maf & x <= 0.8 * maf] <- 0.7 * (
        x[x > 0.4 * maf & x <= 0.8 * maf]
      )
      # high flow months
      efrs[x > 0.8 * maf] <- 0.45 * x[x > 0.8 * maf]
    },
    # "q90q50" - Pastor et al. 2014
    q90q50 = {
      if (dim(x)["year"] == 1) {
        quantiles <- apply(x,
                           c("cell", "month"),
                           stats::quantile,
                           drop = FALSE,
                           probs = c(0.5, 0.9),
                           na.rm = TRUE)
        q90 <- quantiles[2, , ]
        q50 <- quantiles[1, , ]
        # low flow months
        efrs[x <= maf] <- q90[x <= maf]
        # high flow months
        efrs[x > maf] <- q50[x > maf]

      } else {
        stop(
          "Approach \"Q90Q50\" is not supported for time_series_avg",
          " being defined"
        )
      }
    }
  )
  return(efrs)
}
