#' Calculate environmental flow requirements (EFRs)
#'
#' Calculate environmental flow requirements (EFRs) based on the number of years
#' of `dim(x)[3]`.
#'
#' @param x discharge array with `dim(x)=c(cell, month, year)`. If no year
#' dimension is present, the function assumes a single year.
#'
#' @param approach EFR approach to be used , available methods are `c("vmf",
#' "q10q50")` based on a modified version of
#' [Pastor et al. 2014](https://doi.org/10.5194/hess-18-5041-2014)
#' and c(vmf_min", "vmf_max") as suggested by
#' [Gerten et al. 2020](https://doi.org/10.1038/s41893-019-0465-1)
#' and [Steffen et al. 2015](https://doi.org/10.1126/science.1259855)
#'
#' @return EFRs with same unit as `x` (discharge). Cell and month dimensions are
#' preserved, and the year dimension is set to 1 (as EFRs are calculated
#' based on the whole provided period of years).
#'
#' @examples
#' \dontrun{
#' # basic example
#' efrs <- calcEFRs(discharge_30y = discharge, approach = "vmf")
#'
#' dim(efrs)
#' # c(67420, 12, 1)
#'
#'
#' @md
calc_efrs <- function(x,
                      approach = "vmf") {
  # verify available methods
  approach <- match.arg(approach,
    c("vmf", "vmf_min", "vmf_max", "q10q50")
  )

  # calculate mean monthly flow (mmf)
  mmf <- dimnames_year <- NULL
  if (!is.na(dim(x)["year"])) {
    if (dim(x)["year"] > 1) {
      # average year is not clearly defined for an even amount, we round it down
      dimnames_year <- round(mean(as.numeric(dimnames(x)$year)))
      mmf <- apply(x, c("cell", "month"), mean) %>%
        array(dim = c(dim(x)[1:2], year = 1),
              dimnames = list(cell = dimnames(x)$cell,
                              month = dimnames(x)$month,
                              year = dimnames_year)
        )
    }
  } else {
    dimnames_year <- 1
    mmf <- array(x, dim = c(dim(x)[1:2], year = 1),
                 dimnames = list(cell = dimnames(x)$cell,
                                 month = dimnames(x)$month,
                                 year = dimnames_year))
  }

  # function to repeat maf mean for dimension length within apply
  maf_fun <- function(x) {
    xm <- mean(x)
    rep(xm, 12)
  }

  # calculate mean annual flow (maf)
  maf <- apply(x, c("cell"), maf_fun) %>%
    aperm(c("cell", "")) %>% # "" is the month dimension
    array(dim = c(dim(x)[1:2], year = 1),
          dimnames = list(cell = dimnames(x)$cell,
                          month = dimnames(x)$month,
                          year = dimnames_year))

  # initialize efrs array
  efrs <- array(0, dim = dim(mmf), dimnames = dimnames(mmf))

  # apply defined approach
  switch(approach,
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
    # "vmf_min" - Gerten et al. 2020 (vmf minus 15%)
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
    # "vmf_max" - Gerten et al. 2020 (vmf plus 15%)
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
    # "q10q50" - Pastor et al. 2014
    q10q50 = {
      if (dim(x)["year"] == 1) {
        stop(
          "Approach \"Q10Q50\" is not supported for a single year",
          " as quantiles cannot be calculated"
        )
      } else {
        if (dim(x)["year"] < 30) {
          warning(paste0(
            "Quantile calculation for approach \"Q10Q50\" may not be ",
            "meaningful, as the number of years (", dim(x)["year"],
            ") is less than 30."
          ))
        }
        quantiles <- apply(x,
                           c("cell"),
                           stats::quantile,
                           drop = FALSE,
                           probs = c(0.5, 0.1),
                           na.rm = TRUE)

        q10 <- rep(quantiles["10%", ], times = length(dimnames(efrs)$month))
        q50 <- rep(quantiles["50%", ], times = length(dimnames(efrs)$month))

        # set dim and create dimnames
        q10 <- array(q10, dim = dim(efrs),
                     dimnames = dimnames(efrs))
        q50 <- array(q50, dim = dim(efrs),
                     dimnames = dimnames(efrs))

        # low flow months
        efrs[mmf <= maf] <- q10[mmf <= maf]
        # high flow months
        efrs[mmf > maf] <- q50[mmf > maf]

      }
    }
  )
  return(efrs)
}
