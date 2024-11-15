#' Calculate environmental flow requirements (EFRs)
#'
#' Calculate environmental flow requirements (EFRs) based on the number of years
#' of `dim(x)[3]`.
#'
#' @param x discharge array with `dim(x)=c(cell, month, year)`. If no year
#' dimension is present, the function assumes a single year.
#'
#' @param approach EFR approach to be used , available methods are `c("vmf",
#' "q90q50")` based on
#' [Pastor et al. 2014](https://doi.org/10.5194/hess-18-5041-2014)
#' and c(vmf_min", "vmf_max") as modified by
#' [Gerten et al. 2020](https://doi.org/10.1038/s41893-019-0465-1)
#' as well as `"steffen2015"`, a modified version of vmf by
#' [Steffen et al. 2015](https://doi.org/10.1126/science.1259855)
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
    c("vmf", "vmf_min", "vmf_max", "q90q50", "steffen2015")
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
    # "vmf_min" - Gerten et al. 2020
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
    # "vmf_max" - Gerten et al. 2020
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
      if (dim(x)["year"] == 1) {
        stop(
          "Approach \"Q90Q50\" is not supported for a single year",
          " as quantiles cannot be calculated"
        )
      } else {
        if (dim(x)["year"] < 30) {
          warning(paste0(
            "Quantile calculation for approach \"Q90Q50\" may not be ",
            "meaningful, as the number of years (", dim(x)["year"],
            ") is less than 30."
          ))
        }
        quantiles <- apply(x,
                           c("cell", "month"),
                           stats::quantile,
                           drop = FALSE,
                           probs = c(0.5, 0.9),
                           na.rm = TRUE)
        q90 <- quantiles[2, , ]
        q50 <- quantiles[1, , ]

        # FS: I dont think this is necessary - removed
        ## repeat quantiles for each year of x so that dimensions match
        #q90 <- rep(q90, dim(x)["year"]) %>%
        #  array(dim = dim(x), dimnames = dimnames(x))
        #q50 <- rep(q50, dim(x)["year"]) %>%
        #  array(dim = dim(x), dimnames = dimnames(x))
        dim(q90) <- dim(efrs)
        dimnames(q90) <- dimnames(efrs)
        dim(q50) <- dim(efrs)
        dimnames(q50) <- dimnames(efrs)

        # low flow months
        efrs[mmf <= maf] <- q90[mmf <= maf]
        # high flow months
        efrs[mmf > maf] <- q50[mmf > maf]

      }
    }
  )
  return(efrs)
}
