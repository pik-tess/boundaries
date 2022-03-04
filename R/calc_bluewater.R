#' Calculate environmental flow requirements (EFRs)
#'
#' Calculate environmental flow requirements (EFRs) based on the number of years
#' of `dim(x)[3]` or an avg_years interval that is recycled for `dim(x)[3]`.
#'
#' @param x discharge array with `dim(x)=c(ncells, months, years)`
#'
#' @param method EFR method to be used , available methods are `c("VMF",
#'
#' @examples
#' \dontrun{
#' # basic example
#' efrs1 <- calcEFRs(discharge_30y = discharge, method = "VMF")
#' }
#' dim(efrs1)
#' # c(67420, 12)

#'
#' @md
#' @export
calc_bluewater <- function(path,
                           method = "gerten2020",
                           time_span = c(1982, 2011),
                           bin_size = NULL,
                           # to be replaced internally by lpjmlKit::read_output
                           start_year = 1901,
                           end_year = 2011) {
  method <- match.arg(method, c("gerten2020",
                                "steffen2015"))

  # TO BE REPLACED BY lpjmlKit::read_output ...
  #   hardcoded values to be internally replaced
  nstep <- 12
  nbands <- 1
  ncell <- 67420
  size <- 4
  file <- file(paste(path, "discharge.bin", sep = "/"), "rb")
  seek(file,
       where = (time_span[1] - start_year) * nstep * nbands * ncell * size,
       origin = "start")
  discharge_baseline <- readBin(file,
                                double(),
                                n = (ncell * nstep * nbands *
                                     (time_span[2] - time_span[1] + 1)),
                                size = size)
  close(file)
  dim(discharge_baseline) <- c(ncell,
                               nstep,
                               (time_span[2] - time_span[1] + 1))

}