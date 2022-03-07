#' Calculate bluewater planetary boundary status
#'
#' Calculate ...
#'
#' @param path
#'
#' @param method
#'
#' @param time_span
#'
#' @param bin_size
#'
#' @examples
#' \dontrun{
#'  calc_bluewater(path)
#' }
#'
#' @md
#' @export
calc_bluewater <- function(path,
                           method = "gerten2020",
                           time_span = c(1982, 2011),
                           n_avg = NULL,
                           # to be replaced internally by lpjmlKit::read_output
                           start_year = 1901,
                           end_year = 2011) {
  # verify available methods
  method <- match.arg(method, c("gerten2020",
                                "steffen2015"))

  # read baseline discharge
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

  efr_req_uncertainty <- calc_efrs(discharge_baseline, "VMFmin")
  efr_req_safe <- calc_efrs(discharge_baseline, "vmf_max")

  avg_disch_baseline <- apply(discharge_baseline, c(1, 2), mean) # mean of years

}