rm(list = ls(all = TRUE))
gc()

library(lpjmliotools)
library(raster)
library(fields)
library(rgeos)


# paths ---------------------------------------------------------------------- #
run <- "ostberg" # old
dir_inputs <- "/p/projects/open/Vera/synthesepaper/input_data/"
dir_r_out <- "/p/projects/open/Jannes/earth4all/pb_status/R/data/"
dir_lpj_out <- "/p/projects/open/Jannes/earth4all/pb_status/runs/output/"


# arguments & params --------------------------------------------------------- #
ncell <- 67420
timespan <- c(1988, 2017)
ref <- 1901
firstyear <- timespan[1]
nyears <- timespan[2] - timespan[1] + 1
temp_resolution <- 12
nbands <- 1
size <- 4
# 1 hm3/d = 100^3 m3/ 86400 s = 11.574 m3/s
# in order dismiss flows < 0.1 m3/s, translate to 0.00864
cut <- 0.0864 # cut at 1 m3/s

zz <- file(paste(dir_lpj_out, "/", run, "_pnv/grid.bin", sep = ""), "rb")
x <- readBin(zz, integer(), n = 2 * ncell, size = 2) / 100
lon <- x[c(1:ncell) * 2 - 1]
lat <- x[c(1:ncell) * 2]
cellarea <- (111e3 * 0.5) * (111e3 * 0.5) * cos(lat / 180 * pi)
close(zz)


# discharge pnv -------------------------------------------------------------- #

file <- file(paste0(dir_lpj_out, "/", run, "_pnv", "/", "mdischarge.bin"), "rb")
seek(file,
     where = (timespan[1] - ref) * temp_resolution * nbands * ncell * size,
     origin = "start")
disch_baseline_allyears <- readBin(file,
                         double(),
                         n = ncell * temp_resolution * nyears * nbands,
                         size = size)
close(file)
# dim_LPJoutput<- c(ncell,12,nyears)
dim(disch_baseline_allyears) <- c(ncell, temp_resolution, nyears)



calc_efrs <- function(discharge, method = "VMF") {
  ncells <- dim(discharge)[1]
  set_method <- c("VMFmin", "VMF", "VMFmax")
  #method valid?
  if (!method %in% set_method) {
    print(paste(
      "EFR method not recognized, use one of: ",
      paste(set_method, collapse = " ")
    ))
    return(NA)
  }
  #make sure discharge is a ncells,12month,30year array
  if (!all(dim(discharge) == c(ncells, 12, 30))) {
    print(paste(
      "discharge array has wrong dimension,",
      "use c(ncells,nmonths=12,nyears=30)"
    ))
    return(NA)
  }
  efrs <- discharge[, , 1] * 0 #initialize efrs array

  #calculate mean monthly flow and mean annual flow
  mmf <- apply(discharge, c(1, 2), mean)
  maf <- rep(rowMeans(mmf), times = 12)
  dim(maf) <- c(ncells, 12)

  if (method == set_method[1]) { # "VMFmin" - Pastor et al. 2014
    efrs[mmf <= 0.4 * maf] <- 0.45 * mmf[mmf <= 0.4 * maf] # low flow months
    efrs[mmf > 0.4 * maf & mmf <= 0.8 * maf] <- 0.3 * (
      mmf[mmf > 0.4 * maf & mmf <= 0.8 * maf]
    ) # intermediate flow months
    efrs[mmf > 0.8 * maf] <- 0.15 * mmf[mmf > 0.8 * maf] # high flow months
  }else if (method == set_method[2]) { # "VMF" - Pastor et al. 2014
    efrs[mmf <= 0.4 * maf] <- 0.6 * mmf[mmf <= 0.4 * maf] # low flow months
    efrs[mmf > 0.4 * maf & mmf <= 0.8 * maf] <- 0.45 * (
      mmf[mmf > 0.4 * maf & mmf <= 0.8 * maf]
    ) # intermediate flow months
    efrs[mmf > 0.8 * maf] <- 0.3 * mmf[mmf > 0.8 * maf] # high flow months
  }else if (method == set_method[3]) { # "Q90Q50" - Pastor et al. 2014
    # "VMFmin" - Pastor et al. 2014
    efrs[mmf <= 0.4 * maf] <- 0.75 * mmf[mmf <= 0.4 * maf]
    efrs[mmf > 0.4 * maf & mmf <= 0.8 * maf] <- 0.6 * (
      mmf[mmf > 0.4 * maf & mmf <= 0.8 * maf]
    ) # intermediate flow months
    efrs[mmf > 0.8 * maf] <- 0.45 * mmf[mmf > 0.8 * maf] # high flow months
  }
  return(efrs)
}

efr_req_uncertainty <- calc_efrs(disch_baseline_allyears, "VMFmin")
efr_req_safe <- calc_efrs(disch_baseline_allyears, "VMFmax")

disch_baseline <- apply(disch_baseline_allyears, c(1, 2), mean) # mean of years

# discharge landuse ---------------------------------------------------------- #
file <- file(paste(dir_lpj_out, run, "_lu/", "mdischarge.bin", sep = ""), "rb")
seek(file,
     where = (timespan[1] - ref) * temp_resolution * nbands * ncell * size,
     origin = "start")
disch_lu_allyears <- readBin(file,
                             double(),
                             n = ncell * temp_resolution * nyears * nbands,
                             size = size)
close(file)
dim(disch_lu_allyears) <- c(ncell, temp_resolution, nyears) #

disch_lu <- apply(disch_lu_allyears, c(1, 2), mean) # mean over years


# post calculation of EFRs --------------------------------------------------- #

# for (month in 1:12) {
#   # vmfmax (vmf + 15%)
#   efr_req_safe[, month] <- ifelse(
#     efr_req_safe[, month] > disch_baseline[, month],
#     disch_baseline[, month],
#     efr_req_safe[, month]
#   )
#   # vmfmin (vmf - 15%)
#   efr_req_uncertainty[, month] <- ifelse(
#     efr_req_uncertainty[, month] > disch_baseline[, month],
#     disch_baseline[, month],
#     efr_req_uncertainty[, month]
#   )
# }


# uncertainty zone
span <- efr_req_safe - efr_req_uncertainty
# set considered discharge to zero if PNV < 0.1 m^3/s
span[disch_baseline < cut] <- 0

# def_pnv_safe <- efr_req_safe - disch_baseline
# # if disch(PNV) > efr_reqired_safe
# def_pnv_safe[def_pnv_safe < 0] <- 0
# # all cells where PNV discharge is smaller than efr requirements have a value
# def_pnv_safe[disch_baseline < cut] <- 0

# load basin mask
load(paste(dir_r_out, "disch_irrmask_basin_baseline_noefr.RDATA"))

# disch_target - actual-disch <-- efr transgressed if larger 0
def_safe <- efr_req_safe - disch_lu
# substract amout of EFR transgression in PNV --> reduce deficit of scenario by
#   PNV deficit
# def_safe <- def_safe - def_pnv_safe

def_safe_prelim <- def_safe
# remove cells where disch+0.1 > EFR-req
def_safe[def_safe < cut] <- 0
def_safe[disch_pnv < cut] <- 0

# as in Steffen 2015: degree to which EFRs are undermined: expressed as the
#   transgression-to uncertainty ratio (>5% within uncertainty range (yellow),
#   >75% transgressed)
frac <- ifelse(span > 0, def_safe / span, 0)
# EFR deficit smaller 0.1 hm^3/d is not considered
frac[def_safe < cut] <- NA
# ask Jonas: why 0.05?  iguess because <5% transg = safe--> safe is taken out
#   for averaging
frac[frac <= 0.05] <- NA

pbw_frac_mean <- apply(frac, 1, mean, na.rm = T)

cells_maginal_disch <- array(0, ncell)
cells_maginal_disch[which(apply(disch_lu, 1, mean) < cut)] <- (-1)


pb_range <- c(0, 100)
pb_freshwater <- pbw_frac_mean * 100
pb_freshwater[irrmask_basin == 0] <- NA
pb_freshwater[pb_freshwater > pb_range[2]] <- pb_range[2]
pb_freshwater[pb_freshwater < pb_range[1]] <- pb_range[1]
pb_freshwater[cells_maginal_disch < 0] <- (-1)


saveRDS(pb_freshwater,
        file = paste(dir_r_out, "pb_freshwater.rds", sep = ""))
