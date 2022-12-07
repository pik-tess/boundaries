################################################################################
###          Calculation of the status of the planetare boundary            ####
###                           for freshwater use                            ####
################################################################################

# based on the calculation of EFRS with the Variable monthly flow method
# (Pastor et al. 2014), as proposed in Steffen et al. 2015
# boundary status derived from the EFR-transgression-to-uncertainty-ratio

rm(list = ls(all = TRUE))
gc()

library(lpjmliotools)
library(raster)
library(fields)
library(rgeos)

# paths ---------------------------------------------------------------------- #
run <- "old" # ostberg
dir_inputs <- "/p/projects/open/Vera/synthesepaper/input_data/"
dir_r_out <- "/p/projects/open/Jannes/earth4all/pb_status/R/data/"
dir_lpj_out <- "/p/projects/open/Jannes/earth4all/pb_status/runs/output/"
dir_lpj_out_pnv <- paste0(dir_lpj_out, "/", run, "_pnv", "/")
# dir_lpj_out_pnv <- "/p/projects/open/Jannes/lpjml/testing/nitrogen_fixing/output/noirrig_lu/"
dir_lpj_out_lu <- paste0(dir_lpj_out, "/", run, "_lu", "/")


# arguments & params --------------------------------------------------------- #
ncell <- 67420
# timespan for PB status calculation
timespan_results <- c(1988, 2017)
# timespan for EFR calculation (discharge reference)
timespan_ref <- c(1988, 2017)
ref_lu <- 1901 #1950 #1901 # first year of lu simulation
ref <- 1901 #first simulation year of reference run for EFR calculation
nyears_ref <- timespan_ref[2] - timespan_ref[1] + 1
nyears_lu <- timespan_results[2] - timespan_results[1] + 1

temp_resolution <- 12
nbands <- 1
size <- 4
# 1 hm3/d = 100^3 m3/ 86400 s = 11.574 m3/s
# in order dismiss flows < 1 m3/s, translate to 0.0864
cut <- 0.0864 # cut at 1 m3/s

# read in reference discharge for EFR calculation ---------------------------- #

# usually this is a PNV output
file <- file(paste0(dir_lpj_out_pnv, "mdischarge.bin"), "rb")
seek(file,
     where = (timespan_ref[1] - ref) * temp_resolution * nbands * ncell * size,
     origin = "start")
disch_baseline_allyears <- readBin(file,
                         double(),
                         n = ncell * temp_resolution * nyears_ref * nbands,
                         size = size)
close(file)
dim(disch_baseline_allyears) <- c(ncell, temp_resolution, nyears_ref)

disch_baseline <- apply(disch_baseline_allyears, c(1, 2), mean) # mean of years

# calculate EFRs based on reference discharge -------------------------------- #

# function to calculate EFRs based on reference discharge
# adapted from lpjmliotools R package (developed by Fabian)
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
  #if (!all(dim(discharge) == c(ncells, 12, 30))) {
  #  print(paste(
  #    "discharge array has wrong dimension,",
  #    "use c(ncells,nmonths=12,nyears=30)"
  #  ))
  #  return(NA)
  #}
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


# read in discharge landuse -------------------------------------------------- #
file <- file(paste(dir_lpj_out_lu, "mdischarge.bin", sep = ""), "rb")
seek(file,
     where = (timespan_results[1] - ref_lu) * temp_resolution * nbands * ncell *
     size, origin = "start")
disch_lu_allyears <- readBin(file,
                             double(),
                             n = ncell * temp_resolution * nyears_lu * nbands,
                             size = size)
close(file)
dim(disch_lu_allyears) <- c(ncell, temp_resolution, nyears_lu) #

disch_lu <- apply(disch_lu_allyears, c(1, 2), mean) # mean over years

# calculation of EFR transgressions = EFR deficits in LU run ----------------- #

def_safe <- efr_req_safe - disch_lu

# calculation of uncertainty zone -------------------------------------------- #
span <- efr_req_safe - efr_req_uncertainty

# dismiss boundary status calculation if PNV discharge is < 0.1 m^3/s -------- #
span[disch_baseline < cut] <- 0
def_safe[def_safe < cut] <- 0
def_safe[disch_baseline < cut] <- 0

# load basin mask to mask out river basins without irrigation ---------------- #
load(paste(dir_r_out, "disch_irrmask_basin_baseline_noefr.RDATA"))

# calculate boundary status based on transgression to uncertainty ratio ------ #

# as in Steffen 2015: degree to which EFRs are undermined: expressed as the
# transgression-to uncertainty ratio
# if ratio is above >5%: within uncertainty range (yellow)
# if ratio is above >75% transgression (red)

frac <- ifelse(span > 0, def_safe / span, 0)

# set ratio to NA if discharge is smaller than 1 m3/s
frac[def_safe < cut] <- NA

# to average the ratio only over months which are not "safe",
# frac <= 0.05 is set to NA
frac[frac <= 0.05] <- NA

# boundary status is averaged over months with a transgression (yellow or red)
pbw_frac_mean <- apply(frac, 1, mean, na.rm = T)

# to display cells with marginal discharge in other color (grey):
cells_maginal_disch <- array(0, ncell)
cells_maginal_disch[which(apply(disch_lu, 1, mean) < cut)] <- (-1)

# prepare data for plotting
pb_range <- c(0, 100)
pb_freshwater <- pbw_frac_mean * 100
pb_freshwater[irrmask_basin == 0] <- NA
pb_freshwater[pb_freshwater > pb_range[2]] <- pb_range[2]
pb_freshwater[pb_freshwater < pb_range[1]] <- pb_range[1]
pb_freshwater[cells_maginal_disch < 0] <- (-1)


saveRDS(pb_freshwater,
        file = paste(dir_r_out, "pb_freshwater_", run, ".rds", sep = ""))