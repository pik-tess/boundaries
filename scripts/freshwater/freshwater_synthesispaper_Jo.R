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

# remove land use change induced EFR deficits to only calculate EFR deficits
# induced by excessive water withdrawals?
# if this option is enabled, the discharge output from a lu run without
# water withdrawals is needed
remove_luc_efr_deficits <- FALSE

# remove climate change induced EFR deficits?
# enabling this option only has an effect if the reference PNV discharge for
# EFR calcultion was simulated with a different climate than the timeframe of
# interest / the lu run
remove_cc_efr_deficits <- FALSE

# paths ---------------------------------------------------------------------- #
run <- "ostberg" # old
dir_inputs <- "/p/projects/open/Vera/synthesepaper/input_data/"
dir_r_out <- "/p/projects/open/Jannes/earth4all/pb_status/R/data/"
dir_lpj_out <- "/p/projects/open/Jannes/earth4all/pb_status/runs/output/"
dir_lpj_out_pnv <- paste0(dir_lpj_out, "/", run, "_pnv", "/")
dir_lpj_out_lu <- paste0(dir_lpj_out, "/", run, "_lu", "/")
# dir_lpj_out_lu_rf <- paste0(dir_lpj_out, "/", run, "_lu", "/")

# dir_lpj_out <- "/p/projects/open/Vera/synthesepaper/lpj_out/"
# dir_lpj_out_pnv <- paste0(dir_lpj_out, "1950_2009_pnv/")
# dir_lpj_out_lu <- paste0(dir_lpj_out, "1950_2009_baseline_noefr/")

# dir_lpj_out <- "/p/projects/open/Johanna/NETPS_within_PBs/BECCS_nat_veg/LPJmL_runs/"
# dir_lpj_out_pnv <- paste0(dir_lpj_out, "PNV/output/")
# dir_lpj_out_lu <- paste0(dir_lpj_out, "LU_2005/output/")

# arguments & params --------------------------------------------------------- #
ncell <- 67420
timespan <- c(1988, 2017) #c(1980, 2009)
ref <- 1901 #1950 #1901 # first year of simulation
firstyear <- timespan[1]
nyears <- timespan[2] - timespan[1] + 1
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
     where = (timespan[1] - ref) * temp_resolution * nbands * ncell * size,
     origin = "start")
disch_baseline_allyears <- readBin(file,
                         double(),
                         n = ncell * temp_resolution * nyears * nbands,
                         size = size)
close(file)
dim(disch_baseline_allyears) <- c(ncell, temp_resolution, nyears)

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
     where = (timespan[1] - ref) * temp_resolution * nbands * ncell * size,
     origin = "start")
disch_lu_allyears <- readBin(file,
                             double(),
                             n = ncell * temp_resolution * nyears * nbands,
                             size = size)
close(file)
dim(disch_lu_allyears) <- c(ncell, temp_resolution, nyears) #

disch_lu <- apply(disch_lu_allyears, c(1, 2), mean) # mean over years


# read in discharge landuse rainfed ------------------------------------------ #
if (remove_luc_efr_deficits == TRUE) {
  file <- file(paste(dir_lpj_out_lu_rf, "mdischarge.bin", sep = ""), "rb")
  seek(file,
       where = (timespan[1] - ref) * temp_resolution * nbands * ncell * size,
       origin = "start")
  disch_lu__rf_allyears <- readBin(file,
                               double(),
                               n = ncell * temp_resolution * nyears * nbands,
                               size = size)
  close(file)
  dim(disch_lu__rf_allyears) <- c(ncell, temp_resolution, nyears) #

  disch_lu_rf <- apply(disch_lu__rf_allyears, c(1, 2), mean) # mean over years
}

# calculation of EFR transgressions = EFR deficits --------------------------- #


# --- option: take out climate induced EFR deficits -------------------------- #

# set EFRs to PNV discharge if EFRs are below PNV discharge
# only makes sense if PNV derived EFRs were calculated with a different climate
# --> this would then take out climate change induced transgressions to some
# degree

if (remove_cc_efr_deficits == TRUE) {
  for (month in 1:12) {
  # vmfmax (vmf + 15%)
    efr_req_safe[, month] <- ifelse(
     efr_req_safe[, month] > disch_baseline[, month],
     disch_baseline[, month],
     efr_req_safe[, month]
   )
   # vmfmin (vmf - 15%)
   efr_req_uncertainty[, month] <- ifelse(
     efr_req_uncertainty[, month] > disch_baseline[, month],
     disch_baseline[, month],
     efr_req_uncertainty[, month]
   )
 }
}

# calculate EFR deficits in LU run ------------------------------------------- #

def_safe <- efr_req_safe - disch_lu

# if climate induced transgression are to be taken out (see above):
# substract amount of EFR transgression in PNV run
if (remove_cc_efr_deficits == TRUE) {
  def_pnv_safe <- efr_req_safe - disch_baseline
  # if disch(PNV) > efr_req_safe
  def_pnv_safe[def_pnv_safe < 0] <- 0
  # all cells where PNV discharge is smaller than efr requirements have a value
  def_pnv_safe[disch_baseline < cut] <- 0
  def_safe <- def_safe - def_pnv_safe
}

# if land use change induced EFR deficits are to be taken out:
# subtracts EFR deficits from rainfed baseline scenario
if (remove_luc_efr_deficits == TRUE) {
  def_lu_rf <- efr_req_safe - disch_lu_rf
  for (i in 1:67420) {
    for (j in 1:12) {
      if (def_lu_rf[i, j] > 0) {
        def_safe[i, j] <- def_safe[i, j] - def_lu_rf[i, j]
      }
    }
  }
}

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
        file = paste(dir_r_out, "pb_freshwater.rds", sep = ""))


# test

#load("/p/projects/open/Johanna/Glob_Sus_EFRs_and_Dietary_Changes/r_out/plot_settings_Robinson.RDATA")

#range<- c(0,100)
#    brk<- c(-2,0,5,75,100)
#
#    plotvar2 <- toraster_rob(pb_freshwater)
    
#    cols <- c("grey90", "yellowgreen", "yellow1", "red3","darkgrey")
#    mask = array(0,ncell)
#    plot.nat <- toraster_rob(mask)
#
#plotpath <- paste0("/p/projects/open/Johanna/NETPS_within_PBs/BECCS_nat_veg/", 
#        "R_out/", "plots/")
#png(paste(plotpath,"PB_water_status_LPJmL5.1.png",sep=""),width = 8.8, height = 7, pointsize=7, res=600, unit="cm")
#
#      textcex=1
#      par(mar=c(2.5,0,1.5,1.2),xpd=T)
#      
#      image(plot.nat,zlim=c(-1,0),asp=1,xaxt="n",yaxt="n",xlab="",ylab="",col=cols[2],lwd=0.1,bty="n")
#      image(plotvar2,ylim=ylim,xlim=xlim,asp=1,xaxt="n",yaxt="n",xlab="",ylab="",col=cols[1:4],breaks=brk,lwd=0.1,bty="n",add=TRUE)
#      plot(country.outline.rob,add=TRUE,lwd=lwd.cntr, border=transp("#888383",0.7),usePolypath=FALSE)
#      
#      legend(x=-15221874,y=-8046880,cex=0.9,
#             legend=c("safe zone","increasing risk","high risk", "no EFR deficit calculation"),horiz=F, ncol=4,
#             pch=22,pt.bg=c(cols[2:4],"grey90"),pt.lwd=0.0,pt.cex=1.6,box.col=NA,xpd=NA,inset=-0.1,text.width = c(14394230/2.5,14394230/2.6,14394230/2.35,14394230/2.5))
#      dev.off()