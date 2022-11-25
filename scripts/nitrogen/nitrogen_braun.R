# calculate nitrogen PB status based on N in leaching output from LPJmL

rm(list = ls(all = TRUE))
gc()

library(lpjmliotools)
library(raster)
library(RColorBrewer)
library(maptools)
library(fields)

reference_time <- c(2007, 2017) # time period that is considered for PB N calculation

reference_time <- c(2007, 2017) # time period that is considered for PB N calculation

last_year <- 2017
#### directories ####
 run <- "ostberg" # old
# dir_lpj_out <- paste0(
#   "/p/projects/open/Jannes/earth4all/pb_status/runs/output/",
#   run,
#   "_pnv/"
# )
# dir_lpj_out <-"/p/projects/open/Jannes/lpjml/testing/nitrogen_fixing/output/nfixation_fix_pnv_depconst/"
dir_lpj_out <-"/p/projects/open/Jannes/lpjml/testing/nitrogen_fixing/output/nfixation_fixmedian_lu/"

plotpath <- "/p/projects/open/Jannes/earth4all/pb_status/R/plots/"
dir_out <- "/p/projects/open/Jannes/earth4all/pb_status/R/data/"

last_year <- reference_time[2]
#### read in relevant data ####

#------ grid---------
ncell <- 67420
zz <- file(paste(dir_lpj_out, "/grid.bin", sep = ""), "rb")
x <- readBin(zz, integer(), n = 2 * 67420, size = 2) / 100
lon <- x[c(1:ncell) * 2 - 1]
lat <- x[c(1:ncell) * 2]
cellarea2 <- (111e3 * 0.5) * (111e3 * 0.5) * cos(lat / 180 * pi)
close(zz)

#------ runoff-----

runoff <- readMonthly(
  inFile = paste0(dir_lpj_out, "mrunoff.bin", sep = ""),
  startyear = 1901, stopyear = last_year, size = 4,
  headersize = 0, getyearstart = reference_time[1],
  getyearstop = reference_time[2]
)
# unit: mm/month; mm = 1 l/m2

#------ leaching--------
leaching <- readMonthly(
  inFile = paste(dir_lpj_out, "mleaching.bin", sep = ""),
  startyear = 1901, stopyear = last_year, size = 4,
  headersize = 0, getyearstart = reference_time[1],
  getyearstop = reference_time[2]
)
# unit: gN/m2/month
#------ precipitation and pet ------
# no output for the actual run. Taken from another run: p/projects/open/Johanna/NETPS_within_PBs/BECCS_nat_veg/LPJmL_runs/PNV/
mprec <- readMonthly(
  inFile = paste0(dir_lpj_out, "/mprec.bin"),
  startyear = 1979, stopyear = last_year, size = 4, headersize = 0,
  getyearstart = reference_time[1], getyearstop = reference_time[2]
)

mpet <- readMonthly(
  inFile = paste0(dir_lpj_out, "/mpet.bin"),
  startyear = 1979, stopyear = 2018, size = 4, headersize = 0,
  getyearstart = reference_time[1], getyearstop = reference_time[2]
)

#### denitrification factor from Bouwman et al. 2013 #####

# account for denitrification in groundwater and rpiarian zones based on Bouwman et al. 2013

# if long term storage in groundwater is neglected and denitrification from groundwater is assumed not to be covered within LPJmL : (65+9)/(93+11) = 0.71 (used in Gerten et al. 2020)
# if long term storage in groundwater is considered  and denitrification from groundwater is assumed not to be covered within LPJmL : (65)/(93+11) = 0.63
# if long term storage in groundwater is neglected and denitrification is assumed to be covered in LPJmL: (65+9)/(49+5+15+11) = 0.93
# if long term storage in groundwater is considered and denitrification is assumed to be covered in LPJmL: 65/(49+5+15+11) = 0.81

neglect_groundwater_storage <- T
include_groundwater_denit <- T

if (neglect_groundwater_storage == T) {
  numerator <- 65 + 9
} else {
  numerator <- 65
}

if (include_groundwater_denit == T) {
  denominator <- 93 + 11
} else {
  denominator <- 49 + 5 + 15 + 11
}

groundwater_and_denit_factor <- numerator / denominator

#### calculate N concentration in runoff to surface water ####

#------ per month-----

days_per_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

mean_runoff <- apply(runoff, c(1, 2), mean) * cellarea2
# --> has to be converted to l (mm/month = l/m2 --> no conversion necessary)
# cellarea

mean_leaching <- apply(leaching, c(1, 2), mean) * cellarea2 * 1000
# has to be converted to mg N (from gN/m2)
# * cellarea * 1000

N_in_surfaceW <- ifelse(mean_runoff != 0, (mean_leaching * groundwater_and_denit_factor) / mean_runoff, 0)
range(N_in_surfaceW, na.rm = T)


#------ per year-----
mean_prec <- apply(apply(mprec, c(1, 3), sum), c(1), mean)
mean_pet <- apply(apply(mpet, c(1, 3), sum), c(1), mean)

mask_fert <- numeric(67420L)
mask_fert <- ifelse(mean_prec / (mean_pet + 0.000001) < 0.15, 0, 1) #

mean_runoff_year <- apply(mean_runoff, 1, sum)
mean_leaching_year <- apply(mean_leaching, 1, sum)
N_in_surfaceW_year <- ifelse(mean_runoff_year != 0, (mean_leaching_year) / mean_runoff_year, 0)
# N_in_surfaceW_year <- ifelse(mean_runoff_year!=0,(mean_leaching_year*groundwater_and_denit_factor)/mean_runoff_year,0)

N_in_surfaceW_year[mask_fert != 1] <- NA
range(N_in_surfaceW_year, na.rm = T)
save(N_in_surfaceW_year, file = paste(dir_out, "N_surface_water_pnv.RDATA", sep = ""))


#### plot settings (adapted from Vera's scripts for the synthesis paper)####

load("/p/projects/open/Johanna/Glob_Sus_EFRs_and_Dietary_Changes/r_out/plot_settings_Robinson.RDATA")


greens <- brewer.pal(n = 9, name = "BuGn")[c(4, 5, 6, 8, 9)]
blues <- brewer.pal(n = 9, name = "Blues")[c(3, 5, 6, 7)]
reds <- brewer.pal(n = 9, name = "YlOrRd")[c(2, 6, 8, 9)]
for (i in 1:length(greens)) greens[i] <- mix(greens[i], blues[i], 6, 3)
neutral <- mix(mix(mix(reds[1], greens[1], 3, 2), "grey40", 3, 2), "white", 6, 4)

#### mask cells for PB Status ####

# Vera for synthesis Paper:
# thresh <-  which(ann_prec/(ann_pet +0.000001 ) < 0.15)
# (prec = precipitation, pet = potential evapotranspiration; both outputs not available in Werners folder)

# Other options:
# Steffen et al: gray areas are areas where no fertilizer is applied

# option 1 based on fertilization rates (does not make too much sense)
# option 2 based on CFT frac
# option 3 based on fertilization amount




#### yearly average N concentration --> PB status ####

# thresholds from de Vries et al. 2013: critical N concentration in surface water: 1-2.5 mg/l (PB: 1 mg/l; upper uncertainty: 2.5 mg/l)

## prepare plot
range <- c(-1, 4)

var <- N_in_surfaceW_year
var[var > 4] <- 4
var[var == 0] <- 0.001 # cells with runoff = 0 are displayed gray (maybe improve: cells with leaching = 0 also displayed gray) --> dsicarded, therfore set to 0.001 so that 0 is not displayed grey
plotvar1 <- toraster_rob(var)

brk <- c(-1, 0, 1, 2.5, 4)
myyellow <- colorRampPalette(c("yellow", "gold"), bias = 1)(3)[2]
cols <- c("grey90", transp("green4", 0.6), myyellow, "brown3") # brown3, firebrick3
# cols <- c("yellowgreen", "yellow1", "red3")


tiff(paste(plotpath, "PB_N_status_nfix.tiff", sep = ""), width = 17, height = 11, res = 600, unit = "cm", compression = "lzw", family = "sans")
textcex <- 1
par(oma = c(0, 0, 0, 0), xpd = T)
par.settings <- list(fontsize = list(text = 9))

image(plotvar1, zlim = range, ylim = ylim, xlim = c(xlim[1] * 0.8, xlim[2] * 0.92), asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = cols, breaks = brk, lwd = 0.1, bty = "n")
# image(plot.mask, zlim = c(-1, 0), asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", lwd = 0.1, bty = "n", add = TRUE)
plot(country.outline.rob, add = TRUE, lwd = 0.45, border = transp("grey20", 0.6), usePolypath = FALSE)

image.plot(
  legend.only = TRUE, zlim = range, nlevel = 512,
  breaks = seq(from = brk[2], to = range[2], length.out = length(brk) - 1),
  # breaks=seq(from=range[1],to=0,length.out=length(brk)/2+1),
  lab.breaks = c("0", "1", "2.5", ">"), col = cols[2:length(cols)],
  # lab.breaks=c("<",brk[2:(length(brk)/2)],"0"),col=cols[1:(length(brk)/2)],
  smallplot = c(0.35, 0.65, 0.15, 0.18), horizontal = TRUE,
  axis.args = list(cex.axis = 1, mgp = c(2.3, 0.5, 0), tcl = -0.4),
  legend.args = list(text = paste("[N] concentration in runoff in mg/l"), side = 1, adj = 0.5, cex = 1, line = 1.4)
)
dev.off()

#### yearly average N concentration in runoff not translated to PB status #####

## prepare plot
range <- c(0, 10)

N_runoff <- ifelse(mean_runoff_year != 0, (mean_leaching_year) / mean_runoff_year, 0)
range(N_runoff, na.rm = T) # N_in_surface_water

var <- N_runoff
var[var > 10] <- 10
var[var == 0] <- -1
plotvar1 <- toraster_rob(var)

brk <- c(-1, 0, 0.1, seq(1, 10))
goldx <- mix("gold", "white", 3, 2)
cols <- c("grey93", "white", colorRampPalette(c(goldx, "orange", colorRampPalette(reds)(7)[c(3, 5, 7)]))(length(brk) - 3))


tiff(paste(plotpath, "N_conc_yearly_runoff_TEST.tiff", sep = ""), width = 17, height = 11, res = 600, unit = "cm", compression = "lzw", family = "sans")
textcex <- 1
par(oma = c(0, 0, 0, 0), xpd = T)
par.settings <- list(fontsize = list(text = 9))

image(plotvar1, zlim = range, ylim = ylim, xlim = c(xlim[1] * 0.8, xlim[2] * 0.92), asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = cols, breaks = brk, lwd = 0.1, bty = "n")
image(plot.mask, zlim = c(-1, 0), asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", lwd = 0.1, bty = "n", add = TRUE)
plot(country.outline.rob, add = TRUE, lwd = 0.45, border = transp("grey20", 0.6), usePolypath = FALSE)

image.plot(
  legend.only = TRUE, zlim = range, nlevel = 512,
  breaks = seq(from = 0, to = range[2], length.out = length(brk) - 1),
  # breaks=seq(from=range[1],to=0,length.out=length(brk)/2+1),
  lab.breaks = c(0, 0.1, seq(1, 9, 1), ">"), col = cols[2:length(cols)],
  # lab.breaks=c("<",brk[2:(length(brk)/2)],"0"),col=cols[1:(length(brk)/2)],
  smallplot = c(0.2, 0.8, 0.15, 0.18), horizontal = TRUE,
  axis.args = list(cex.axis = 1, mgp = c(2.3, 0.5, 0), tcl = -0.4),
  legend.args = list(text = paste("[N] concentration in runoff/leaching in mg/l"), side = 1, adj = 0.5, cex = 1, line = 1.4)
)
dev.off()



#------ leaching--------
dir_lpj_out <- paste0(
  "/p/projects/open/Jannes/earth4all/pb_status/runs/output/",
  run,
  "_pnv/"
)
pnv_leaching <- readMonthly(
  inFile = paste(dir_lpj_out, "mleaching.bin", sep = ""),
  startyear = 1901, stopyear = last_year, size = 4,
  headersize = 0, getyearstart = reference_time[1],
  getyearstop = reference_time[2]
)

dir_lpj_out <- paste0(
  "/p/projects/open/Jannes/earth4all/pb_status/runs/output/",
  run,
  "_lu/"
)

lu_leaching <- readMonthly(
  inFile = paste(dir_lpj_out, "mleaching.bin", sep = ""),
  startyear = 1901, stopyear = last_year, size = 4,
  headersize = 0, getyearstart = reference_time[1],
  getyearstop = reference_time[2]
)

mean_pnv = apply(pnv_leaching,c(1,2), mean, na.rm=T)
mean_lu = apply(lu_leaching,c(1,2), mean, na.rm=T)

a_leach_pnv = apply(mean_pnv,c(1), sum, na.rm=T)
a_leach_lu = apply(mean_lu,c(1), sum, na.rm=T)

a_leach_compare = a_leach_pnv/a_leach_lu

png(paste(plotpath, "leaching_reldiff_PNV-LU_scatterplot.png", sep = ""), width = 9, height = 9, units = "cm", res = 150, pointsize = 7)
par(mar=c(6,6,2,2))
plot(a_leach_lu, a_leach_pnv,
     ylim= c(0,14),
     xlim=c(0,14), 
     ylab=expression(paste("N leaching for LPJmL 5.3 - PNV (gN ", m^-2,"",yr^-1,")")),
     xlab=expression(paste("N leaching for LPJmL 5.3 - LU (gN ", m^-2,"",yr^-1,")")))
abline(coef = c(0,1), col="red")
dev.off()

plot.ras <- raster(ncols = 720, nrows = 360)
plot.ras[cellFromXY(plot.ras, cbind(lon, lat))] <- a_leach_compare
plot.ras <- crop(plot.ras, extent(-180,180,-56,84))
plot.ras[is.nan(plot.ras) | is.nan(plot.ras)] <- NA
plot.ras[plot.ras < 0.6] <- 0.6
plot.ras[plot.ras > 1.4] <- 1.4

png(paste(plotpath, "leaching_reldiff_PNV-LU.png", sep = ""), width = 14, height = 8, units = "cm", res = 150, pointsize = 7)
plot(plot.ras, zlim=c(0.6,1.4), col=colorRampPalette(rev(c("darkblue","#d4eef7","#f5de0c","darkred")),space="rgb")(64), main="Nleaching ratio (PNV/LU)")
map(add=T,lwd=0.3)
dev.off()