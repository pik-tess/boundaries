################################################################################
###          Status of 4 terrestrial planetary boundaries is                ####
###                   calculated and plotted                                ####
################################################################################

# based on scripts from Vera Heck for Gerten et al. 2020

rm(list = ls(all = TRUE))
olson = TRUE # use olson biome classification or LPJmL derived?

#========== libraries ==========================================================
load(RColorBrewer)
library(raster)
library(lpjmliotools) # LPJmL package developed by Fabian
library(abind)
#library(ncdf4)
#library(rgdal)

#========== directories ========================================================
dir <- c("/p/projects/open/Johanna/NETPS_within_PBs/")
dir_lpjml_out <-  paste0(dir,"lpjml_runs/output/")
r_out <- paste0(dir, "R/r_out/r_data/")
if (olson == TRUE) {
  plotpath <- paste0(dir, "R/r_out/", "plots/BECCS/7_PB_status/olson/")
} else {
  plotpath <- paste0(dir, "R/r_out/", "plots/BECCS/7_PB_status/lpjml/")
}
lu_input_lpjml <- paste0("/p/projects/open/Jannes/inputgen/input_CLM2/",
                        "cftfrac_fullhist_1500-2017_64bands_f2o.clm")

#========== set variables ======================================================
ncell <- 67420
timeframe <- c(2015, 2054)
firstyear <- timeframe[1]
lastyear <- timeframe[2]
nyears <- lastyear - firstyear + 1
ref <- 2015 # first year of simulation
landuse_year <- 2015 # year of constant land use scenario

#========== load relevent data from other scripts ==============================
# for plot settings
load(paste0("/p/projects/open/Johanna/Glob_Sus_EFRs_and_Dietary_Changes/r_out/",
            "plot_settings_Robinson.RDATA"))
# for PB Water: efr_safe, efr_uncertainty, mdisch_cell
load(file = paste0(r_out, "ann_efr_vmf_", firstyear, "_", lastyear, ".RDATA"))
if (olson == TRUE) {
  # for PB BI: current_LU_derived_BII
  load(file = paste0(r_out, "BII_data_newbold_olson.RDATA"))
  # for PB LSC:
  load(file = paste0(r_out, firstyear, "-", lastyear, "_biomes_pnv_olson",
        ".RDATA"))
} else {
  # for PB BI: current_LU_derived_BII
  load(file = paste0(r_out, "BII_data_newbold_lpjml.RDATA"))
  # for PB LSC:
  load(file = paste0(r_out, firstyear, "-", lastyear, "_biomes_pnv_lpjml",
        ".RDATA"))
}
# for PB N: N_runoff_year_lu
load(file = paste(r_out, "data_PB_N.RDATA", sep = ""))

#========== PB Water ===========================================================
# think about whether to change PB water script to also save discharge / efrs in
# hm3 / day --> better comparable to synthesis paper / slighlty changes the
# results if calculations are done instead with km3/month as done here
days_per_months <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

disch_baseline <- mdisch_cell # PNV discharge
# 1 hm3/d = 100^3 m3/ 86400 s = 11.574 m3/s
# in order dismiss flows < 1 m3/s, translate to 0.0864
cut <- 0.0864 # cut at 1 m3/s
# convert to km3/months
cut <- cut * mean(days_per_months) / 10^3

#---------------read in lu discharge--------------------------------------------
discharge_lu <- readMonthly(inFile = paste0(dir_lpjml_out,
                             "lu/mdischarge.bin"),
                             startyear = firstyear, stopyear = lastyear,
                             size = 4, headersize = 51, getyearstart = firstyear,
                             getyearstop = lastyear)
disch_lu_hm3day <- apply(discharge_lu, c(1, 2), mean)
# unit: hm3/day, conversion to km3/months
for (m in 1:12) discharge_lu[, m, ] <- discharge_lu[, m, ] *
                                      days_per_months[m] * 10^-3

# mean discharge_lu in km3/months (not needed for EFR calculation, but later)
disch_lu <- apply(discharge_lu, c(1, 2), mean)

# calculation of EFR transgressions = EFR deficits in LU run ----------------- #

def_safe <- efr_safe - disch_lu

# calculation of uncertainty zone -------------------------------------------- #
span <- efr_safe - efr_uncertainty

# dismiss boundary status calculation if PNV discharge is < 0.1 m^3/s -------- #
span[disch_baseline < cut] <- 0
def_safe[def_safe < cut] <- 0
def_safe[disch_baseline < cut] <- 0

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
pb_freshwater[pb_freshwater > pb_range[2]] <- pb_range[2]
pb_freshwater[pb_freshwater < pb_range[1]] <- pb_range[1]
pb_freshwater[cells_maginal_disch < 0] <- (-1)

#========== PB Land system Change ==============================================
# read CFT input
#dir_LPJmL_LU <- paste0("/p/projects/lpjml/input/historical/input_VERSION2/")
#cft_all <- readCFTinput(inFile = paste0(dir_LPJmL_LU,
#                        "cft1700_2005_irrigation_systems_64bands.bin"),
#                        startyear = 1700, stopyear = 2005, bands = 64, size = 2,
#                        headersize = 43, dtype = "integer",
#                        landuse_year, landuse_year) / 1000

ref <- 1500 # first year of lu input
  end <- 2017 # last year of lu input
  headersize <- 51
  size <- 4
  dtype <- "double"
  divide <- 1
  year <- 2015
  ncfts <- 16
nirr <- 4 # number of bands per cfts
nbands <- ncfts * nirr # bands in lu input

cft_all <- readCFTinput(inFile = lu_input_lpjml, startyear = ref,
            stopyear = end, bands = nbands,  size = size, dtype = dtype,
            headersize = headersize, getyearstart = year, getyearstop = year) /
            divide

lu_sum <- apply(cft_all, c(2), sum) # all cft bands summed up

biomes <- c("trop_forest_mask_oceania", "trop_forest_mask_africa",
            "trop_forest_mask_america", "trop_forest_mask_asia",
            "temp_forest_mask_northam", "temp_forest_mask_southam",
            "temp_forest_mask_europe",
            "temp_forest_mask_africa",  "temp_forest_mask_asia",
            "temp_forest_mask_oceania", "boreal_forest_mask_america",
            "boreal_forest_mask_eurasia")

lu_def <- array(NA, length(biomes))
LSC <- array(NA, ncell)
transgressed_biomes <- vector()
for (i in 1:length(biomes)) { #loop over biomes
  ii <- which(get(biomes[i])  == 1) #cells that belong to forest biome i
  forest_area <- get(biomes[i]) * cellarea #forest area of biome i
  deforestation <- sum(lu_sum * forest_area) / sum(forest_area)
  assign(paste0("state_", biomes[i]), 100 - round(deforestation * 100, 0))
  if (i <= 4 | i >= 11)  { #tropical and boreal
    if (deforestation > 0.400001) {
      LSC[ii] <-  2.5
      transgressed_biomes <- abind(transgressed_biomes, biomes[i])
    } else if (deforestation <= 0.400001 &  deforestation > 0.15) {
      LSC[ii]  <- 1.5
      transgressed_biomes <- abind(transgressed_biomes,biomes[i] )
    } else if (deforestation <= 0.15) {
      LSC[ii] <- 0.5
    }
  } else if (i >= 5 & i <= 10) { #temperate
    if (deforestation > 0.7) {
      LSC[ii] <- 2.5
      transgressed_biomes <- abind(transgressed_biomes, biomes[i])
    } else if (deforestation <= 0.7 &  deforestation > 0.5) {
      LSC[ii] <- 1.5
      transgressed_biomes <- abind(transgressed_biomes, biomes[i])
    } else if (deforestation <= 0.5) {
      LSC[ii] <- 0.5
    }
  }
  lu_def[i] <- deforestation
}

# 50% forest according to Ingo  40% frgmented forest, 60% dense forest
# deforested_cells <- which((lu_sum * forest_mask) > 0.5)
# LSC[deforested_cells] <- 3.5

#========== PB N ===============================================================

# N concentration in runoff in LU run
# conversion from g/ l to mg/l

N_runoff_year_lu <- ifelse(runoff_year_lu != 0,
                          N_runoff_year_lu / runoff_year_lu * 1000,
                          0)

# N concentration in runoff in PNV run with constant depostion 
# first 30 years of input (to include transgressions due to higher
# N depostion)
N_runoff_year_PNV <- ifelse(runoff_year_pnv_constdep != 0,
                          N_runoff_year_pnv_constdep /
                          runoff_year_pnv_constdep * 1000, 0)

# test - better results
N_runoff_year_lu <- ifelse((N_runoff_year_lu -
                                  N_runoff_year_PNV) > 0,
                                  N_runoff_year_lu -
                                  N_runoff_year_PNV, 0)

# test2
# nur in cells mit Überschreitungen in PNV: Abzug der Überschreitung??
#PNV_trans <- numeric(length = ncell)
#for (i in 1:ncell) {
#  if (N_runoff_year_PNV[i] > 1) {
#    # only the amount that is above the threshold
#    PNV_trans[i] <- N_runoff_year_PNV[i] - 1
#    N_runoff_year_lu[i] <- ifelse((N_runoff_year_lu[i] - PNV_trans[i]) > 0,
#                              N_runoff_year_lu[i] - N_runoff_year_PNV[i], 0)
#  }
#}

# mask for dry regions
mask <- mask_n
plot.mask <- toraster_rob(mask)

#========== save ===============================================================
# transgression status for each boundary: 0.5 = safe, 1.5 = uncertainty,
# 2.5 = transgressed, -0.5, boundary definition does not apply
LSC_data <- LSC
LSC_data[is.na(LSC_data)] <- (-0.5)

N_data <- BI_data <- W_data <- rep(NA, ncell)
N_data[N_runoff_year_lu >= 2.5] <- 2.5
N_data[N_runoff_year_lu >= 1 & N_runoff_year_lu < 2.5] <- 1.5
N_data[N_runoff_year_lu < 1] <- 0.5
N_data[mask == 0] <- (-0.5)

BI_data[mean_bii_biomes < 30] <- 2.5
BI_data[mean_bii_biomes >= 30 & mean_bii_biomes < 90] <- 1.5
BI_data[mean_bii_biomes >= 90] <- 0.5

W_data[pb_freshwater < 0] <- (-0.5)
W_data[pb_freshwater >= 0 & pb_freshwater < 5] <- 0.5
W_data[pb_freshwater >= 5 & pb_freshwater < 75] <- 1.5
W_data[pb_freshwater >= 75] <- 2.5
W_data[is.na(pb_freshwater)] <- 0.5

# to remove cells without lu_frac
W_data[is.na(BI_data)] <- N_data[is.na(BI_data)] <- 
                                                LSC_data[is.na(BI_data)] <- NA


# matrix with all transgression stati
PB_stati <- rbind(LSC_data, N_data, BI_data, W_data)
rownames(PB_stati) <- c("Land System Change", "Nitrogen", "Biosphere Integrity",
                        "Water")
# vector with number of boundary transgressions
pb_transgressions <- pb_transgressions_wo_water <- rep(0, ncell)
for (i in 1:ncell) {
  pb_transgressions[i] <- length(which(PB_stati[, i] == 1.5 |
                                         PB_stati[, i] == 2.5))
  pb_transgressions_wo_water[i] <- length(which(PB_stati[1:3, i] == 1.5 |
                                         PB_stati[1:3, i] == 2.5))
}
save(PB_stati, LSC_data, N_data, BI_data, W_data, pb_transgressions,
        pb_transgressions_wo_water, file = paste0(r_out, "pb_stati.RDATA"))


#========== plotting ===========================================================
# a boundary status
# b number of transgressed boundaries 

# a boundary status ------------------------------------------------------------
mask <- rep(0, 67420)
plot.nat <- toraster_rob(mask)

# covert to robinson raster
plotvar_w <- toraster_rob(W_data)
plotvar_b <- toraster_rob(BI_data)
plotvar_l <- toraster_rob(LSC_data)
plotvar_n <- toraster_rob(N_data)

# cols & lines
lwd.cntr <- 0.3
darkgrey <- "grey20"
myyellow <- colorRampPalette(c("yellow", "gold"), bias = 1)(3)[2]
cols <- c("grey90", transp("green4", 0.6), myyellow, "brown3",
            "darkgrey")

#pdf(paste0(plotpath, "PB_status_", landuse_year, ".pdf"),
#                            width = 17/2.54, height = 7.5/2.54, pointsize = 7)
png(paste0(plotpath, "PB_status_", landuse_year, "+.png"), width = 17,
                height = 7.5, units = "cm", res = 600, pointsize = 7)
{layout(matrix(c(1:4), 2, 2, byrow = TRUE))
 textcex <- 1
 par(mar = rep(0, 4), xpd = T)
 brk <- c(-1:4)
 # PB B ======
#image(plot.nat, zlim = c(-1, 0), asp = 1, xaxt = "n", yaxt = "n", xlab = "",
#            ylab = "", col = cols[2], lwd = 0.1, bty = "n")
 image(plotvar_b, asp = 1, xaxt = "n", yaxt = "n",
            xlab = "", ylab = "", col = cols, breaks = brk, lwd = 0.1,
            bty = "n")
 plot(country.outline.rob, add = TRUE, lwd = lwd.cntr,
                            border = transp(darkgrey, 0.4), usePolypath = FALSE)
 mtext("a", 3, -2.7, cex = textcex, font = 1, adj = 0.1)
 mtext("Biosphere integrity", 1, -3.5, cex = textcex, font = 1, adj = 0.57)
 
 # PB Land ======
 #image(plot.nat, zlim = c(-1, 0), asp = 1, xaxt = "n", yaxt = "n", xlab = "",
#            ylab = "", col = cols[1], lwd = 0.1, bty = "n")
 image(plotvar_l, asp = 1, xaxt = "n", yaxt = "n",
            xlab = "", ylab = "", col = cols, breaks = brk, lwd = 0.1,
            bty = "n")
 plot(country.outline.rob, add = TRUE, lwd = lwd.cntr,
                            border = transp(darkgrey, 0.4), usePolypath = FALSE)
 mtext("b", 3, -2.7, cex = textcex, font = 1, adj = 0.1)
 mtext("Land system change", 1, -3.5, cex = textcex, font = 1, adj = 0.57)

# PB W ======
 image(plotvar_w, asp = 1, xaxt = "n", yaxt = "n",
            xlab = "", ylab = "", col = cols, breaks = brk, lwd = 0.1,
            bty = "n")
 plot(country.outline.rob, add = TRUE, lwd = lwd.cntr,
                            border = transp(darkgrey, 0.4), usePolypath = FALSE)
 mtext("c", 3, -2.7, cex = textcex, font = 1, adj = 0.1)
 mtext("Freshwater use", 1, -3.5, cex = textcex, font = 1, adj = 0.57)

# PB N=======
 image(plotvar_n, asp = 1, xaxt = "n", yaxt = "n",
            xlab = "", ylab = "", col = cols, breaks = brk, lwd = 0.1,
            bty = "n")
 #image(plot.mask, zlim = c(-1, 0), asp = 1, xaxt = "n", yaxt = "n", xlab = "",
 #           ylab = "", col = cols[1], lwd = 0.1, bty = "n", add = TRUE)
 plot(country.outline.rob, add = TRUE, lwd = lwd.cntr,
                            border = transp(darkgrey, 0.4), usePolypath = FALSE)
 mtext("d", 3, -2.7, cex = textcex, font = 1, adj = 0.1)
 mtext("Nitrogen flows", 1, -3.5, cex = textcex, font = 1, adj = 0.57)

# Legend
 par(fig = c(0.4, 0.49, 0.1, 0.17), new = TRUE)
 legend(x = -180, y = -30, cex = 1,
       legend = c("safe zone", "increasing risk", "high risk"), horiz = F,
       pch = 22, pt.bg = cols[2:4], pt.lwd = 0.0, pt.cex = 1.6, box.col = NA,
       xpd = NA, inset = -0.1)

 dev.off()
}

# b number of boundary transgressions ------------------------------------------

pb_trans <- toraster_rob(pb_transgressions)
brk <- seq(-0.5, 4.5, 1)
cols <- c("white", "#ffae00da", "#ff5100", "#ad0808", "black")

tiff(paste(plotpath, "n_PB_transgressions_", landuse_year, ".tiff", sep = ""),
  width = 17, height = 11, res = 600, unit = "cm", compression = "lzw",
  family = "sans")
  textcex <- 1
  par(oma <- c(0, 0, 0, 0), xpd = T)
  par.settings <- list(fontsize = list(text = 9))

  image(pb_trans, zlim = brk, ylim = ylim, xlim = c(xlim[1] * 0.8,
      xlim[2] * 0.92), asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
      col = cols, breaks = brk, lwd = 0.1, bty = "n")
  plot(country.outline.rob, add = TRUE, lwd = lwd.cntr + 0.2,
                            border = transp(darkgrey, 0.4), usePolypath = FALSE)
  mtext("Number of PB transgressions", 1, -14.5, cex = textcex, font = 1,
      adj = 0.45)

  legend(x = -6221874, y = -8046880, cex = 0.9,
             legend = c("1", "2", "3", "4"), horiz = F, ncol = 4,
             pch = 22, pt.bg = c(cols[2:5]), pt.lwd = 0.0, pt.cex = 1.6,
             box.col = NA, xpd = NA, inset = -0.1,
             text.width = rep(14394230 / 6, 4))

  dev.off()

# plot number of PB transgressions without water
pb_trans_wo_water <- toraster_rob(pb_transgressions_wo_water)
tiff(paste(plotpath, "n_PB_transgressions_wo_water_", landuse_year, ".tiff", sep = ""),
  width = 17, height = 11, res = 600, unit = "cm", compression = "lzw",
  family = "sans")
  textcex <- 1
  par(oma <- c(0, 0, 0, 0), xpd = T)
  par.settings <- list(fontsize = list(text = 9))

  image(pb_trans_wo_water, zlim = brk, ylim = ylim, xlim = c(xlim[1] * 0.8,
      xlim[2] * 0.92), asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
      col = cols, breaks = brk, lwd = 0.1, bty = "n")
  plot(country.outline.rob, add = TRUE, lwd = lwd.cntr + 0.2,
                            border = transp(darkgrey, 0.4), usePolypath = FALSE)
  mtext("Number of PB transgressions without water", 1, -14.5, cex = textcex,
      font = 1, adj = 0.45)

  legend(x = -6221874, y = -8046880, cex = 0.9,
             legend = c("1", "2", "3", "4"), horiz = F, ncol = 4,
             pch = 22, pt.bg = c(cols[2:5]), pt.lwd = 0.0, pt.cex = 1.6,
             box.col = NA, xpd = NA, inset = -0.1,
             text.width = rep(14394230 / 6, 4))

dev.off()

# plot areas without PB transgressions (without water)

no_pb_transgressions <- rep(1, ncell)
no_pb_transgressions[which(pb_transgressions_wo_water != 0)] <- 0
no_pb_trans <- toraster_rob(no_pb_transgressions)
no_pb_transgressions_inkl_w <- rep(1, ncell)
no_pb_transgressions_inkl_w[which(pb_transgressions != 0)] <- 0
no_pb_trans_inkl_w <- toraster_rob(no_pb_transgressions_inkl_w)


brk <- c(-0.5, 0.5, 1.5)
cols <- rev(c("#85e285", "black"))

tiff(paste(plotpath, "no_PB_transgressions_wo_water_", landuse_year, ".tiff", sep = ""),
  width = 17, height = 11, res = 600, unit = "cm", compression = "lzw",
  family = "sans")
  textcex <- 1
  par(oma <- c(0, 0, 0, 0), xpd = T)
  par.settings <- list(fontsize = list(text = 9))

  image(no_pb_trans, zlim = brk, ylim = ylim, xlim = c(xlim[1] * 0.8,
      xlim[2] * 0.92), asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
      col = cols, breaks = brk, lwd = 0.1, bty = "n")
  plot(country.outline.rob, add = TRUE, lwd = lwd.cntr + 0.2,
                            border = transp(darkgrey, 0.4), usePolypath = FALSE)
  mtext("Number of PB transgressions without water", 1, -14.5, cex = textcex,
      font = 1, adj = 0.45)

  legend(x = -6221874, y = -8046880, cex = 0.9,
             legend = c("no PB transgressions"), horiz = F, ncol = 4,
             pch = 22, pt.bg = cols[1], pt.lwd = 0.0, pt.cex = 1.6,
             box.col = NA, xpd = NA, inset = -0.1,
             text.width = 14394230 / 6)

dev.off()


# plot area availability for plantations considering PB transgressions

load(file = paste0(r_out, "pot_fracs", ".RDATA"))
bt_irr_lf <- toraster_rob(bt_irr_lowFert_pot_frac)
bt_irr_hf <- toraster_rob(bt_irr_highFert_pot_frac)
bt_rf_lf <- toraster_rob(bt_rf_lowFert_pot_frac)
bt_rf_hf <- toraster_rob(bt_rf_highFert_pot_frac)
bg_irr_lf <- toraster_rob(bg_irr_lowFert_pot_frac)
bg_irr_hf <- toraster_rob(bg_irr_highFert_pot_frac)
bg_rf_lf <- toraster_rob(bg_rf_lowFert_pot_frac)
bg_rf_hf <- toraster_rob(bg_rf_highFert_pot_frac)

# multiply area availability by PB mask

# a - seperate for plantation configurations
 
 # for rainfed with 3 PBs included, for irrigated with 4 PBs included
library(scales)
BFT_runs_rf <- c("bt_rf_lf", "bt_rf_hf", "bg_rf_lf", "bg_rf_hf")
BFT_runs_irr <- c("bt_irr_lf", "bt_irr_hf", "bg_irr_lf", "bg_irr_hf")
BFT_runs <- c("bg_rf_hf", "bg_rf_lf", "bg_irr_hf", "bg_irr_lf",
              "bt_rf_hf", "bt_rf_lf", "bt_irr_hf", "bt_irr_lf")

for (bioname in BFT_runs_rf) {
  assign(paste0("plotvar_", bioname), get(bioname) * no_pb_trans)
}
for (bioname in BFT_runs_irr) {
  assign(paste0("plotvar_", bioname), get(bioname) * no_pb_trans_inkl_w)
}
png(paste0(plotpath, "BFT_areas_for_optimization.png"), width = 20, height = 7.5,
                                         units = "cm", res = 600, pointsize = 7)
{layout(matrix(c(1:8), 2, 4, byrow = F))
  textcex <- 1
  par(mar = rep(0, 4), xpd = T)
  brk <- seq(0, 1, 0.1)
  cols <- alpha("#1ab65b", seq(0, 1, length = 10))
  for (bioname in BFT_runs) {
  image(get(paste0("plotvar_", bioname)), ylim = ylim, xlim = xlim, asp = 1, xaxt = "n",
        yaxt = "n", xlab = "", ylab = "", col = cols, breaks = brk, lwd = 0.1,
        bty = "n")
  plot(country.outline.rob, add = TRUE, lwd = lwd.cntr,
        border = "grey", usePolypath = FALSE)
  mtext(bioname, 1, -4.0, cex = textcex, font = 1, adj = 0.57)
  }
 dev.off()
}

# b - one combined map, always including the maximumg fraction from all BP masks


max_frac <- mosaic(plotvar_bt_irr_lf, plotvar_bt_irr_hf, plotvar_bt_rf_lf,
                    plotvar_bt_rf_hf, plotvar_bg_irr_lf, plotvar_bg_irr_hf,
                    plotvar_bg_rf_lf, plotvar_bg_rf_hf, fun = "max")

png(paste0(plotpath, "areas_for_opti_overlay.png"), width = 15.5, height = 7.5,
                                         units = "cm", res = 600, pointsize = 7)
{textcex <- 1
  par(mar = rep(1, 4), xpd = T)
  brk <- seq(0, 1, 0.1)
  cols <- alpha("#1ab65b", seq(0, 1, length = 10))
  image(max_frac, ylim = ylim, xlim = xlim, asp = 1, xaxt = "n",
        yaxt = "n", xlab = "", ylab = "", col = cols, breaks = brk, lwd = 0.1,
        bty = "n")
  plot(country.outline.rob, add = TRUE, lwd = lwd.cntr,
        border = "grey", usePolypath = FALSE)
dev.off()
}