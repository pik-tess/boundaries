################################################################################
###          calculates necessary data for the calculation                  ####
###      of the PB for Landsystem change (forest biome classification)      ####
################################################################################


rm(list = ls(all = TRUE))
save <- TRUE
plot <- FALSE

#========== libraries ==========================================================
library(lpjmliotools) # LPJmL package developed by Fabian
library(raster)
library(RColorBrewer)
library(maptools)
library(fields)

#========== directories ========================================================
dir <- c("/p/projects/open/Johanna/NETPS_within_PBs/")
dir_LPJmL_out <-  paste0(dir,"lpjml_runs/output/pnv/")

r_out <- paste0(dir, "R/r_out/r_data/")
plotpath <- paste0(dir, "R/r_out/", "plots/BECCS/5_PB_L/")

#========== set variables ======================================================

ncell <- 67420
timeframe <- c(2015, 2054)
firstyear <- timeframe[1]
lastyear <- timeframe[2]
nyears <- lastyear - firstyear + 1
ref <- 2015 # first year of simulation
end <- 2054 # last year of simulation
npft <- 12

#========== definition of Biomes  ==============================================
load(paste0(r_out, "biome_classes_lpjml.RDATA"))

#========== create Biome masks  ================================================

# LPJ forest biomes mask
forest_mask <- array(0, ncell)
forest_mask[biome_classes == 1 |
            biome_classes == 2 |
            biome_classes == 3 |
            biome_classes == 4 |
            biome_classes == 5 |
            biome_classes == 6 |
            biome_classes == 7 |
            biome_classes == 8] <- 1

# trop forets mask
trop_forest_mask <- array(0, ncell)
trop_forest_mask[biome_classes == 1 |
                 biome_classes == 2] <- 1

# temperate
temp_forest_mask <- array(0, ncell)
temp_forest_mask[biome_classes == 3 |
                 biome_classes == 4 |
                 biome_classes == 5 |
                 biome_classes == 6] <- 1

# boreal
boreal_forest_mask <- array(0, ncell)
boreal_forest_mask[biome_classes == 7 |
                   biome_classes == 8] <- 1

#test
#lpj.ras         <- raster(ncols = 720, nrows = 360)
#lpj.ras[cellFromXY(lpj.ras, cbind(lon, lat))] <- temp_forest_mask_america
#plot(lpj.ras)

if (plot == TRUE) {
  load("/p/projects/open/Johanna/Glob_Sus_EFRs_and_Dietary_Changes/r_out/plot_settings_Robinson.RDATA")

  var <- trop_forest_mask
  plotvar1 <- toraster_rob(var)

  var <- temp_forest_mask
  plotvar2 <- toraster_rob(var)

  var <- boreal_forest_mask
  plotvar3 <- toraster_rob(var)

  var <- array(1, ncell)
  plotvar0 <- toraster_rob(var)

  brk   <- c(0.01, 1)
  cols0 <- c("#ebebeb")
  cols1 <- c("#4f994f")
  cols2 <- c("#a6e284")
  cols3 <- c("#929cf3")


  tiff(paste0(plotpath,"forest_biomes_LPJmL5.1_1979_2018_new_biomes.tiff"),width = 17, height = 11, res=600, unit="cm",compression = 'lzw', family = "sans")
  textcex = 1
  par(oma = c(0, 0, 0, 0), xpd = T)
  par.settings = list(fontsize = list(text = 9))

  image(plotvar0,zlim=brk,ylim=ylim,xlim=c(xlim[1]*0.8,xlim[2]*0.92),asp=1,xaxt="n",yaxt="n",xlab="",ylab="",col=cols0,breaks=brk,lwd=0.1,bty="n")
  image(plotvar1,zlim=brk,ylim=ylim,xlim=c(xlim[1]*0.8,xlim[2]*0.92),asp=1,xaxt="n",yaxt="n",xlab="",ylab="",col=cols1,breaks=brk,lwd=0.1,bty="n", add=TRUE)
  image(plotvar2,zlim=brk,ylim=ylim,xlim=c(xlim[1]*0.8,xlim[2]*0.92),asp=1,xaxt="n",yaxt="n",xlab="",ylab="",col=cols2,breaks=brk,lwd=0.1,bty="n", add=TRUE)
  image(plotvar3,zlim=brk,ylim=ylim,xlim=c(xlim[1]*0.8,xlim[2]*0.92),asp=1,xaxt="n",yaxt="n",xlab="",ylab="",col=cols3,breaks=brk,lwd=0.1,bty="n", add=TRUE)

  dev.off()

  #========== plot tree cover ==================================================
   var <- fpc.tree.total
  var[is.na(var)] <- -2
  plotvar1 <- toraster_rob(var)

  brk <- c(-2, -1, 0.001, seq(0.1, 1, 0.1))
  cols1 <- c("#ebebeb", "white",
              colorRampPalette(brewer.pal(9, "YlGn"))(length(brk) - 3))



  tiff(paste(plotpath,"tree_cover_LPJmL5.1_1979_2018_LU_2015.tiff",sep=""),width = 17, height = 11, res=600, unit="cm",compression = 'lzw', family = "sans")
  textcex=1
  par(oma=c(0,0,0,0),xpd=T)
  par.settings = list(fontsize = list(text = 9))

  image(plotvar1,zlim=brk,ylim=ylim,xlim=c(xlim[1]*0.8,xlim[2]*0.92),asp=1,xaxt="n",yaxt="n",xlab="",ylab="",col=cols1,breaks=brk,lwd=0.1,bty="n")
  plot(country.outline.rob,add=TRUE,lwd=0.45, border=transp("grey20",0.6),usePolypath=FALSE)

  image.plot(legend.only=TRUE,zlim=c(0,1),nlevel=512,
           breaks=seq(0,1,0.1),
           lab.breaks=seq(0,1,0.1) * 100,col=cols1[3:length(cols1)],
           smallplot=c(0.35,0.65,0.15,0.18),horizontal=TRUE,
           axis.args=list(cex.axis=1,mgp=c(2.3,0.5,0),tcl=-0.4),
           legend.args=list(text=paste("tree cover in %"),side=1,adj=0.5,cex=1,line=1.4))
  dev.off()

}
#========== create Biome masks per continent ===================================

load(paste0(r_out, "continent_masks.RDATA"))
# to do: build loop over continents and write

## tropical
trop_forest_mask_america <- trop_forest_mask * america_mask
trop_forest_mask_africa  <- trop_forest_mask * africa_mask
trop_forest_mask_oceania <- trop_forest_mask * oceania_mask
trop_forest_mask_asia    <- trop_forest_mask * asia_mask

diff_trop <-(trop_forest_mask_oceania + trop_forest_mask_africa + trop_forest_mask_asia+ trop_forest_mask_america -trop_forest_mask)
if(sum(diff_trop)< (0 - length(trop_forest_mask[trop_forest_mask> 0])/100)) print("trop is fishy") #no more than 1% of the biome should occur outside of the defined masks

### temperate
temp_forest_mask_america <- temp_forest_mask * america_mask
temp_forest_mask_northam <- temp_forest_mask * northamerica_mask
temp_forest_mask_southam <- temp_forest_mask * southamerica_mask
temp_forest_mask_africa  <- temp_forest_mask * africa_mask
temp_forest_mask_oceania <- temp_forest_mask * oceania_mask
temp_forest_mask_asia    <- temp_forest_mask * asia_mask
temp_forest_mask_europe  <- temp_forest_mask * europe_mask

diff_temp <- temp_forest_mask_america + temp_forest_mask_europe +
 temp_forest_mask_africa + temp_forest_mask_asia + temp_forest_mask_oceania -temp_forest_mask
#if(sum(diff_temp)< (0 - length(temp_forest_mask[temp_forest_mask> 0])/100)) print("temp is fishy") #no more than 1% of the biome should occur outside of the defined masks


### boreal
boreal_forest_mask_america <- boreal_forest_mask * america_mask
boreal_forest_mask_asia    <-  boreal_forest_mask * asia_mask
boreal_forest_mask_europe  <-  boreal_forest_mask * europe_mask
boreal_forest_mask_eurasia <-  boreal_forest_mask_europe + boreal_forest_mask_asia

diff_bor <- boreal_forest_mask_america +boreal_forest_mask_asia + boreal_forest_mask_europe -boreal_forest_mask
#if(sum(diff_bor)< (0 - length(boreal_forest_mask[boreal_forest_mask> 0])/100)) print("boreal is fishy") #no more than 1% of the biome should occur outside of the defined masks

#========== save data ==========================================================

if (save == TRUE) {
  save(forest_mask,
     temp_forest_mask, trop_forest_mask, boreal_forest_mask,
     trop_forest_mask_oceania, trop_forest_mask_africa,
     trop_forest_mask_america, trop_forest_mask_asia,
     temp_forest_mask_america, temp_forest_mask_northam,
     temp_forest_mask_southam, temp_forest_mask_europe,
     temp_forest_mask_africa,  temp_forest_mask_asia, temp_forest_mask_oceania,
     boreal_forest_mask_america, boreal_forest_mask_asia,
     boreal_forest_mask_europe, boreal_forest_mask_eurasia,
     file = paste0(r_out, firstyear, "-", lastyear, "_biomes_pnv_lpjml.RDATA"))
}