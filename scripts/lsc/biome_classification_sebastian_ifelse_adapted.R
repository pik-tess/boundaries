# ifelse biome classification from Sebastian
# with the following modifications: 

# biome classification after Sebastian Ostberg (Ostberg et al. 2013 &
# Ostberg et al 2015) and Vera Heck for Gerten et al. 2020

# adaptations to account for two additional pfts: 
# polar grass & boreal deciduous needleleaved
# boreal decidious needleleaved: no new forest biome added, but boreal 
# decidious forest refers to both broadleaved and needleleaved deciduous 
# boreal forests
# polar grass: pft indices adapted accordingly, but no explicit consideration
# of polar grass within the conditions of the biome classification
# reason: polar grass very dominant over temperate grass
# --> artic tundra in temperate regions

# thresholds for woodland and open shrubland from Ostberg et al. 
# 2015 (>30 and 10 % tree cover respectively)

# additionally:
# additional condition for temperate grass to avoid temperate grass in very
# northern regions: latitude < 55

# gives (almost) the same results as
# biome_classification_sebastian_for_adapted


rm(list = ls(all = TRUE))

#========== libraries ==========================================================
library(lpjmliotools) # LPJmL package developed by Fabian
library(raster)
library(RColorBrewer)
library(maptools)
library(fields)

#========== directories ========================================================
dir <- c("/p/projects/open/Johanna/NETPS_within_PBs/BECCS_nat_veg/")
#dir_LPJmL_out <-  paste0(dir, "LPJmL_runs/PNV/output/")
dir_LPJmL_out = "/p/projects/open/Jannes/earth4all/pb_status/runs/output/pnv_newphen/"
dir_LPJmL_in <- paste0("/p/projects/lpjml/input/historical/")
r_out <- paste0(dir, "R_out/R_data/")
plotpath <- paste0(dir, "R_out/", "plots/biome_classification/")

#========== set variables ======================================================
ncell <- 67420
timeframe <- c(1979, 2014)
firstyear <- timeframe[1]
lastyear <- timeframe[2]
nyears <- lastyear - firstyear + 1
ref <- 1979 # first year of simulation
end <- 2014 # last year of simulation
fpc.unit  <- "[%]"
vegc.unit <- "[gC/m2]"
npft <- 12


#---- read in fpc, atemp and vegc ---------------------------------------------#

#----- fpc.bin -----------------------------------------------------------------
    fpc  <- readCFToutput(inFile = paste0(dir_LPJmL_out, "fpc.bin"),
                  startyear = ref, stopyear = end, bands = npft,
                  size = 4, headersize = 0, getyearstart = firstyear,
                  getyearstop = lastyear)
    fpc  <- apply(fpc, c(1, 2), mean, na.rm = TRUE) # mean over reference period

#----- vegc.bin ----------------------------------------------------------------
    vegc <- readYearly(inFile = paste0(dir_LPJmL_out, "vegc.bin"),
                  startyear = firstyear, stopyear = lastyear,
                  size = 4, headersize = 0, getyearstart = firstyear,
                  getyearstop = lastyear)
    vegc    <- apply(vegc, 1, mean, na.rm = TRUE) # mean over reference period

#----- atemp -------------------------------------------------------------------
    temp <- autoReadInput(inFile = paste0(dir_LPJmL_in, "CRU_TS4.03/",
                        "cru_ts4.03.1901.2018.tmp.clm"),
                        getyearstart = firstyear, getyearstop = lastyear,
                        manu = TRUE, msize = 2, mheadersize = 43)
    # has to be checked: read in correctly ?!
    # compare to Vera: days per months??
    temp <- apply(temp, c(1, 2), mean, na.rm = TRUE)
    #average value over all years
    atemp <- apply(temp, 2, mean, na.rm = TRUE)





#--- biome definitions --------------------------------------------------------#
{
  biome.class.names.lpj <- c("Tropical Rainforest", # 1
                             "Tropical Seasonal & Deciduous Forest", # 2
                             "Temperate Broadleaved Evergreen Forest", # 3
                             "Temperate Broadleaved Deciduous Forest", # 4
                             "Mixed Forest", # 6
                             "Temperate Coniferous Forest", # 5
                             "Boreal Evergreen Forest",	# 7
                             "Boreal Deciduous Forest", # 8
                             "Warm Woody Savanna, Woodland & Shrubland", # 9
                             "Warm Savanna & Open Shrubland", #10
                             "Warm Grassland", #11
                             "Temperate Woody Savanna, Woodland & Shrubland", # 12
                             "Temperate Savanna & Open Shrubland", #13
                             "Temperate Grassland",	#14
                             "Arctic Tundra", #15
                             "Desert", #16
                             "Rocks and Ice",
                             "Water"
  )

  fpc.names <- c("natvegfrac", 					#1
                 "Tropical Broadleaved Evergreen Tree",		#2
                 "Tropical Broadleaved Raingreen Tree",		#3
                 "Temperate Needleleaved Evergreen Tree",	#4
                 "Temperate Broadleaved Evergreen Tree",	#5
                 "Temperate Broadleaved Summergreen Tree",  #6
                 "Boreal Needleleaved Evergreen Tree",	    #7
                 "Boreal Broadleaved Summergreen Tree",     #8
                 "Boreal Needleleaved Summergreen Tree",    #9
                 "Tropical C4 grass",                       #10
                 "Temperate C3 grass",			    #11
                 "Polar C3 grass"		      		    #12
  )
  

  fpc.temperate <- c(4, 5, 6)
  fpc.tropical  <- c(2, 3)
  fpc.boreal    <- c(7, 8, 9)
  fpc.needle    <- c(4, 7, 9)
  fpc.grass     <- c(10, 11, 12)
  fpc.evergreen <- c(2, 4, 5,7)
  fpc.trees     <- sort(unique(c(fpc.temperate, fpc.tropical, fpc.boreal,
                        fpc.needle)))

}


#---- biome classification ----------------------------------------------------#
biome.classification <- function (in.fpc, in.atemp, in.vegc, use.vegc=TRUE) {
    if (dim(in.fpc)[1] != length(in.atemp)) {
      stop("Error: in.fpc and in.atemp do not match")
    }
    if (dim(in.fpc)[2] != length(fpc.names)) {
      stop("Error: in.fpc has wrong number of bands")
    }
    if (length(dim(in.fpc)) != 2) {
      stop("Error: in.fpc must be array of dimensions (pft, cell)")
    }
    fpc.tree.total     <- apply(in.fpc[, fpc.trees], 1, sum, na.rm = TRUE)
    fpc.tree.temperate <- apply(in.fpc[, fpc.temperate], 1, sum, na.rm = TRUE)
    fpc.tree.tropical  <- apply(in.fpc[, fpc.tropical], 1, sum, na.rm = TRUE)
    fpc.tree.boreal    <- apply(in.fpc[, fpc.boreal], 1, sum, na.rm = TRUE)
    fpc.tree.needle    <- apply(in.fpc[, fpc.needle], 1, sum, na.rm = TRUE)
    fpc.tree.evergreen <- apply(in.fpc[, fpc.evergreen], 1, sum, na.rm = TRUE)
    fpc.grass.total    <- apply(in.fpc[, fpc.grass], 1, sum, na.rm = TRUE)
    fpc.total          <- apply(in.fpc[, -1], 1, sum, na.rm = TRUE)
    max.share          <- apply(in.fpc[, -1], 1, max, na.rm = TRUE)
    # -1 cuts the first entry off
    fpc.tree.broadleaf <- fpc.tree.total - fpc.tree.needle

    biome.class <- integer(length(in.atemp))
    biome.class[] <- NA

    if (!use.vegc) {
      in.vegc[] <- 7500 # vegc test always positive
    }

biome.class <- ifelse(fpc.total <= 0.05,
                          ifelse(in.atemp < (-2), which(biome.class.names.lpj == "Arctic Tundra"), which(biome.class.names.lpj == "Desert")),
                          ifelse(fpc.tree.total >= 0.6,
                                 ifelse(in.fpc[,7]==max.share & (fpc.tree.broadleaf < (0.4*fpc.tree.total)),
                                        which(biome.class.names.lpj=="Boreal Evergreen Forest"),
                                        ifelse((fpc.tree.evergreen < (0.4*fpc.tree.total) & (in.fpc[,8]==max.share | in.fpc[,9]== max.share)),
                                               which(biome.class.names.lpj=="Boreal Deciduous Forest"),
                                               ifelse((in.fpc[,4]==max.share) & (fpc.tree.broadleaf < (0.4*fpc.tree.total)),# (in.fpc[,4]> (in.fpc[,5]+in.fpc[,6])),
                                                      which(biome.class.names.lpj=="Temperate Coniferous Forest"),
                                                      ifelse((fpc.tree.tropical<0.4*fpc.tree.total) & (fpc.tree.needle < (0.4*fpc.tree.total)), #in.fpc[,5]+in.fpc[,6]>=0.6,
                                                             ifelse(in.fpc[,5]>in.fpc[,6],
                                                                    which(biome.class.names.lpj=="Temperate Broadleaved Evergreen Forest"),
                                                                    which(biome.class.names.lpj=="Temperate Broadleaved Deciduous Forest")),
                                                             ifelse((fpc.tree.boreal+fpc.tree.temperate)<(0.4*fpc.tree.total),
                                                                    ifelse(in.fpc[,2]>in.fpc[,3] & in.vegc >= 7500,
                                                                           which(biome.class.names.lpj=="Tropical Rainforest"),
                                                                           ifelse(in.vegc >= 7500,
                                                                                  which(biome.class.names.lpj=="Tropical Seasonal & Deciduous Forest"),
                                                                                  which(biome.class.names.lpj=="Warm Woody Savanna, Woodland & Shrubland")
                                                                           )
                                                                    ), # end tropical forests/savanna
                                                                    which(biome.class.names.lpj=="Mixed Forest")
                                                             ) # end tropical or mixed
                                                      )# end Temperate Broadleaved
                                               )# end Temperate
                                        )# end Boreal Deciduous
                                 ), # end Boreal Evergreen and end trees
                                 ifelse(fpc.tree.total>=(0.3) & in.atemp >= (-2) & lat < 55,# & fpc.tree.boreal==0,
                                        ifelse(in.fpc[,11] > in.fpc[,10],
                                               which(biome.class.names.lpj=="Temperate Woody Savanna, Woodland & Shrubland"), #12
                                               which(biome.class.names.lpj=="Warm Woody Savanna, Woodland & Shrubland") #9
                                        ),
                                        ifelse(fpc.tree.total>=(0.1) & in.atemp >= (-2) & lat < 55,# & fpc.tree.boreal==0,
                                               ifelse(in.fpc[,11]>in.fpc[,10],
                                                      which(biome.class.names.lpj=="Temperate Savanna & Open Shrubland"), #13
                                                      which(biome.class.names.lpj=="Warm Savanna & Open Shrubland") #10
                                               ),
                                               ifelse(in.atemp >= (-2) & lat < 55,
                                                      ifelse(in.fpc[,11]>in.fpc[,10],
                                                             which(biome.class.names.lpj=="Temperate Grassland"), #14
                                                             which(biome.class.names.lpj=="Warm Grassland") #11
                                                      ),
                                                      which(biome.class.names.lpj=="Arctic Tundra") #15
                                               ) # end grasslands
                                        ) # end dry savanna
                                 ) # end moist savanna
                          ) # end forest or savanna
    ) # end all ifelse
      
return(biome.class)
}
print("Computing biome classes for reference period")
biome_classes <- biome.classification(fpc, atemp, vegc)


#---- plotting ----------------------------------------------------------------#

 biome.class.cols <- c(rev(brewer.pal(6, "YlOrBr")), # warm
    rev(brewer.pal(9, "YlGn")[c(3,5,7,9)]),
    rev(brewer.pal(9, "GnBu"))[c(2:4,6,8,9)], # cold below forest
    "grey", #Rocks & Ice
    "lightblue" # Water
  )

 biome.class.cols <-  biome.class.cols[c(1,2,7,8,9,10,13,12,3,4,5,14,15,16,11,6,17,18)]
 #biome.class.cols <- c(rev(brewer.pal(9, "YlGn"))[-9], # all trees
 #   rev(brewer.pal(9, "Oranges"))[4:6], # warm below forest
 #   rev(brewer.pal(9, "GnBu"))[2:5], # cold below forest
 #   "cornsilk4", # desert
 #   "grey", #Rocks & Ice
 #   "lightblue" # Water
 # )


#test

png(file= paste0(plotpath, "ostberg_adapted_ifelse_5.3fix.png"), width = 6.8, height = 3.6, units = "in", res = 300, 
            pointsize = 6, type = "cairo")
par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
par(cex = 1.2, cex.main = 1.2)
lpj.ras         <- raster(ncols=720,nrows=360)
lpj.ras[cellFromXY(lpj.ras,cbind(lon,lat))] <- biome_classes
plot(lpj.ras, breaks = seq(0.5, 18.5, 1), col = biome.class.cols,  legend = FALSE)
dev.off()


legend("bottom", legend=biome.class.names.lpj,
       fill=biome.class.cols, cex=0.8)

test<- numeric(length=ncell)
test[which(atemp <= 0)] <- 2

lpj.ras         <- raster(ncols=720,nrows=360)
lpj.ras[cellFromXY(lpj.ras,cbind(lon,lat))] <- test
plot(lpj.ras)

# ice is missing in classification