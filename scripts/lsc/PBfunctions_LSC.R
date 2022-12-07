####### BIOME Classification #################################################
## necessary for the assessment of the land-system change planetary boundary
## Arne Tobian
## 211207

# list of PFTs (LPJ5.1)
fpc.names <- c("natvegfrac", 														#1 
               "Tropical Broadleaved Evergreen tree",		#2;  ID: TROPICAL_BROADLEAVED_EVERGREEN_TREE
               "Tropical Broadleaved Raingreen tree",		#3;  ID: TROPICAL_BROADLEAVED_RAINGREEN_TREE
               "Temperate Needleleaved Evergreen tree",	#4;  ID: TEMPERATE_NEEDLELEAVED_EVERGREEN_TREE
               "Temperate Broadleaved Evergreen tree",	#5;  ID: TEMPERATE_BROADLEAVED_EVERGREEN_TREE
               "Temperate Broadleaved Eummergreen tree",#6;  ID: TEMPERATE_BROADLEAVED_SUMMERGREEN_TREE
               "Boreal Needleleaved Evergreen tree",		#7;  ID: BOREAL_NEEDLELEAVED_EVERGREEN_TREE
               "Boreal Broadleaved Summergreen tree",   #8;  ID: BOREAL_BROADLEAVED_SUMMERGREEN_TREE 
               "Boreal Needleleaved Summergreen Tree",	#9;  ID: BOREAL_NEEDLELEAVED_SUMMERGREEN_TREE -> wie broadleaved
               "Tropical C4 grass",										  #10; ID: TROPICAL_HERBACEOUS
               "Temperate C3 grass",										#11; ID: TEMPERATE_HERBACEOUS
               "Polar C3 grass"											    #12; ID: POLAR_HERBACEOUS -> wie temperate C3 oder artic tundra?
)

#define input and output directories
indir      <- ("/../../") 
outdir      <- ("/../../") 

#define timespan (for averaging biome masks over a period)
timespan <- c(startyear,stopyear) #corresponds to 1850-1859

####### Calc FPC covers ##########################################
#stack fpc bands and calculate mean
for (i in 1:12) {
  nam <- paste0("fpc_band_",i)
  fpc_temp <- brick(paste0(indir,"fpc.nc"), level = i) #obs! other routine for bindary data
  fpc_temp <- fpc_temp[[timespan[1]:timespan[2]]]    
  fpc_temp <- calc(fpc_temp,mean)
  assign(nam, fpc_temp)
}

fpc.r.stack <- stack(fpc_band_1,fpc_band_2,fpc_band_3,
                     fpc_band_4,fpc_band_5,fpc_band_6,
                     fpc_band_7,fpc_band_8,fpc_band_9,
                     fpc_band_10,fpc_band_11,fpc_band_12)

fpc.tree.temp      <- (fpc.r.stack[[4]] + fpc.r.stack[[5]] +
                         fpc.r.stack[[6]])*fpc.r.stack[[1]] 
fpc.tree.trop      <- (fpc.r.stack[[2]] + fpc.r.stack[[3]])*fpc.r.stack[[1]]
fpc.tree.boreal    <- (fpc.r.stack[[7]] + fpc.r.stack[[8]]  + fpc.r.stack[[9]])
                        *fpc.r.stack[[1]]  #c(7,8,9) 
fpc.tree.total     <- fpc.tree.temp + fpc.tree.trop + fpc.tree.boreal
fpc.total          <- (sum(fpc.r.stack) - fpc.r.stack[[1]])*fpc.r.stack[[1]]

####### Forest biome classification ################################

# more than 50% of total tree cover are by biome-PFTs 
# more than 5% of the cell are covered by vegetation
# tree coverage is more than 60% per cell -> FAo ab 10%

# simplified biome classification after Sebastian Ostberg (2013) and Vera Heck
fpc.tree.temp [fpc.tree.temp < (0.5*fpc.tree.total) | 
                 fpc.tree.total <= 0.6 |fpc.total <= 0.05 ] <- 0
fpc.tree.trop [fpc.tree.trop < (0.5*fpc.tree.total) | 
                 fpc.tree.total <= 0.6 |fpc.total <= 0.05 ] <- 0
fpc.tree.boreal[fpc.tree.boreal < (0.5*fpc.tree.total) | 
                  fpc.tree.total <= 0.6 | fpc.total <= 0.05] <- 0 

####### BIOME masks #################################################

#LPJ forest biomes mask

empty_mask <- fpc.r.stack[[1]]
empty_mask <- empty_mask*0

forest_mask <- empty_mask
forest_mask[fpc.tree.trop > 0 | fpc.tree.temp > 0 | fpc.tree.boreal > 0] <-1

#trop forets mask 
trop_forest_mask_pnv <- empty_mask
trop_forest_mask_pnv[fpc.tree.trop > 0] <-1

#temp
temp_forest_mask_pnv <- empty_mask
temp_forest_mask_pnv[fpc.tree.temp > 0] <-1

#boreal
boreal_forest_mask_pnv <- empty_mask
boreal_forest_mask_pnv[fpc.tree.boreal > 0] <-1

africa <-continents[continents$CONTINENT=="Africa",]
southamerica <-continents[continents$CONTINENT=="South America",]
northamerica <-continents[continents$CONTINENT=="North America",]
europe <-continents[continents$CONTINENT=="Europe",]
asia <-continents[continents$CONTINENT=="Asia",]
oceania <- continents[continents$CONTINENT=="Oceania",]
australia <- continents[continents$CONTINENT=="Australia",]

america <- bind(southamerica,northamerica)
oceania <- bind(oceania,australia)

biome.list <- c("trop","temp","boreal")
cont.list <- c("africa","america","europe","oceania","asia")

for (i in 1:length(biome.list)){
  forestmask.name <- paste0(biome.list[i],"_forest_mask_pnv")
  forestmask <- get(forestmask.name)
  for (j in 1: length(cont.list)){
    cont <- get(cont.list[j])
    # mask.raster <- crop(mask(forestmask, cont),cont) #adjust LU pattern
    mask.raster <- mask(forestmask, cont) #adjust LU pattern
    if(maxValue(mask.raster) == 0) next
    mask.name <- paste0(forestmask.name,"_",cont.list[j])
    assign(mask.name,mask.raster)
    #plot(mask.raster, main = cont.list[j])
  }
}

####### Export masks as Rdata ###########################################


# saves 
save(list = ls(.GlobalEnv, pattern = "trop_forest_mask|temp_forest_mask|boreal_forest_mask"), file = paste0(outdir,"biomes.Rdata"))
