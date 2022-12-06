library(lpjmliotools)
library(raster)
library(maptools)

dir_input_data   <- "/p/projects/open/Vera/synthesepaper/input_data/"
dir <- c("/p/projects/open/Johanna/NETPS_within_PBs/")
dir_r_out <- paste0(dir, "R/r_out/r_data/")

#========== load data ==========================================================

continents <- readShapeSpatial(paste(dir_input_data,
                   "WorldContinents/world_continents.shp",sep=""),verbose=FALSE)
  
africa <- continents[continents$CONTINENT == "Africa", ]
southamerica <- continents[continents$CONTINENT == "South America", ]
northamerica <- continents[continents$CONTINENT == "North America", ]
europe <- continents[continents$CONTINENT == "Europe", ]
asia <- continents[continents$CONTINENT == "Asia", ]
oceania <- continents[continents$CONTINENT == "Oceania", ]
australia <- continents[continents$CONTINENT == "Australia", ]
  
id_v           <- 1:ncell
id.ras         <- raster(ncols = 720, nrows = 360)
id.ras[cellFromXY(id.ras,cbind(lon, lat))] <- id_v
  
extr_afr     <- extract(id.ras, africa, small = T, weights = T)
extr_southam <- extract(id.ras, southamerica, small = T, weights = T)
extr_northam <- extract(id.ras, northamerica, small = T, weights = T)
extr_europe  <- extract(id.ras, europe, small = T, weights = T)
extr_asia    <- extract(id.ras, asia, small = T, weights = T)
extr_oceania    <- extract(id.ras, oceania, small = T, weights = T)
extr_australia    <- extract(id.ras, australia, small = T, weights = T)
  
  
africa_mask <-  southam_mask <-  northam_mask <- europe_mask <- asia_mask <-
                            oceania_mask <- australia_mask <- array(0, ncell)
for (i in 1:ncell) {
    if (any(extr_afr[[1]][, 1] == i, na.rm=T) == T) africa_mask[i] <- 1
    else if (any(extr_southam[[1]][,1]==i,na.rm=T)==T)  southam_mask[i] <- 1
    else if (any(extr_northam[[1]][,1]==i,na.rm=T)==T)  northam_mask[i] <- 1
    else if (any(extr_europe[[1]][,1]==i,na.rm=T)==T)   europe_mask[i] <- 1
    else if (any(extr_asia[[1]][,1]==i,na.rm=T)==T)     asia_mask[i] <- 1
    else if (any(extr_oceania[[1]][,1] == i, na.rm = T) == T)  oceania_mask[i] <- 1
    else if (any(extr_australia[[1]][,1] == i, na.rm = T) == T)  australia_mask[i] <- 1
}
  oceania_mask <- australia_mask + oceania_mask
  america_mask <- northam_mask + southam_mask
  
  e <- europe_mask
  
  europe_mask <- e
  #### crop continent biomes #####
  crop_vec <- function(vec, extent) {
    # extent must be x1,x2,y1,y2
    vec[lon < extent[1]] <- 0
    vec[lon > extent[2]] <- 0
    vec[lat < extent[3]] <- 0
    vec[lat > extent[4]] <- 0
    return(vec)
  }
  
 test <- rep(0,67420)
 #extent <- c(-121,-25,-90,5)
 extent <- c(-151,-35,12,80)
 test[lon < extent[1] | lon > extent[2] | lat < extent[3] | lat > extent[4]] <- 2
 lpj.ras         <- raster(ncols=720,nrows=360)
 lpj.ras[cellFromXY(lpj.ras,cbind(lon,lat))] <- test
 plot(lpj.ras)


  ## tropical
  america2 <- crop_vec(array(1,ncell),c(-121,-25,-90,30))
  america3 <- crop_vec(array(1,ncell),c(-151,-35,0,80))
  america_south <- crop_vec(array(1,ncell),c(-121,-25,-90,5))
  america_north <- crop_vec(array(1,ncell),c(-151,-35,12,80))
  africa2  <- crop_vec(array(1,ncell),c(-30, 85, -47, 5))
  oceania2 <- crop_vec(array(1,ncell),c(-190, -121, -39, 31))
  oceania3 <- crop_vec(array(1,ncell),c(150, 180, -39, 31))
  oceania4 <- crop_vec(array(1,ncell),c(135, 180, -10, 21))
  oceania5 <- crop_vec(array(1,ncell),c(100, 180, -50, -12))
  asia2  <- crop_vec(array(1,ncell),c(53, 150, -10, 90))
  europe2  <- crop_vec(array(1,ncell),c(-20, 40, 38, 65))
  europe3  <- crop_vec(array(1,ncell),c(-35, -15, 30, 50))
  
  oceania_mask <- oceania_mask +oceania2 +oceania3 + oceania4 + oceania5
  oceania_mask[oceania_mask>1]=1
  america_mask <- america_mask + america2+ america3
  america_mask[america_mask>1]=1
  southamerica_mask <- southam_mask + america_south
  southamerica_mask[southamerica_mask > 1] <- 1
  northamerica_mask <- northam_mask + america_north
  northamerica_mask[northamerica_mask > 1] <- 1
  asia_mask <- asia_mask + asia2
  asia_mask[asia_mask>1]=1
  africa_mask <- africa_mask + africa2
  africa_mask[africa_mask>1]=1
  europe_mask <- europe_mask + europe2 +europe3
  europe_mask[europe_mask>1]=1
  
  length(which(oceania_mask+america_mask>1))
  oceania_mask[oceania_mask+america_mask>1]=0
  
  length(which(oceania_mask+asia_mask>1))
  oceania_mask[oceania_mask+asia_mask>1]=0
  
  length(which(europe_mask+asia_mask>1))
  europe_mask[europe_mask+asia_mask>1]=0
  
  length(which(europe_mask+africa_mask>1))
  europe_mask[europe_mask+africa_mask>1]=0
  
  length(which(asia_mask+africa_mask>1))
  asia_mask[asia_mask+africa_mask>1]=0
  
  length(which(america_mask+africa_mask>1))
  america_mask[america_mask+africa_mask>1]=0
  
  
  plot.ras         <- raster(ncols=720,nrows=360)
  plot.ras[cellFromXY(plot.ras,cbind(lon,lat))] <-  northam_mask
  plot(plot.ras)
  save(africa_mask,america_mask,europe_mask, asia_mask, oceania_mask,
  southamerica_mask, northamerica_mask,
       file=paste(dir_r_out,"continent_masks.RDATA",sep=""))
