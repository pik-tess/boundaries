# written by Fabian Stenzel
# 2022 - stenzel@pik-potsdam.de

## todo: documentation, HANPP find scaling factor for grassland harvest

netcdfPFT2lpjarray <- function(ncInFile,var,lon,lat){
  #ncInfile = "/p/projects/amazonas/drueke/POEM/work/fabian_LU/history/output-20040101/pft_npp.nc"
  #var="NPP"
  require(raster)
  require(ncdf4)
  # NetCDF file
  file_nc <- nc_open(ncInFile)
  # get spatial extent and resolution
  # this will give a warning if the NetCDF has more than one data field, e.g. crop bands or time axis
  file_raster <- raster(ncInFile)
  cells <- cellFromXY(file_raster, cbind(lon,lat))
  time <- ncvar_get(file_nc,"time")
  years <- length(time)
  data <- ncvar_get(file_nc, var)
  #flip latitude according to netcdf
  if ((yFromRow(file_raster, 1) < yFromRow(file_raster,  2))
      != (file_nc$dim$lat$vals[1] < file_nc$dim$lat$vals[2])) {
    if (years == 1) data <- data[, file_nc$dim$lat$len:1, ]
    else data <- data[, file_nc$dim$lat$len:1,, ]
  }
  nbands <- dim(data)[3]
  nat <- apply(data[,,1:11,],c(1,2,4),sum)
  hum <- apply(data[,,12:43,],c(1,2,4),sum)
  rm(data)
  #print(paste0(length(lon),nbands))
  if (years == 1)   outdata <- array(0,dim = c(length(lon)))
  else outdata <- array(0,dim = c(length(lon),years))
  # loop over individual days and years in NetCDF
  for (year in 1:years) {
    for (band in 1:nbands) {
      if (years == 1) dat <- data[,,band]
      else dat <- data[,,band,year]
      # dissolve lon and lat dimensions into one
      dim(dat) <- c(file_nc$dim$lat$len*file_nc$dim$lon$len)
      # extract data for LPJmL coordinates using cellFromXY()
      # resulting array has dimensions c(ncell, nbands)
      if (years == 1) outdata[,band] <- dat[cells]
      else outdata[,band,year] <- dat[cells]

    }
  }
  nc_close(file_nc)
  return(outdata)
}
plotNPPmap <- function(data, file, title, legendtitle, eps = FALSE){
  brks <- seq(0,1000,100)
  palette <- c("orangered4","orangered","orange","gold","greenyellow","limegreen","green4","darkcyan","darkslategrey","navy")
  data[data<brks[1]] <- brks[1]
  data[data>brks[length(brks)]] <- brks[length(brks)]
  if (eps) {
    file = strsplit(file,".",fixed=TRUE)[[1]]
    file = paste(c(file[1:(length(file) - 1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width = 22,  height = 8.5,paper="special")
  }else{
    png(file, width=7.25,  height = 3.5, units = "in", res = 300, pointsize = 6,type="cairo")
  }
  ra <- raster::raster(ncols=720, nrows=360)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  extent <- raster::extent(c(-180, 180, -60, 90))
  par(bty="n",oma=c(0,0,0,0),mar=c(0,0,0,3),xpd=T)
  raster::plot(ra,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes = FALSE)
  title(title,line=-1)
  fields::image.plot(legend.only=TRUE,zlim=range(brks),col = palette,useRaster=FALSE,breaks=brks,lab.breaks=brks,legend.shrink = 0.8,
                     legend.args=list(legendtitle,side=3, font=2, line=1))
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))
  dev.off()
}
plotHANPPmap <- function(data, file, title, legendtitle, zeroThreshold = 0.1, eps = FALSE){
  brks <- c(-400,-200,-100,-50,-zeroThreshold,zeroThreshold,10,20,30,40,50,60,70,80,100)
  classes <- c("<-200","-200 - -100","-100 - -50",paste0("-50 - -",zeroThreshold),paste0("-",zeroThreshold," - ",zeroThreshold),paste0(zeroThreshold," - 10"),"10 - 20","20 - 30","30 - 40","40 - 50","50 - 60","60 - 70","70 - 80","80 - 100")
  palette <- c("navy","royalblue3","royalblue1","skyblue1","grey80","yellowgreen","greenyellow","yellow","gold","orange","orangered","orangered4","brown4","black")
  data[data < brks[1]] <- brks[1]
  data[data > brks[length(brks)]] <- brks[length(brks)]
  if (eps) {
    file = strsplit(file,".",fixed = TRUE)[[1]]
    file = paste(c(file[1:(length(file) - 1)],"eps"), collapse = ".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width = 22,  height = 8.5, paper = "special")
  }else{
    png(file, width = 7.25, height = 3.5, units = "in", res = 300, pointsize = 6, type = "cairo")
  }
  ra <- raster::raster(ncols = 720, nrows = 360)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  extent <- raster::extent(c(-180, 180, -60, 90))
  par(bty = "n", oma = c(0,0,0,0), mar = c(0,0,0,0), xpd = T)
  raster::plot(ra, ext = extent, breaks = brks, col = palette, main = "", legend = FALSE, axes = FALSE)
  title(title, line = -2)
  maps::map('world', add = TRUE, res = 0.4, lwd = 0.25,ylim = c(-60,90))
  legend(x = -180, y = 50, fill = palette, border = palette, legend = classes, title = legendtitle)
  dev.off()
}
plotNPPovertime <- function(file, dataAct, dataEco, dataPot, dataBE, dataHANPP,
                            dataHarvNPP, dataHANPPperc, firstyr, plotyrs, highlightyrs = 2000, minVal = 0,
                            maxVal = 100, legendpos="topleft", ext = FALSE, eps = FALSE){
  #maxVal=max(c(dataPot,dataAct,dataEco,dataHANPP))
  lastyr = firstyr + length(dataAct) - 1
  colz = c("slateblue","gold","green3","red3","darkorange","black")
  if (eps) {
    file = strsplit(file,".",fixed=TRUE)[[1]]
    file = paste(c(file[1:(length(file) - 1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width = 22,  height = 8.5,paper="special")
  }else{
    png(file, width=3.5,  height = 3, units = "in", res = 300, pointsize = 6,type="cairo")
  }
  par(bty="o",oma=c(0,0,0,0),mar=c(4,5,1,3))
  plot(x=seq(firstyr,lastyr,1),y=dataPot,ylab="GtC/yr",xlab="Year",xlim=plotyrs,
       ylim=c(minVal, maxVal),type = "l",col=colz[1],xaxs="i",yaxs="i")
  lines(x=seq(firstyr,lastyr,1),y=dataAct,type = "l",col=colz[2])
  lines(x=seq(firstyr,lastyr,1),y=dataEco,type = "l",col=colz[3])

  par(bty="n",oma=c(0,0,0,0),mar=c(4,5,1,3), new = T)
  plot(x=seq(firstyr,lastyr,1),y=dataHANPPperc,ylab="",xlab="",xlim=plotyrs,
       ylim=c(10, 30),type = "l",col=colz[4],xaxs="i", yaxs="i", axes = F)
  axis(side = 4, col = colz[4],col.axis = colz[4])
  mtext(text = "%", col=colz[4], side = 4,line = 2)

  if (ext){
    lines(x=seq(firstyr,lastyr,1),y=dataBE,type = "l",col=colz[5])
    HANPPperc=dataHANPP/dataPot*100
    lines(x=seq(firstyr,lastyr,1),y=HANPPperc,type = "l",col="black")
  }
  if (!is.null(highlightyrs)){
    for (y in highlightyrs){
      lines(x=c(y,y),y=c(minVal,maxVal),col="grey40")
    }
  }

  if (ext) legend(legendpos,legend = c("NPPpot (PNV)","NPPact (landuse)","NPPeco","HANPP","Bioenergy","total NPP on cropland"),col=colz ,lty=1,cex=1)
  else legend(legendpos,legend = c("NPPpot (PNV)","NPPact (landuse)","NPPeco","M-COL [% NPPpi]"),col=colz[1:4] ,lty=1,cex=1)
  dev.off()
}
readHANPPdata <- function(inFol_lu, inFol_pnv, startyr, stopyr, gridbased = T,
                          pftbands = 11, cftbands = 32, headersize = 0, dataFile,
                          ncells = 67420, fileType = "clm", varnames, saveDataFile = F) {
  totalbands = (pftbands + cftbands)
  if (fileType == "clm") {
    nppFile = paste0(inFol_lu, varnames["npp","outname"])
    print(paste0("reading NPP file: ",nppFile))
    if (varnames["npp","timestep"] == "Y") {
      ynpp = readYearly(inFile = nppFile,startyear = startyr,stopyear = stopyr,size = 4,
                        headersize = headersize, getyearstart = startyr,getyearstop = stopyr)# gC/m2
    }else if (varnames["npp","timestep"] == "M") {
      npp = readMonthly(inFile = nppFile,startyear = startyr,stopyear = stopyr,size = 4,
                        headersize = headersize, getyearstart = startyr,getyearstop = stopyr)# gC/m2
      ynpp = apply(npp, c(1,3), sum) #gC/m2
    }
    pft_nppFile = paste0(inFol_lu,varnames["pft_npp","outname"])
    print(paste0("reading PFT_NPP file: ",pft_nppFile))
    pftnpp = readCFToutput(inFile = pft_nppFile,startyear = startyr,stopyear = stopyr,size = 4, headersize = headersize, bands=totalbands,getyearstart = startyr,getyearstop = stopyr)
    harvfile = paste0(inFol_lu,varnames["pft_harvest","outname"])
    rharvfile = paste0(inFol_lu,varnames["pft_rharvest","outname"])
    print(paste0("reading harvest file: ",harvfile))
    harvest = readCFToutput(inFile = harvfile,startyear = startyr,stopyear = stopyr,size = 4, headersize = headersize, bands=cftbands,getyearstart = startyr,getyearstop = stopyr)
    print(paste0("reading rharvest file: ",rharvfile))
    rharvest = readCFToutput(inFile = rharvfile,startyear = startyr,stopyear = stopyr,size = 4, headersize = headersize, bands=cftbands,getyearstart = startyr,getyearstop = stopyr)
    timberFile = paste0(inFol_lu,varnames["timber_harvest","outname"])
    print(paste0("reading timber file: ",timberFile))
    timber = readYearly(inFile = timberFile,startyear = startyr,stopyear = stopyr,size = 4, headersize = headersize, getyearstart = startyr,getyearstop = stopyr,ncells = ncells)
    fireFile = paste0(inFol_lu,varnames["firec","outname"])
    print(paste0("reading fire file: ",fireFile))
    fire = readYearly(inFile = fireFile,startyear = startyr,stopyear = stopyr,size = 4, headersize = headersize, getyearstart = startyr,getyearstop = stopyr,ncells = ncells)
    cftFile = paste0(inFol_lu,varnames["cftfrac","outname"])
    print(paste0("reading cftfrac file: ",cftFile))
    cftfrac = readCFToutput(inFile = cftFile,startyear = startyr,stopyear = stopyr,size = 4,bands = cftbands, headersize = headersize, getyearstart = startyr,getyearstop = stopyr)
    pnv_nppFile = paste0(inFol_pnv,varnames["npp","outname"])
    print(paste0("reading PNV NPP file: ",pnv_nppFile))
    if (varnames["npp","timestep"] == "Y") {
      ynpp_potential = readYearly(inFile = pnv_nppFile,startyear = startyr,stopyear = stopyr,size = 4, headersize = headersize, getyearstart = startyr,getyearstop = stopyr)# gC/m2
    }else if (varnames["npp","timestep"] == "M") {
      npp_potential = readMonthly(inFile = pnv_nppFile,startyear = startyr,stopyear = stopyr,size = 4, headersize = headersize, getyearstart = startyr,getyearstop = stopyr)# gC/m2
      ynpp_potential = apply(npp_potential, c(1,3), sum) #gC/m2
    }
    fpcFile = paste0(inFol_lu,varnames["fpc","outname"])
    print(paste0("reading FPC file: ",fpcFile))
    fpc = readCFToutput(inFile = fpcFile,startyear = startyr,stopyear = stopyr,size = 4,bands = (pftbands + 1), headersize = headersize, getyearstart = startyr,getyearstop = stopyr)
  }else if (fileType == "nc") {
    nppFile = paste0(inFol_lu,varnames["pft_npp","outname"])
    print(paste0("reading NPP file: ",nppFile))
    npp = netcdfMonthly2lpjarray(ncInFile = nppFile,var = "NPP",timeVar="time",lon = lon,lat=lat)# gC/m2
    print(paste0("reading PFT_NPP file: ",pft_nppFile))
    pftnpp = netcdfCFT2lpjarray(ncInFile = pft_nppFile,var = "NPP",lon = lon,lat=lat)
    print(paste0("reading cftfrac file: ",cftFile))
    cftfrac = netcdfCFT2lpjarray(ncInFile = cftFile,var = "cftfrac",lon = lon,lat=lat)
    print(paste0("reading PNV NPP file: ",pnv_nppFile))
    npp_potential = netcdfMonthly2lpjarray(ncInFile = pnv_nppFile,var = "NPP",timeVar="time",lon = lon,lat=lat)# gC/m2
    print(paste0("reading FPC file: ",fpcFile))
    fpc = netcdfCFT2lpjarray(ncInFile = fpcFile,var = "fpc",lon = lon, lat = lat)

    ynpp = apply(npp, c(1,3), sum) #gC/m2
    ynpp_potential = apply(npp_potential, c(1,3), sum) #gC/m2
  }else if (fileType == "meta") {
      nppFile = paste0(inFol_lu, varnames["npp","outname"])
      if (varnames["npp","timestep"] == "Y") {
        ynpp = autoReadMetaOutput (metaFile = nppFile, getyearstart = startyr, getyearstop = stopyr,suppressReadPrint = T,headersize = headersize)[,1,]
      }else if (varnames["npp","timestep"] == "M") {
        npp = autoReadMetaOutput (metaFile = nppFile, getyearstart = startyr, getyearstop = stopyr,suppressReadPrint = T,headersize = headersize)[,1,,]
        ynpp = apply(npp, c(1,3), sum) #gC/m2
      }
      pft_nppFile = paste0(inFol_lu,varnames["pft_npp","outname"])
      pftnpp = autoReadMetaOutput(metaFile = pft_nppFile,getyearstart = startyr, getyearstop = stopyr,suppressReadPrint = T,headersize = headersize)
      harvfile = paste0(inFol_lu,varnames["pft_harvest","outname"])
      rharvfile = paste0(inFol_lu,varnames["pft_rharvest","outname"])
      harvest = autoReadMetaOutput(metaFile = harvfile,getyearstart = startyr, getyearstop = stopyr,suppressReadPrint = T,headersize = headersize)
      rharvest = autoReadMetaOutput(metaFile = rharvfile,getyearstart = startyr, getyearstop = stopyr,suppressReadPrint = T,headersize = headersize)
      timberFile = paste0(inFol_lu,varnames["timber_harvest","outname"])
      timber = autoReadMetaOutput(metaFile = timberFile,getyearstart = startyr, getyearstop = stopyr,suppressReadPrint = T,headersize = headersize)
      fireFile = paste0(inFol_lu,varnames["firec","outname"])
      fire = autoReadMetaOutput(metaFile = timberFile,getyearstart = startyr, getyearstop = stopyr,suppressReadPrint = T,headersize = headersize)
      cftFile = paste0(inFol_lu,varnames["cftfrac","outname"])
      cftfrac = autoReadMetaOutput(metaFile = timberFile,getyearstart = startyr, getyearstop = stopyr,suppressReadPrint = T,headersize = headersize)
      pnv_nppFile = paste0(inFol_pnv,varnames["npp","outname"])
      if (varnames["npp","timestep"] == "Y") {
        ynpp_potential = autoReadMetaOutput(metaFile = pnv_nppFile,getyearstart = startyr, getyearstop = stopyr,suppressReadPrint = T,headersize = headersize)[,1,]
      }else if (varnames["npp","timestep"] == "M") {
        npp_potential = autoReadMetaOutput(metaFile = pnv_nppFile,getyearstart = startyr, getyearstop = stopyr,suppressReadPrint = T,headersize = headersize)[,1,,]
        ynpp_potential = apply(npp_potential, c(1,3), sum) #gC/m2
      }
      fpcFile = paste0(inFol_lu,varnames["fpc","outname"])
      fpc = autoReadMetaOutput(metaFile = fpcFile,getyearstart = startyr, getyearstop = stopyr,suppressReadPrint = T,headersize = headersize)

  }else{
    print(paste0("Unrecognized file type (",fileType,")"))
    break
  }
  if (gridbased) {
    pftnpp_grasslands = apply(pftnpp[, c(totalbands - 18, totalbands - 2), ],c(1,3),sum) #gC/m2
    pftnpp_cft = apply(pftnpp[,-c(1:(pftbands), totalbands - 18, totalbands - 2), ], c(1,3), sum) #gC/m2
    pftnpp_bioenergy = apply(pftnpp[, c(totalbands - 17, totalbands - 16, totalbands - 1, totalbands), ], c(1,3), sum) #gC/m2
    #harvest_grasslands=apply(harvest[,c(totalbands-18,totalbands-2),],c(1,3),sum) #gC/m2
    #harvest_bioenergy=apply(harvest[,c(totalbands-17,totalbands-16,totalbands-1,totalbands),],c(1,3),sum) #gC/m2
    harvest_cft = apply(harvest[,,], c(1,3), sum) #gC/m2
    rharvest_cft = apply(rharvest[,,], c(1,3), sum) #gC/m2
    pftnpp_nat = apply(pftnpp[,1:(pftbands),], c(1,3), sum) #gC/m2
  }else{# todo: complete this part
    pftnpp_cft = apply(pftnpp[,(totalbands - cftbands + 1):totalbands,]*cftfrac,c(1,3),sum) #gC/m2
    pftnpp_nat = apply(pftnpp[,1:(pftbands),],c(1,3),sum)*fpc[,1,] #gC/m2
  }
  if (saveDataFile) {
    if (!file.exists(dataFile) ) {
      print(paste0("Writing data file: ",dataFile))
      save(ynpp_potential,ynpp,pftnpp_cft,pftnpp_nat,pftnpp_grasslands,pftnpp_bioenergy,harvest_cft,rharvest_cft,fire,timber,cftfrac,fpc,file = dataFile)
    }else{
      print(paste0("Data file (",dataFile,") already exists, old file renamed to: ",dataFile,"_sav"))
      file.rename(dataFile, paste0(dataFile,"_sav"))
      save(ynpp_potential,ynpp,pftnpp_cft,pftnpp_nat,pftnpp_grasslands,pftnpp_bioenergy,harvest_cft,rharvest_cft,fire,timber,cftfrac,fpc,file = dataFile)
    }
  }
  # return variables to calcHANPP
  return(list(ynpp_potential = ynpp_potential, ynpp = ynpp, pftnpp_cft = pftnpp_cft,
              pftnpp_nat = pftnpp_nat, pftnpp_grasslands = pftnpp_grasslands,
              pftnpp_bioenergy = pftnpp_bioenergy, harvest_cft = harvest_cft,
              rharvest_cft = rharvest_cft, fire = fire, timber = timber,
              cftfrac = cftfrac, fpc = fpc))

} # end of readHANPPdata

#' Calculate the ecosystem change metric gamma between 2 simulations/timesteps
#'
#' Function to calculate the ecosystem change metric gamma according
#' to Sykes (1999), Heyder (2011), and Ostberg (2015,2018).
#' This is a reformulated version in R, not producing 100% similar values
#' than the C/bash version from Ostberg 2018, but following their methodology.
#'
#' @param inFol_pnv folder of pnv reference run
#' @param inFol_lu folder of landuse scenario run
#' @param startyr first year of simulations
#' @param stopyr last year of simulations
#' @param gridbased logical are pft outputs gridbased or pft-based?
#' @param pftbands number of natural plant functional types (== bands in fpc - 1)
#' @param cftbands number of crop functional types
#' @param p fraction of pasture band to take
#' @param readPreviouslySavedData flag whether to read previously saved data
#'        instead of reading it in from output files (default F)
#' @param saveDataFile whether to save input data to file (default T)
#' @param dataFile file to save computed hanpp data to (default NULL)
#' @param ncells number of cells in lpjml grid
#' @param fileType type of output files - one of "clm", "nc", "meta" (default: "clm")
#' @param headersize headersize of the output files (default 0)
#' @param varnames data.frame with names of output files -- can be specified to account for variable file names (default NULL -- standard names are used)
#'
#' @return list data object containing arrays of ...
#'
#' @examples
#' \dontrun{
#' }
#' @export
calcHANPP <- function(inFol_lu, inFol_pnv, startyr, stopyr, gridbased = T,
                      pftbands = 11, cftbands = 32, p = 1, readPreviouslySavedData = F,
                      saveDataFile = T, dataFile, ncells = 67420, fileType = "clm",
                      headersize = 0, varnames){
  # reading required data
  if (readPreviouslySavedData) {
    if (file.exists(dataFile)) {
      print(paste0("Reading in data from previously saved data file"))
      load(dataFile)
    }else{
      stop(paste0("dataFile: '",dataFile,"' does not exist but is required since reading is set to FALSE."))
    }
  }else{
    print(paste0("Reading in data from outputs"))
    readVars <- readHANPPdata(inFol_lu = inFol_lu,inFol_pnv = inFol_pnv,startyr = startyr,
                              stopyr = stopyr, gridbased = gridbased, pftbands = pftbands,
                              cftbands = cftbands,dataFile = dataFile, ncells = ncells,
                              fileType = fileType, headersize = headersize,
                              varnames = varnames, saveDataFile = saveDataFile)
    ynpp_potential <- readVars$ynpp_potential
    ynpp <- readVars$ynpp
    pftnpp_cft <- readVars$pftnpp_cft
    pftnpp_nat <- readVars$pftnpp_nat
    pftnpp_grasslands <- readVars$pftnpp_grasslands
    pftnpp_bioenergy <- readVars$pftnpp_bioenergy
    harvest_cft <- readVars$harvest_cft
    rharvest_cft <- readVars$rharvest_cft
    fire <- readVars$fire
    timber <- readVars$timber
    cftfrac <- readVars$cftfrac
    fpc <- readVars$fpc
  }

  print(paste0("Calculating data"))
  npp_act_overtime <- colSums(ynpp*cellarea)/10^15 # from gC/m2 to GtC
  lastyr <- startyr + length(npp_act_overtime) - 1
  npp_pot_overtime <- colSums(ynpp_potential*cellarea)/10^15 # from gC/m2 to GtC
  npp_eco_overtime <- colSums(pftnpp_nat*cellarea)/10^15 # from gC/m2 to GtC
  harvest_cft_overtime <- colSums(harvest_cft*cellarea)/10^15 # from gC/m2 to GtC
  rharvest_cft_overtime <- colSums(rharvest_cft*cellarea)/10^15 # from gC/m2 to GtC
  timber_harvest_overtime <- colSums(timber*cellarea)/10^15 # from gC/m2 to GtC
  fire_overtime <- colSums(fire*cellarea)/10^15 # from gC/m2 to GtC
  #harvest_grassland_overtime=colSums(harvest_grasslands*cellarea)/10^15 # from gC/m2 to GtC
  npp_harv_overtime <- colSums((pftnpp_cft + pftnpp_grasslands*p)*cellarea)/10^15 #from gC/m2 to GtC
  npp_bioenergy_overtime <- colSums(pftnpp_bioenergy*cellarea)/10^15 #from gC/m2 to GtC
  npp_luc_overtime <- npp_pot_overtime - npp_act_overtime
  hanpp_overtime <- harvest_cft_overtime + rharvest_cft_overtime + timber_harvest_overtime + fire_overtime + npp_luc_overtime
  hanpp_overtime_perc_piref <- hanpp_overtime/mean(npp_pot_overtime[1:10])*100
  hanpp <- harvest_cft + rharvest_cft + timber + fire + ynpp_potential - ynpp
  hanpp_perc <- hanpp/ynpp_potential*100 #actual NPPpot as ref
  hanpp_perc_piref <- hanpp/rowMeans(ynpp_potential[,1:10])*100 # NPPpi as ref
  return(list(hanpp_overtime = hanpp_overtime, hanpp = hanpp, hanpp_perc = hanpp_perc,
              hanpp_overtime_perc_piref = hanpp_overtime_perc_piref,
              hanpp_perc_piref = hanpp_perc_piref, ynpp_potential = ynpp_potential,
              npp_act_overtime = npp_act_overtime, npp_pot_overtime = npp_pot_overtime,
              npp_eco_overtime = npp_eco_overtime, harvest_cft_overtime = harvest_cft_overtime,
              rharvest_cft_overtime = rharvest_cft_overtime, fire_overtime = fire_overtime,
              timber_harvest_overtime = timber_harvest_overtime,
              npp_bioenergy_overtime = npp_bioenergy_overtime, npp_luc_overtime = npp_luc_overtime))
}
plotHANPP <- function(hanppData, outFol, plotyears, minVal, maxVal, legendpos, startyr, mapyear, highlightyear, eps = FALSE){
  mapindex <- mapyear - startyr
  print(paste0("Plotting HANPP figures"))
  dir.create(file.path(outFol),showWarnings = F)
  plotGlobalWlin(data = hanppData$hanpp[,mapindex], file = paste0(outFol,"HANPP",mapyear,"_absolute.png"),
                 title = paste0("HANPP in ",mapyear), min = 0, max = 1000, legendtitle = "GtC", legYes = T, onlyPos = F, eps = eps)
  plotNPPovertime(file = paste0(outFol,"HANPP_overtime_LPJmL_",plotyears[1],"-",plotyears[2],".png"),
                  dataAct = hanppData$npp_act_overtime, dataEco = hanppData$npp_eco_overtime,
                  dataPot = hanppData$npp_pot_overtime, dataBE = hanppData$npp_bioenergy_overtime,
                  dataHANPP = hanppData$hanpp_overtime, dataHarvNPP = hanppData$npp_harv_overtime,
                  dataHANPPperc = hanppData$hanpp_overtime_perc_piref,
                  firstyr = startyr, plotyrs = plotyears, minVal = minVal,
                  legendpos = legendpos,maxVal = maxVal, eps = eps,highlightyrs = highlightyear)
  plotHANPPmap(data = hanppData$hanpp_perc[,mapindex], file = paste0(outFol,"HANPP",mapyear,"_LPJmL.png"),
               title = "",legendtitle = "% of NPPpot", eps = eps)
  plotHANPPmap(data = hanppData$hanpp_perc_piref[,mapindex], file = paste0(outFol,"HANPP_piref_",mapyear,"_LPJmL.png"),
               title = "",legendtitle = "% of NPPpi", eps = eps)
  plotNPPmap(data = hanppData$ynpp[,mapindex],file = paste0(outFol,"NPP",mapyear,"_LPJmL.png"),
             title = paste0("NPP ",mapyear),legendtitle = "gC/m2", eps = eps)
} # end of plotHANPP
