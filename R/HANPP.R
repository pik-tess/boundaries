# written by Fabian Stenzel
# 2022 - stenzel@pik-potsdam.de

# todo: optional reading of monthly netcdf CFT outputs

#' Plot global map of MCOL to file
#'
#' Plot global map of MCOL to file with legend colors similar to Haberl et al. 2007
#'
#' @param data array containing MCOL percentage value for each gridcell
#' @param file character string for location/file to save plot to
#' @param plotyears range of years to plot over time
#' @param title character string title for plot
#' @param legendtitle character string legend title
#' @param zeroThreshold smallest value to be distinguished from 0 in legend,
#'        both for negative and positive values (default: 0.1)
#' @param eps write eps file instead of PNG (boolean) - (default: FALSE)
#'
#' @return none
#'
#' @examples
#' \dontrun{
#' }
#' @export
plotMCOLmap <- function(data, file, title, legendtitle, zeroThreshold = 0.1, eps = FALSE){
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

#' Plot absolute MCOL, overtime, maps, and npp into given folder
#'
#' Plot to file a comparison over time of global sums of MCOL, NPPpot, NPPeco,
#' and NPPact, with legend similar to Krausmann et al. 2013
#'
#' @param mcolData mcol data list object (returned from calcMCOL) containing
#'                  mcol, npp_eco_overtime, npp_act_overtime, npp_pot_overtime,
#'                  npp_bioenergy_overtime, mcol_overtime, npp_harv_overtime,
#'                  mcol_overtime_perc_piref, mcol_perc, mcol_perc_piref, ynpp
#' @param file character string for location/file to save plot to
#' @param firstyr first year of mcol object
#' @param plotyears range of years to plot over time
#' @param minVal y-axis minimum value for plot over time (default: 0)
#' @param maxVal y-axis maximum value for plot over time (default: 100)
#' @param legendpos position of legend (default: "topleft")
#' @param mapyear year to plot mcol map for
#' @param mapyear_buffer +- years around mapyear to average mcol
#'                       (make sure these years exist in mcol data)
#' @param highlightyrs year(s) that should be highlighted in overtime plot (default: 2000)
#' @param ref reference period for mcol ("pi" or "act"), to either use
#'        mcolData$mcol_overtime_perc_piref or mcolData$mcol_overtime
#' @param eps write plots as eps, instead of png (default = FALSE)
#'
#' @return none
#'
#' @examples
#' \dontrun{
#' }
#' @export
plotMCOLovertime <- function(mcolData, file, firstyr, plotyrs, highlightyrs = 2000, minVal = 0,
                            maxVal = 100, legendpos = "topleft", ext = FALSE, eps = FALSE){
  lastyr = firstyr + length(mcolData$npp_act_overtime) - 1
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
  plot(x=seq(firstyr,lastyr,1),y=mcolData$npp_pot_overtime,ylab="GtC/yr",xlab="Year",xlim=plotyrs,
       ylim=c(minVal, maxVal),type = "l",col=colz[1],xaxs="i",yaxs="i")
  lines(x=seq(firstyr,lastyr,1),y=mcolData$npp_act_overtime,type = "l",col=colz[2])
  lines(x=seq(firstyr,lastyr,1),y=mcolData$npp_eco_overtime,type = "l",col=colz[3])

  par(bty="n",oma=c(0,0,0,0),mar=c(4,5,1,3), new = T)
  if (ref == "pi") {
    plot(x=seq(firstyr,lastyr,1),y=mcolData$mcol_overtime_perc_piref,ylab="",xlab="",xlim=plotyrs,
         ylim=c(10, 30),type = "l",col=colz[4],xaxs="i", yaxs="i", axes = F)
  } else if (ref == "act") {
    plot(x=seq(firstyr,lastyr,1),y=mcolData$mcol_overtime,ylab="",xlab="",xlim=plotyrs,
         ylim=c(10, 30),type = "l",col=colz[4],xaxs="i", yaxs="i", axes = F)
  }else stop(paste0("Unknown value for parameter ref: ",ref," - Aborting."))

  axis(side = 4, col = colz[4],col.axis = colz[4])
  mtext(text = "%", col=colz[4], side = 4,line = 2)

  if (!is.null(highlightyrs)){
    for (y in highlightyrs){
      lines(x=c(y,y),y=c(minVal,maxVal),col="grey40")
    }
  }
  legend(legendpos,legend = c("NPPpot (PNV)","NPPact (landuse)","NPPeco","M-COL [% NPPpi]"),col=colz[1:4] ,lty=1,cex=1)
  dev.off()
}
#' Calculate the ecosystem change metric gamma between 2 simulations/timesteps
#'
#' Function to calculate the ecosystem change metric gamma according
#' to Sykes (1999), Heyder (2011), and Ostberg (2015,2018).
#' This is a reformulated version in R, not producing 100% similar values
#' than the C/bash version from Ostberg 2018, but following their methodology.
#'
#' @param inFol_lu folder of landuse scenario run
#' @param inFol_pnv folder of pnv reference run
#' @param startyr first year of simulations
#' @param stopyr last year of simulations
#' @param gridbased logical are pft outputs gridbased or pft-based?
#' @param pftbands number of natural plant functional types (== bands in fpc - 1)
#' @param cftbands number of crop functional types
#' @param headersize headersize of the output files (default 0)
#' @param saveDataFile whether to save input data to file (default T)
#' @param dataFile file to save read in input data to
#' @param ncells number of cells in lpjml grid
#' @param fileType type of output files - one of "clm", "nc", "meta" (default: "clm")
#' @param varnames data.frame with names (outname) and timestep of output files
#'                 see code of calcMCOL() for an example with default values
#'
#' @return list data object containing arrays of ...
#'
#' @examples
#' \dontrun{
#' }
#' @export
readMCOLdata <- function(inFol_lu, inFol_pnv, startyr, stopyr, gridbased = T,
                          pftbands = 11, cftbands = 32, headersize = 0,
                          saveDataFile = FALSE, dataFile, ncells = 67420,
                          fileType = "clm", varnames) {
  totalbands = (pftbands + cftbands)
  if (fileType == "clm") {
    nppFile = paste0(inFol_lu, varnames["npp","outname"])
    print(paste0("reading NPP file: ",nppFile))
    if (varnames["npp","timestep"] == "Y") {
      ynpp = drop(lpjmliotools::readYearly(inFile = nppFile,startyear = startyr, ncells = ncells,
                    stopyear = stopyr,size = 4, headersize = headersize,
                    getyearstart = startyr,getyearstop = stopyr))# gC/m2
    }else if (varnames["npp","timestep"] == "M") {
      npp = lpjmliotools::readMonthly(inFile = nppFile, startyear = startyr, stopyear = stopyr,
                size = 4, ncells = ncells, headersize = headersize,
                getyearstart = startyr,getyearstop = stopyr)# gC/m2
      ynpp = apply(npp, c(1,4), sum, drop = T) #gC/m2
    }
    pft_nppFile = paste0(inFol_lu,varnames["pft_npp","outname"])
    print(paste0("reading PFT_NPP file: ",pft_nppFile))
    if (varnames["pft_npp","timestep"] == "Y") {
      pftnpp = drop(lpjmliotools::readCFToutput(inFile = pft_nppFile,startyear = startyr,
                  stopyear = stopyr,size = 4, headersize = headersize, ncells = ncells,
                  bands=totalbands,getyearstart = startyr,getyearstop = stopyr))
    }else if (varnames["pft_npp","timestep"] == "M") {
      mpftnpp = lpjmliotools::readMonthlyCFToutput(inFile = pft_nppFile,startyear = startyr,
                    ncells = ncells, stopyear = stopyr,size = 4,
                    headersize = headersize, bands=totalbands,
                    getyearstart = startyr,getyearstop = stopyr)
      pftnpp = apply(mpftnpp, c(1,2,4), sum, drop = T) #gC/m2
    }
    harvfile = paste0(inFol_lu,varnames["pft_harvest","outname"])
    rharvfile = paste0(inFol_lu,varnames["pft_rharvest","outname"])
    print(paste0("reading harvest file: ",harvfile))
    if (varnames["pft_harvest","timestep"] == "Y") {
      harvest = drop(lpjmliotools::readCFToutput(inFile = harvfile,startyear = startyr,
                                  stopyear = stopyr,size = 4, headersize = headersize, ncells = ncells,
                                  bands=cftbands,getyearstart = startyr,getyearstop = stopyr))
    }else if (varnames["pft_harvest","timestep"] == "M") {
      mharvest = lpjmliotools::readMonthlyCFToutput(inFile = harvfile,startyear = startyr,
                                     ncells = ncells, stopyear = stopyr,size = 4,
                                     headersize = headersize, bands=cftbands,
                                     getyearstart = startyr,getyearstop = stopyr)
      harvest = apply(mharvest, c(1,2,4), sum, drop = T) #gC/m2
    }
    print(paste0("reading rharvest file: ",rharvfile))
    if (varnames["pft_rharvest","timestep"] == "Y") {
      rharvest = drop(lpjmliotools::readCFToutput(inFile = rharvfile,startyear = startyr,
                  stopyear = stopyr,size = 4, headersize = headersize, ncells = ncells,
                  bands=cftbands,getyearstart = startyr,getyearstop = stopyr))
    }else if (varnames["pft_rharvest","timestep"] == "M") {
      mrharvest = lpjmliotools::readMonthlyCFToutput(inFile = rharvfile,startyear = startyr,
                     ncells = ncells, stopyear = stopyr,size = 4, headersize = headersize,
                     bands=cftbands, getyearstart = startyr,getyearstop = stopyr)
      rharvest = apply(mrharvest, c(1,2,4), sum, drop = T) #gC/m2
    }
    timberFile = paste0(inFol_lu,varnames["timber_harvest","outname"])
    print(paste0("reading timber file: ",timberFile))
    if (varnames["timber_harvest","timestep"] == "Y") {
      timber = drop(lpjmliotools::readYearly(inFile = timberFile,startyear = startyr, ncells = ncells,
                             stopyear = stopyr,size = 4, headersize = headersize,
                             getyearstart = startyr,getyearstop = stopyr))# gC/m2
    }else if (varnames["timber_harvest","timestep"] == "M") {
      mtimber = lpjmliotools::readMonthly(inFile = timberFile, startyear = startyr, stopyear = stopyr,
                        size = 4, ncells = ncells, headersize = headersize,
                        getyearstart = startyr,getyearstop = stopyr)# gC/m2
      timber = apply(mtimber, c(1,4), sum, drop = T) #gC/m2
    }
    fireFile = paste0(inFol_lu,varnames["firec","outname"])
    print(paste0("reading fire file: ",fireFile))
    if (varnames["firec","timestep"] == "Y") {
      fire = drop(lpjmliotools::readYearly(inFile = fireFile,startyear = startyr, ncells = ncells,
                               stopyear = stopyr,size = 4, headersize = headersize,
                               getyearstart = startyr,getyearstop = stopyr))# gC/m2
    }else if (varnames["firec","timestep"] == "M") {
      mfire = lpjmliotools::readMonthly(inFile = fireFile, startyear = startyr, stopyear = stopyr,
                            size = 4, ncells = ncells, headersize = headersize,
                            getyearstart = startyr,getyearstop = stopyr)# gC/m2
      fire = apply(mfire, c(1,4), sum, drop = T) #gC/m2
    }
    cftFile = paste0(inFol_lu,varnames["cftfrac","outname"])
    print(paste0("reading cftfrac file: ",cftFile))
    if (varnames["cftfrac","timestep"] == "Y") {
      cftfrac = drop(lpjmliotools::readCFToutput(inFile = cftFile,startyear = startyr,
                                    stopyear = stopyr,size = 4, headersize = headersize, ncells = ncells,
                                    bands=cftbands, getyearstart = startyr,getyearstop = stopyr))
    }else if (varnames["cftfrac","timestep"] == "M") {
      mcftfrac = lpjmliotools::readMonthlyCFToutput(inFile = cftFile,startyear = startyr,
                                       ncells = ncells, stopyear = stopyr,size = 4, headersize = headersize,
                                       bands=cftbands, getyearstart = startyr,getyearstop = stopyr)
      cftfrac = apply(mcftfrac, c(1,2,4), sum, drop = T) #gC/m2
    }
    pnv_nppFile = paste0(inFol_pnv,varnames["npp","outname"])
    print(paste0("reading PNV NPP file: ",pnv_nppFile))
    if (varnames["npp","timestep"] == "Y") {
      ynpp_potential = drop(lpjmliotools::readYearly(inFile = pnv_nppFile,startyear = startyr,
                          stopyear = stopyr,size = 4, headersize = headersize,
                          getyearstart = startyr,getyearstop = stopyr))# gC/m2
    }else if (varnames["npp","timestep"] == "M") {
      npp_potential = lpjmliotools::readMonthly(inFile = pnv_nppFile,startyear = startyr,
                        stopyear = stopyr,size = 4, headersize = headersize,
                        getyearstart = startyr,getyearstop = stopyr)# gC/m2
      ynpp_potential = apply(npp_potential, c(1,4), sum, drop = T) #gC/m2
    }
    fpcFile = paste0(inFol_lu,varnames["fpc","outname"])
    print(paste0("reading FPC file: ",fpcFile))
    if (varnames["fpc","timestep"] == "Y") {
      fpc = drop(lpjmliotools::readCFToutput(inFile = fpcFile,startyear = startyr,
                                   stopyear = stopyr,size = 4, headersize = headersize, ncells = ncells,
                                   bands = (pftbands + 1),getyearstart = startyr,getyearstop = stopyr))
    }else if (varnames["fpc","timestep"] == "M") {
      mfpc = lpjmliotools::readMonthlyCFToutput(inFile = fpcFile,startyear = startyr, ncells = ncells,
                  stopyear = stopyr,size = 4, headersize = headersize,
                  bands = (pftbands + 1), getyearstart = startyr,getyearstop = stopyr)
      fpc = apply(mfpc, c(1,2,4), sum, drop = T) #gC/m2
    }
  }else if (fileType == "nc") {
    nppFile = paste0(inFol_lu,varnames["pft_npp","outname"])
    print(paste0("reading NPP file: ",nppFile))
    npp = lpjmliotools::netcdfMonthly2lpjarray(ncInFile = nppFile, lon = lon, lat = lat)# gC/m2
    print(paste0("reading PFT_NPP file: ",pft_nppFile))
    pftnpp = lpjmliotools::netcdfCFT2lpjarray(ncInFile = pft_nppFile,lon = lon, lat = lat)
    print(paste0("reading cftfrac file: ",cftFile))
    cftfrac = lpjmliotools::netcdfCFT2lpjarray(ncInFile = cftFile,lon = lon, lat = lat)
    print(paste0("reading PNV NPP file: ",pnv_nppFile))
    npp_potential = lpjmliotools::netcdfMonthly2lpjarray(ncInFile = pnv_nppFile,lon = lon, lat = lat)# gC/m2
    print(paste0("reading FPC file: ",fpcFile))
    fpc = lpjmliotools::netcdfCFT2lpjarray(ncInFile = fpcFile,var = "fpc",lon = lon, lat = lat)

    ynpp = apply(npp, c(1,3), sum) #gC/m2
    ynpp_potential = apply(npp_potential, c(1,3), sum) #gC/m2

    nppFile = paste0(inFol_lu, varnames["npp","outname"])
    print(paste0("reading NPP file: ",nppFile))
    if (varnames["npp","timestep"] == "Y") {
      ynpp = drop(lpjmliotools::netcdfYearly2lpjarray(ncInFile = nppFile, lon = lon, lat = lat))# gC/m2
    }else if (varnames["npp","timestep"] == "M") {
      npp = lpjmliotools::netcdfMonthly2lpjarray(ncInFile = nppFile, lon = lon, lat = lat)# gC/m2
      ynpp = apply(npp, c(1,4), sum, drop = T) #gC/m2
    }
    pft_nppFile = paste0(inFol_lu,varnames["pft_npp","outname"])
    print(paste0("reading PFT_NPP file: ",pft_nppFile))
    if (varnames["pft_npp","timestep"] == "Y") {
      pftnpp = drop(lpjmliotools::netcdfCFT2lpjarray(inFile = pft_nppFile, lon = lon, lat = lat))
    }else if (varnames["pft_npp","timestep"] == "M") {
      stop("Monthly netcdf output for pft_npp currently not supported.")
    }
    harvfile = paste0(inFol_lu,varnames["pft_harvest","outname"])
    rharvfile = paste0(inFol_lu,varnames["pft_rharvest","outname"])
    print(paste0("reading harvest file: ",harvfile))
    if (varnames["pft_harvest","timestep"] == "Y") {
      harvest = drop(lpjmliotools::netcdfCFT2lpjarray(inFile = harvfile, lon = lon, lat = lat))
    }else if (varnames["pft_harvest","timestep"] == "M") {
      stop("Monthly netcdf output for pft_harvest currently not supported.")
    }
    print(paste0("reading rharvest file: ",rharvfile))
    if (varnames["pft_rharvest","timestep"] == "Y") {
      rharvest = drop(lpjmliotools::netcdfCFT2lpjarray(inFile = rharvfile, lon = lon, lat = lat))
    }else if (varnames["pft_rharvest","timestep"] == "M") {
      stop("Monthly netcdf output for pft_rharvest currently not supported.")
    }
    timberFile = paste0(inFol_lu,varnames["timber_harvest","outname"])
    print(paste0("reading timber file: ",timberFile))
    if (varnames["timber_harvest","timestep"] == "Y") {
      timber = drop(lpjmliotools::netcdfYearly2lpjarray(ncInFile = timberFile, lon = lon, lat = lat))# gC/m2
    }else if (varnames["timber_harvest","timestep"] == "M") {
      mtimber = lpjmliotools::netcdfMonthly2lpjarray(ncInFile = timberFile, lon = lon, lat = lat)# gC/m2
      timber = apply(mtimber, c(1,4), sum, drop = T) #gC/m2
    }
    fireFile = paste0(inFol_lu,varnames["firec","outname"])
    print(paste0("reading fire file: ",fireFile))
    if (varnames["firec","timestep"] == "Y") {
      fire = drop(lpjmliotools::netcdfYearly2lpjarray(ncInFile = fireFile, lon = lon, lat = lat))# gC/m2
    }else if (varnames["firec","timestep"] == "M") {
      mfire = lpjmliotools::netcdfMonthly2lpjarray(ncInFile = fireFile, lon = lon, lat = lat)# gC/m2
      fire = apply(mfire, c(1,4), sum, drop = T) #gC/m2
    }
    cftFile = paste0(inFol_lu,varnames["cftfrac","outname"])
    print(paste0("reading cftfrac file: ",cftFile))
    if (varnames["cftfrac","timestep"] == "Y") {
      cftfrac = drop(lpjmliotools::netcdfCFT2lpjarray(inFile = cftFile, lon = lon, lat = lat))
    }else if (varnames["cftfrac","timestep"] == "M") {
      stop("Monthly netcdf output for cftfrac currently not supported.")
    }
    pnv_nppFile = paste0(inFol_pnv,varnames["npp","outname"])
    print(paste0("reading PNV NPP file: ",pnv_nppFile))
    if (varnames["npp","timestep"] == "Y") {
      ynpp_potential = drop(lpjmliotools::netcdfYearly2lpjarray(ncInFile = pnv_nppFile, lon = lon, lat = lat))# gC/m2
    }else if (varnames["npp","timestep"] == "M") {
      npp_potential = lpjmliotools::netcdfMonthly2lpjarray(ncInFile = pnv_nppFile, lon = lon, lat = lat)# gC/m2
      ynpp_potential = apply(npp_potential, c(1,4), sum, drop = T) #gC/m2
    }
    fpcFile = paste0(inFol_lu,varnames["fpc","outname"])
    print(paste0("reading FPC file: ",fpcFile))
    if (varnames["fpc","timestep"] == "Y") {
      fpc = drop(lpjmliotools::netcdfCFT2lpjarray(inFile = fpcFile, lon = lon, lat = lat))
    }else if (varnames["fpc","timestep"] == "M") {
      stop("Monthly netcdf output for fpc currently not supported.")
    }
  }else{
    print(paste0("Unrecognized file type (",fileType,")"))
    break
  }
  if (gridbased) {
    pftnpp_grasslands <- apply(pftnpp[, c(totalbands - c(18,2)), ],c(1,3),sum) #gC/m2
    pftnpp_cft <- apply(pftnpp[,-c(1:(pftbands), totalbands - 18, totalbands - 2), ], c(1,3), sum) #gC/m2
    pftnpp_bioenergy <- apply(pftnpp[, c(totalbands - c(17,16,1,0)), ], c(1,3), sum) #gC/m2
    harvest_grasslands <- apply(harvest[,c(cftbands - c(18,2)),],c(1,3),sum) #gC/m2 only from grassland bands
    harvest_bioenergy <- apply(harvest[,c(cftbands - c(17,16,1,0)),],c(1,3),sum) #gC/m2 only from bioenergy bands
    harvest_cft <- apply(harvest[,-c(cftbands - c(18,17,16,2,1,0)),], c(1,3), sum) #gC/m2 not from grassland and bioenergy bands
    rharvest_cft <- apply(rharvest[,-c(cftbands - c(18,17,16,2,1,0)),], c(1,3), sum) #gC/m2 not from grassland and bioenergy bands
    pftnpp_nat <- apply(pftnpp[,1:(pftbands),], c(1,3), sum) #gC/m2
  }else{# todo: complete this part
    pftnpp_cft <- apply(pftnpp[,(totalbands - cftbands + 1):totalbands,]*cftfrac,c(1,3),sum) #gC/m2
    pftnpp_nat <- apply(pftnpp[,1:(pftbands),],c(1,3),sum)*fpc[,1,] #gC/m2
  }
  if (saveDataFile) {
    if (!file.exists(dataFile) ) {
      print(paste0("Writing data file: ",dataFile))
      save(ynpp_potential,ynpp,pftnpp_cft,pftnpp_nat,pftnpp_grasslands,pftnpp_bioenergy,
           harvest_cft,rharvest_cft,fire,timber,cftfrac,fpc,harvest_grasslands,
           harvest_bioenergy,file = dataFile)
    }else{
      print(paste0("Data file (",dataFile,") already exists, old file renamed to: ",dataFile,"_sav"))
      file.rename(dataFile, paste0(dataFile,"_sav"))
      save(ynpp_potential,ynpp,pftnpp_cft,pftnpp_nat,pftnpp_grasslands,pftnpp_bioenergy,
           harvest_cft,rharvest_cft,fire,timber,cftfrac,fpc,harvest_grasslands,
           harvest_bioenergy,file = dataFile)
    }
  }
  # return variables to calcMCOL
  return(list(ynpp_potential = ynpp_potential, ynpp = ynpp, pftnpp_cft = pftnpp_cft,
              pftnpp_nat = pftnpp_nat, pftnpp_grasslands = pftnpp_grasslands,
              pftnpp_bioenergy = pftnpp_bioenergy, harvest_cft = harvest_cft,
              rharvest_cft = rharvest_cft, fire = fire, timber = timber,
              cftfrac = cftfrac, fpc = fpc, harvest_grasslands = harvest_grasslands,
              harvest_bioenergy = harvest_bioenergy))

} # end of readMCOLdata

#' Calculate the ecosystem change metric gamma between 2 simulations/timesteps
#'
#' Function to calculate the ecosystem change metric gamma according
#' to Sykes (1999), Heyder (2011), and Ostberg (2015,2018).
#' This is a reformulated version in R, not producing 100% similar values
#' than the C/bash version from Ostberg 2018, but following their methodology.
#'
#' @param inFol_lu folder of landuse scenario run
#' @param inFol_pnv folder of pnv reference run
#' @param startyr first year of simulations
#' @param stopyr last year of simulations
#' @param gridbased logical are pft outputs gridbased or pft-based?
#' @param npp_threshold lower threshold for npp (to mask out non-lu areas
#'        according to Haberl et al. 2007). Below MCOL will be set to 0.
#'        (default: 20 gC/m2)
#' @param pftbands number of natural plant functional types (== bands in fpc - 1)
#' @param cftbands number of crop functional types
#' @param readPreviouslySavedData flag whether to read previously saved data
#'        instead of reading it in from output files (default FALSE)
#' @param saveDataFile whether to save input data to file (default FALSE)
#' @param dataFile file to save/read input data to/from (default NULL)
#' @param ncells number of cells in lpjml grid
#' @param fileType type of output files - one of "clm", "nc" (default: "clm")
#' @param headersize headersize of the output files (default 0)
#' @param varnames data.frame with names (outname) and timestep of output files
#'        -- can be specified to account for variable file names
#'        (default NULL -- standard names are used)
#' @param grass_scaling whether to scale pasture harvest according to
#'        data given via grass_harvest_file (default T)
#' @param grass_harvest_file file containing grazing data to rescale the
#'        grassland harvests according to Herrero et al. 2013. File contains:
#'        grazing_data list object with $name and $id of 29 world regions, and
#'        $Herrero_2000_kgDM_by_region containing for each of these regions and
#'        mapping_lpj67420_to_grazing_regions array with a mapping between 67420
#'        LPJmL cells and the 29 regions
#'
#' @return list data object containing arrays of ...
#'
#' @examples
#' \dontrun{
#' }
#' @export
calcMCOL <- function(inFol_lu, inFol_pnv, startyr, stopyr, gridbased = T, npp_threshold = 20,
                      pftbands = 11, cftbands = 32, readPreviouslySavedData = FALSE,
                      saveDataFile = FALSE, dataFile = NULL, ncells = 67420, fileType = "clm",
                      headersize = 0, varnames = NULL, grass_scaling = T, grass_harvest_file = "/p/projects/open/Fabian/LPJbox/grazing_data.RData"){
  if (grass_scaling && !file.exists(grass_harvest_file)) stop(paste0("Grass harvest scaling enabled, but grass_harvest_file does not exist in: ",grass_harvest_file))

  if (is.null(varnames)) {
    print("Varnames not given, using standard values, which might not fit this specific configuration. Please check!")
    varnames <- data.frame(row.names = c("npp",    "pft_npp",    "pft_harvest",     "pft_rharvest",     "firec",    "timber_harvest",     "cftfrac",    "fpc"),
                       outname = c(  "npp.bin","pft_npp.bin","pft_harvestc.bin","pft_rharvestc.bin","firec.bin","timber_harvestc.bin","cftfrac.bin","fpc.bin"),
                       timestep = c( "Y",      "Y",          "Y",               "Y",                "Y",        "Y",                  "Y",          "Y"))
  }

  # reading required data
  if (readPreviouslySavedData) {
    if (file.exists(dataFile)) {
      print(paste0("Reading in data from previously saved data file"))
      load(dataFile)
    }else{
      stop(paste0("dataFile: '",dataFile,"' does not exist but is required since reading is set to FALSE."))
    }
    if (saveDataFile) {
      saveDataFile <- FALSE
      print("Both readPreviouslySavedData and saveDataFile have been set to TRUE. Overwriting with the same data does not make sense, saving disabled. ")
    }
  }else{
    print(paste0("Reading in data from outputs"))
    readVars <- readMCOLdata(inFol_lu = inFol_lu,inFol_pnv = inFol_pnv,startyr = startyr,
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
    harvest_bioenergy <- readVars$harvest_bioenergy
    harvest_grasslands <- readVars$harvest_grasslands
    rharvest_cft <- readVars$rharvest_cft
    fire <- readVars$fire
    timber <- readVars$timber
    cftfrac <- readVars$cftfrac
    fpc <- readVars$fpc
  }

  print(paste0("Calculating data"))
  if (grass_scaling) {
    load(grass_harvest_file)
    nregs <- length(grazing_data$name)
    lpj_grass_harvest_region <- array(0,dim = nregs)
    lpj_grass_harvest_2000 <- rowMeans(harvest_grasslands[,(1995-startyr+1):(2005-startyr+1)])*cellarea/1000*2 # from gC/m2 to kgDM
    grassland_scaling_factor_cellwise <- array(1,dim = ncells)
    for (r in 1:nregs) {
      lpj_grass_harvest_region[r] <- sum(lpj_grass_harvest_2000[which(mapping_lpj67420_to_grazing_regions==r)])
    }
    scaling_factor <- grazing_data$Herrero_2000_kgDM_by_region/lpj_grass_harvest_region
    for (r in 1:nregs) {
      grassland_scaling_factor_cellwise[which(mapping_lpj67420_to_grazing_regions==r)] <- scaling_factor[r]
    }
    #plotGlobalMan(data = grassland_scaling_factor_cellwise,file = "~/scaling_factor_lpj_grazing.png",
    #         title = "grassland scaling factor",brks = c(seq(0,1,length.out = 11),seq(2,10,length.out = 9),60),
    #         palette = "RdBu",legendtitle = "",legYes = T,eps=F)
    harvest_grasslands <- harvest_grasslands*rep(grassland_scaling_factor_cellwise,
                                                 times = length(harvest_grasslands[1,]))
  }

  npp_act_overtime <- colSums(ynpp*cellarea)/10^15 # from gC/m2 to GtC
  # lastyr <- startyr + length(npp_act_overtime) - 1
  npp_pot_overtime <- colSums(ynpp_potential*cellarea)/10^15 # from gC/m2 to GtC
  npp_eco_overtime <- colSums(pftnpp_nat*cellarea)/10^15 # from gC/m2 to GtC
  #npp_harv_overtime <- colSums((pftnpp_cft + pftnpp_grasslands)*cellarea)/10^15 #from gC/m2 to GtC
  #npp_bioenergy_overtime <- colSums(pftnpp_bioenergy*cellarea)/10^15 #from gC/m2 to GtC
  npp_luc_overtime <- npp_pot_overtime - npp_act_overtime

  harvest_cft_overtime <- colSums(harvest_cft*cellarea)/10^15 # from gC/m2 to GtC
  rharvest_cft_overtime <- colSums(rharvest_cft*cellarea)/10^15 # from gC/m2 to GtC
  harvest_grasslands_overtime <- colSums(harvest_grasslands*cellarea)/10^15 # from gC/m2 to GtC
  harvest_bioenergy_overtime <- colSums(harvest_bioenergy*cellarea)/10^15 # from gC/m2 to GtC

  timber_harvest_overtime <- colSums(timber*cellarea)/10^15 # from gC/m2 to GtC
  fire_overtime <- colSums(fire*cellarea)/10^15 # from gC/m2 to GtC


  mcol_overtime <- harvest_cft_overtime + rharvest_cft_overtime +
                    harvest_grasslands_overtime + harvest_bioenergy_overtime +
                    timber_harvest_overtime + fire_overtime + npp_luc_overtime
  mcol_overtime_perc_piref <- mcol_overtime/mean(npp_pot_overtime[1:10])*100
  mcol <- harvest_cft + rharvest_cft + harvest_grasslands + harvest_bioenergy + timber + fire + ynpp_potential - ynpp
  mcol[ynpp_potential<npp_threshold] <- 0 # set to 0 below lower threshold of NPP
  mcol_perc <- mcol/ynpp_potential*100 #actual NPPpot as ref
  mcol_perc_piref <- mcol/rowMeans(ynpp_potential[,1:10])*100 # NPPpi as ref

  if (grass_scaling){
    return(list(mcol_overtime = mcol_overtime, mcol = mcol, mcol_perc = mcol_perc,
                mcol_overtime_perc_piref = mcol_overtime_perc_piref,
                mcol_perc_piref = mcol_perc_piref, ynpp_potential = ynpp_potential,
                npp_act_overtime = npp_act_overtime, npp_pot_overtime = npp_pot_overtime,
                npp_eco_overtime = npp_eco_overtime, #npp_bioenergy_overtime = npp_bioenergy_overtime,
                harvest_cft_overtime = harvest_cft_overtime, npp_luc_overtime = npp_luc_overtime,
                rharvest_cft_overtime = rharvest_cft_overtime, fire_overtime = fire_overtime,
                timber_harvest_overtime = timber_harvest_overtime, harvest_cft = harvest_cft,
                grassland_scaling_factor_cellwise = grassland_scaling_factor_cellwise))
  }else{
    return(list(mcol_overtime = mcol_overtime, mcol = mcol, mcol_perc = mcol_perc,
                mcol_overtime_perc_piref = mcol_overtime_perc_piref,
                mcol_perc_piref = mcol_perc_piref, ynpp_potential = ynpp_potential,
                npp_act_overtime = npp_act_overtime, npp_pot_overtime = npp_pot_overtime,
                npp_eco_overtime = npp_eco_overtime, #npp_bioenergy_overtime = npp_bioenergy_overtime,
                harvest_cft_overtime = harvest_cft_overtime, npp_luc_overtime = npp_luc_overtime,
                rharvest_cft_overtime = rharvest_cft_overtime, fire_overtime = fire_overtime,
                timber_harvest_overtime = timber_harvest_overtime, harvest_cft = harvest_cft))
  }

}

#' Plot absolute MCOL, overtime, maps, and npp into given folder
#'
#' Wrapper function to plot absolute mcol, overtime, maps, and npp into given folder
#'
#' @param mcolData mcol data list object (returned from calcMCOL) containing
#'                  mcol, npp_eco_overtime, npp_act_overtime, npp_pot_overtime,
#'                  npp_bioenergy_overtime, mcol_overtime, npp_harv_overtime,
#'                  mcol_overtime_perc_piref, mcol_perc, mcol_perc_piref, ynpp
#' @param outFol folder to write into
#' @param plotyears range of years to plot over time
#' @param minVal y-axis minimum value for plot over time
#' @param maxVal y-axis maximum value for plot over time
#' @param legendpos position of legend
#' @param startyr first year of mcolData object
#' @param mapyear year to plot mcol map for
#' @param mapyear_buffer +- years around mapyear to average mcol
#'                       (make sure these years exist in mcolData)
#' @param highlightyear year(s) that should be highlighted in overtime plot
#' @param eps write plots as eps, instead of png (default = FALSE)
#'
#' @return none
#'
#' @examples
#' \dontrun{
#' }
#' @export
plotMCOL <- function(mcolData, outFol, plotyears, minVal, maxVal, legendpos,
                  startyr, mapyear, mapyear_buffer = 5, highlightyear, eps = FALSE){
  mapindex <- mapyear - startyr
  print(paste0("Plotting MCOL figures"))
  dir.create(file.path(outFol),showWarnings = F)
  lpjmliotools::plotGlobal(data = rowMeans(mcolData$mcol[,(mapindex-mapyear_buffer):(mapindex+mapyear_buffer)]),
             file = paste0(outFol,"MCOL",mapyear,"_absolute.png"),
             title = paste0("MCOL in ",mapyear), min = 0, max = 1000,
             legendtitle = "GtC", legYes = T, onlyPos = F, eps = eps, type = "lin")
  plotMCOLovertime(mcolData = mcolData, file = paste0(outFol,"MCOL_overtime_LPJmL_",plotyears[1],"-",plotyears[2],".png"),
                  firstyr = startyr, plotyrs = plotyears, minVal = minVal, ref = "pi",
                  legendpos = legendpos, maxVal = maxVal, eps = eps, highlightyrs = highlightyear)
  plotMCOLmap(data = rowMeans(mcolData$mcol_perc[,(mapindex-mapyear_buffer):(mapindex+mapyear_buffer)]),
               file = paste0(outFol,"MCOL",mapyear,"_LPJmL.png"),legendtitle = "% of NPPpot", eps = eps,
               title = paste0(mapindex-mapyear_buffer, " - ",mapindex+mapyear_buffer) )
  plotMCOLmap(data = rowMeans(mcolData$mcol_perc_piref[,(mapindex-mapyear_buffer):(mapindex+mapyear_buffer)]),
               file = paste0(outFol,"MCOL_piref_",mapyear,"_LPJmL.png"),
               title = paste0(mapindex-mapyear_buffer, " - ",mapindex+mapyear_buffer),legendtitle = "% of NPPpi", eps = eps)
  lpjmliotools::plotGlobalMan(data = rowMeans(mcolData$ynpp[,(mapindex-mapyear_buffer):(mapindex+mapyear_buffer)]),file = paste0(outFol,"NPP",mapyear,"_LPJmL.png"), brks = seq(0,1000,100),
                palette = c("orangered4","orangered","orange","gold","greenyellow","limegreen","green4","darkcyan","darkslategrey","navy"),
                title = paste0("NPP average ",mapindex-mapyear_buffer, "-",mapindex+mapyear_buffer),legendtitle = "gC/m2",legYes = T)
} # end of plotMCOL
