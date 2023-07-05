library(raster)
library(RColorBrewer)
library(maps)
library(scales) #needed for alpha()
library(fields) #for image.plot
library(zoo) #for running mean
library(lpjmliotools)
# ===== Settings =======
fplotspath <- "/home/stenzel/ownCloud/5_PBwater_reworked/plots/"
ncells=67420
nDays=c(31,28,31,30,31,30,31,31,30,31,30,31)
readin=F
basinExportCSV=F
# ===== Functions =======
#source("/home/stenzel/Nextcloud_PIK/7_irrigation_History/data/lpjmliotools_20201104.R")
#show routing-trace
showRoute <- function(ind,routingTable){
  if (routingTable[ind]<1){return(ind)} #can be 0 or -8 -> endcell or nacell
  else{return(c(ind,showRoute(routingTable[ind],routingTable)))}
}
upstreamCells <- function(ind,routingTable){
  return(which(routingTable==ind))
}
nUpstreamCells <- function(ind,routingTable){
  return(length(which(routingTable==ind)))
}
showNetwork <- function(ind,routingTable){#ind is vector
  lind=sapply(ind,nUpstreamCells,routingTable)
  if (sum(lind)<1){
    print(unlist(ind))
    lv[[i]]<<-unlist(ind)
    i<<-0
    } #all branches ended
  else{
    p=unlist(ind)
    if (length(dim(p))>1){ dim(p)=c(length(p[,1]))}
    print(p)
    lv[[i]]<<-p
    i<<-i+1
    showNetwork(unlist(sapply(ind,upstreamCells,routingTable)),routingTable)
  }
}
plotBasin <- function(cell,endcellList,zoom=F,off=20){#ind is vector
  data=array(0,dim=c(ncells))
  e=endcellList[cell]
  data[which(endcellList==e)]=1
  data[cell]=2
  palette <- c("white","blue","red")
  brks <- c(-0.5,0.5,1.5,2.5)
  ra <- raster::raster(ncols=720, nrows=360)
  range <- range(data)
  ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
  ext=c(-180, 180, -60, 90)
  if (zoom){ext=c(lon[cell]-off,lon[cell]+off,lat[cell]-off,lat[cell]+off)}
  extent <- raster::extent(ext)
  par(bty="n",mar=c(0,0,0,0),oma=c(0,0,0,0))
  raster::plot(ra,ext=extent,breaks=brks,col=palette,main="",legend=FALSE,axes=FALSE)
  maps::map('world',add=TRUE,res=0.4, lwd=0.25,ylim=c(-60,90))
}
calcEFRs <- function(discharge,method="VMF"){
  ncells=dim(discharge)[1]
  efrMethodsSet = c("PBpaper","VMF","Q90Q50")
  #method valid?
  if (!method %in% efrMethodsSet){
    print(paste("EFR method not recognized, use one of: ",paste(efrMethodsSet,collapse = " ")))
    return(NA)
  }
  #make sure discharge is a ncells,30year,12month array
  if (!all(dim(discharge) == c(ncells,12,30))){
    print(paste("discharge array has wrong dimension, use c(ncells,nmonths=12,nyears=30)"))
    return(NA)
  }
  efrs=discharge[,,1]*0 #initialize efrs array
  
  #calculate mean monthly flow and mean annual flow
  MMF=apply(discharge,c(1,2),mean)
  quantiles=apply(MMF,c(1),quantile,probs=c(0.5,0.9,0.05,0.95))
  MAF=rep(rowMeans(MMF),times=12)
  Q90=rep(quantiles[2,],times=12)
  Q50=rep(quantiles[1,],times=12)
  dim(MAF)=c(ncells,12)
  dim(Q90)=c(ncells,12)
  dim(Q50)=c(ncells,12)
  remove(quantiles)
  
  if (method==efrMethodsSet[1]){ #"PBpaper" - Steffen et al. 2015
    efrs[MMF<=0.4*MAF]=0.75*MMF[MMF<=0.4*MAF] # low flow months
    efrs[MMF>0.4*MAF & MMF<=0.8*MAF]=0.7*MMF[MMF>0.4*MAF & MMF<=0.8*MAF] # intermediate flow months
    efrs[MMF>0.8*MAF]=0.45*MMF[MMF>0.8*MAF] # high flow months
  }else if (method==efrMethodsSet[2]){ # "VMF" - Pastor et al. 2014
    efrs[MMF<=0.4*MAF]=0.6*MMF[MMF<=0.4*MAF] # low flow months
    efrs[MMF>0.4*MAF & MMF<=0.8*MAF]=0.45*MMF[MMF>0.4*MAF & MMF<=0.8*MAF] # intermediate flow months
    efrs[MMF>0.8*MAF]=0.3*MMF[MMF>0.8*MAF] # high flow months
  }else if (method==efrMethodsSet[3]){ # "Q90Q50" - Pastor et al. 2014
    efrs[MMF<=MAF]=Q90[MMF<=MAF] # low flow months
    efrs[MMF>MAF]=Q50[MMF>MAF] # high flow months
  }else if (method==efrMethodsSet[4]){ # "newPB"
    efrs=Q95-Q05 # low flow months
  }
  
  return(efrs)
}
flowSeason <- function(discharge,method="VMF"){
  ncells=dim(discharge)[1]
  efrMethodsSet = c("PBpaper","VMF","Q90Q50")
  #method valid?
  if (!method %in% efrMethodsSet){
    print(paste("EFR method not recognized, use one of: ",paste(efrMethodsSet,collapse = " ")))
    return(NA)
  }
  #make sure discharge is a ncells,12month array
  if (!all(dim(discharge) == c(ncells,12))){
    print(paste("discharge array has wrong dimension, use c(ncells,nmonths=12)"))
    return(NA)
  }
  fs=discharge*0 #initialize fs array
  
  #calculate mean monthly flow and mean annual flow
  MMF=discharge
  MAF=rep(rowMeans(MMF),times=12)
  dim(MAF)=c(ncells,12)

  if (method==efrMethodsSet[1]){ #"PBpaper" - Steffen et al. 2015
    fs[MMF<=0.4*MAF]="l"
    fs[MMF>0.4*MAF & MMF<=0.8*MAF]="i"
    fs[MMF>0.8*MAF]="h"
  }else if (method==efrMethodsSet[2]){ # "VMF" - Pastor et al. 2014
    fs[MMF<=0.4*MAF]=1#"l"
    fs[MMF>0.4*MAF & MMF<=0.8*MAF]=2#"i"
    fs[MMF>0.8*MAF]=3#"h"
  }else if (method==efrMethodsSet[3]){ # "Q90Q50" - Pastor et al. 2014
    fs[MMF<=MAF]="l"
    fs[MMF>MAF]="h"
  }
  
  return(fs)
}
sumUpPBandEFRs_old <- function(endcell,routing,cellindex,efr,discharge,remotes=F){
  basins=which(routing==0) # indices of ocean draining cells (due to adding 1, -1 becomes 0)
  basinPB=array(0,dim = c(length(basins),2))# [basin#,c(PB,EFR)]
  cellEFR=routing*0
  cellPB=routing*0
  totalEFR=routing*0
  a=1:length(basins)
  if (length(remotes>1)){ 
    efr[remotes==2]=discharge[remotes==2] #outflowcells of remote areas: remotes==2
    efr[remotes==1]=discharge[remotes==1] #remote areas: remotes==1
    #discharge[remotes==1]=0 #remote areas: remotes==1
  }
  for (b in a){# for all basins (basins)
    #print(paste("processing basin:",basins[b]))
    basinCells=which(endcell==basins[b]) #get all cell indices from the basin with endcell basins[b]
    basinCellIndices=cellindex[basinCells] #get all their level indices (how far upstream from endcell?)
    efrsum=0 #reset sum of basin efrs
    pbsum=0 #reset sum of basin pb
    for (i in max(basinCellIndices):1){ #for each level from the basin
      parallelCells=basinCells[which(cellindex[basinCells]==i)]
      for (c in parallelCells){ # for each cell from the same level
        upCells=upstreamCells(c,routing) # get all direct upstream cells
        upCellsEFR=sum(efr[upCells]) # sum of EFRs of direct upstream cells
        upCellsPB=sum(discharge[upCells]-efr[upCells]) # for each cell, the PB is calc. as discharge minus EFRs
        totalEFR[c]=max(efr[c],upCellsEFR)
        cellBalancePB=discharge[c]-efr[c]-upCellsPB # PB (available water) for given cell as discharge minus own EFRs minus water to be used upstream to avoid double counting
        cellBalanceEFR=efr[c]-upCellsEFR # EFRs for given cell as own EFRs minus EFRs already protected upstream to avoid double counting
        if (basins[b]==24) print(paste("Basin:",b,"Level:",i,"Cell:",c,"upCellsEFR:",upCellsEFR,"upCellsPB:",upCellsPB,"discharge:",discharge[c],"efr:",efr[c]))
        cellEFR[c]=cellBalanceEFR
        cellPB[c]=cellBalancePB
        if (cellBalancePB>0){pbsum=pbsum+cellBalancePB}# only add if positive
        if (cellBalanceEFR>0){efrsum=efrsum+cellBalanceEFR}# only add if positive
      }#end for each cell
    }#end for each index level
    basinPB[b,2]=efrsum
    basinPB[b,1]=pbsum
  }#end for all basins
  ef=sum(basinPB[,2])
  pb=sum(basinPB[,1])
  print(paste0("EFRs: ",ef," -- PB: ",pb," -- Total: ",ef+pb," -- Ratio EFR/Total: ",ef/(ef+pb)))
  return(basinPB)
}#end function
sumUpPBandEFRs <- function(endcell,routing,cellindex,efr,discharge,remotes=F,flowseason=F){
  basins=which(routing==0) # indices of ocean draining cells (due to adding 1, -1 becomes 0)
  totalEFR=routing*0
  totalPB=routing*0
  a=1:length(basins)
  if (length(remotes>1)){ 
    efr[remotes==2]=discharge[remotes==2] #outflowcells of remote areas: remotes==2
    efr[remotes==1]=discharge[remotes==1] #remote areas: remotes==1
  }
  if (length(flowseason>1)){ 
    #fs=flowSeason(discharge=discharge,method = "VMF")
    #efr[fs=="h"]=discharge[fs=="h"] #outflowcells of remote areas: remotes==2
  }
  for (b in a){# for all basins (basins)
    basinCells=which(endcell==basins[b]) #get all cell indices from the basin with endcell basins[b]
    basinCellIndices=cellindex[basinCells] #get all their level indices (how far upstream from endcell?)
    for (i in max(basinCellIndices):1){ #for each level from the basin
      parallelCells=basinCells[which(cellindex[basinCells]==i)]
      for (c in parallelCells){ # for each cell from the same level
        upCells=upstreamCells(c,routing) # get all direct upstream cells
        upCellsEFR=sum(totalEFR[upCells]) # sum of EFRs of direct upstream cells
        upCellsPB=sum(totalPB[upCells]) # same for PB
        totalEFR[c]=max(efr[c],upCellsEFR) # EFRs for given cell as own EFRs minus EFRs already protected upstream to avoid double counting
        totalPB[c]=max(0,discharge[c]-totalEFR[c]) # PB (available water) for given cell as discharge minus totally to save EFRs (own+upstream)
      }#end for each cell
    }#end for each index level
  }#end for all basins
  ret=cbind(totalEFR,totalPB)
  colnames(ret)=c("EFR","PB")
  #print(paste("Summary: ","global EFR sum:",round(sum(ret[basins,1]))," -- global PB:",round(sum(ret[basins,2]))," -- Approx. global water: ",round(sum(ret[basins,]))," -- Ratio (EFR/total): ",round(sum(ret[basins,1])/sum(ret[basins,]),3)))
  return(ret)
}#end function
pbEFRstats <- function(efrPBarray,endcells){
  if (length(dim(efrPBarray))==3){
    efr=sum(efrPBarray[endcells,1,])
    pb=sum(efrPBarray[endcells,2,])
  }else if (length(dim(efrPBarray))==2){
    efr=sum(efrPBarray[endcells,1])
    pb=sum(efrPBarray[endcells,2])
  }else{
    print("Dimensions of efrPBarray not suitable - needs to be c(ncells,2,12) or c(ncells,2)")
  }
  total=efr+pb
  print(paste("Summary: ","global PB:",round(pb)," -- global EFR sum:",round(efr)," -- Approx. global water: ",round(total)," -- Ratio (EFR/total): ",round(efr/total,3)))
}
autoReadInput <- function(inFile,getyearstart=-1,getyearstop=-1,manu=F,msize=4,mheadersize=43){
  hdr=readheader(filename=inFile)$header
  print(hdr)
  startyear=hdr[3]
  stopyear=hdr[3]+hdr[4]-1
  ncells=hdr[6]
  nbands=hdr[7]
  if (getyearstart==-1){
    getyearstart=startyear
  }
  if (getyearstop==-1){
    getyearstop=stopyear
  }
  if (hdr[1]==1){#header version 1 
    headersize=36
  }else if (hdr[1]==2){#header version 2
    headersize=43
  }else{ #header version 3
    headersize=51
  }
  if (length(hdr)>10){
    if (hdr[11]==0){
      size=1
      inputType="char"
    }else if(hdr[11]==1){
      size=2
      inputType="integer"
    }else if(hdr[11]==2){
      size=4
      inputType="integer"
    }else if(hdr[11]==3){
      size=4
      inputType="double"
    }else if(hdr[11]==4){
      size=8
      inputType="double"
    }
  }
  if (manu){
    size=msize
    headersize=mheadersize
  }
  if (getyearstop>stopyear){
    stop(paste("unexpected usage: getyearstop (",getyearstop,") larger than stopyear (",stopyear,") -- stopping"))
  }
  if (getyearstart<startyear){
    stop(paste("unexpected usage: getyearstart (",getyearstart,") smaller than startyear (",startyear,") -- stopping"))
  }
  nyears=getyearstop-getyearstart+1
  input <- file(inFile,"rb")
  seek(input, where = headersize+(getyearstart-startyear)*ncells*size*nbands, origin="start")
  if (inputType == "integer"){
    dataIn <- readBin(input,integer(),n = nyears*ncells*nbands, size=size)
  }else if (inputType == "double"){
    dataIn <- readBin(input,double(),n = nyears*ncells*nbands, size=size)
  }else{
    dataIn <- readBin(input,character(),n = nyears*ncells*nbands, size=size)
  }
  close(input)      #remove to save space
  print(paste("nyears=",nyears," nbands=",nbands," ncells=",ncells," size=",size," headersize=",headersize))
  if (nyears==1){
    dim(dataIn) <- c(nbands,ncells)
  }else{
    dim(dataIn) <- c(nbands,ncells,nyears)
  }
  return(dataIn*hdr[["scalar"]])
}
daily2monthly <- function(daily,method="sum"){
  nDays=c(31,28,31,30,31,30,31,31,30,31,30,31)
  if (!length(daily)==365) print("Error - input is not in daily format")
  monthly=array(0,dim=c(12))
  for (m in 1:12){
    if (m==1) fdom=1
    else fdom=sum(nDays[1:(m-1)])+1
    ldom=fdom+nDays[m]-1
    if (method=="sum"){
      monthly[m]=sum(daily[fdom:ldom])
    }else if (method=="mean"){
      monthly[m]=mean(daily[fdom:ldom])
    }else print("Error - method not recognized")
  }
  return(monthly)
}
removeHighFlows <- function(ddischarge,threshold=0.9,period="all"){
  years=dim(ddischarge)[3]
  #make sure discharge is a ncells,365,X array
  if (!all(dim(ddischarge)[1:2] == c(ncells,365))){
    print(paste("discharge array has wrong dimension, use c(ncells,ndays=365)"))
    return(NA)
  }
  if (period=="all"){
    daily_quantiles=apply(ddischarge,c(1),quantile,probs=c(threshold))
    dailyQ90=rep(daily_quantiles,years*365)
    dim(dailyQ90)=c(ncells,365,years)
  }else if (period=="year"){
    daily_quantiles=apply(ddischarge,c(1,3),quantile,probs=c(threshold))
    dailyQ90=rep(daily_quantiles,365)
    dim(dailyQ90)=c(ncells,years,365)
    dailyQ90=aperm(dailyQ90,c(1,3,2))
  }else if (period=="raindays"){
    a=ddischarge
    a[a==0]=NA
    daily_quantiles=apply(a,c(1),quantile,probs=c(threshold),na.rm=T)
    dailyQ90=rep(daily_quantiles,years*365)
    dim(dailyQ90)=c(ncells,365,years)
  }else{
    print(paste("parameter period needs to be set as 'all', 'raindays' or 'year'"))
    return(NA)
  }
  ddischarge[ddischarge>dailyQ90]=0
  return(ddischarge)
}

# ===== Data prep =====  ####
oFol=fplotspath
iFol="/home/stenzel/Nextcloud_PIK/4_PB_transgressions_BE/plots/"
outFile=paste0(oFol,"dataEFR.RData")
if (readin){
  iFol="/p/projects/open/Fabian/runs/IRRHIST/cru_ts3.21/"
  load(file=paste0("/p/projects/open/Fabian/calc/data/lat.RData"))
  load(file=paste0("/p/projects/open/Fabian/calc/data/lon.RData"))
  cellarea <- (111194.9/2)*(111194.9/2)*cos(lat/180*pi) # cellarea in m2
  discharge_2000=rowSums(apply(readMonthly(inFile = paste0(iFol,"historic_noefr_ilim_landuse/mdischarge.bin"),startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1996,headersize = 0,getyearstop = 2005),c(1,2),mean)*rep(nDays,each=67420))/10^3 # 10 year mean from hm3/d to km3/year
  discharge_1685=rowSums(apply(readMonthly(inFile = paste0(iFol,"efrcalc_efr_ilim_nolu/mdischarge.bin"),startyear = 1670,stopyear = 1699,size = 4,getyearstart = 1670,headersize = 0,getyearstop = 1699),c(1,2),mean)*rep(nDays,each=67420))/10^3 # 30 year mean from hm3/d to km3/year
  efrs_1685=rowSums(apply(readMonthly(inFile = paste0(iFol,"historic_efr_ilim_landuse/mefr_threshold.bin"),startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1996,headersize = 0,getyearstop = 2005),c(1,2),mean)*rep(nDays,each=67420))/10^3 # 10 year mean from hm3/d to km3/year
  efr_deficit_2000=rowSums(apply(readMonthly(inFile = paste0(iFol,"historic_efr_ilim_landuse/mefr_deficit.bin"),startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1996,headersize = 0,getyearstop = 2005),c(1,2),mean)*rep(nDays,each=67420))/10^3 # 10 year mean from hm3/d to km3/year
  drainage=autoReadInput(inFile = "/p/projects/open/Fabian/input_VERSION2/drainage.bin",manu=T,msize=4,mheadersize=43)
  routing=drainage[1,]+1 #add 1 since in C indexing starts at 0 but in R at 1
  endcells=which(routing==0) # indices of ocean draining cells (due to adding 1, -1 becomes 0)
  nacells=which(routing==-8) # cells that were not part of STN, LPJmL treats all values below 0 as outflow cells
  #creating cellindex (increasing index upstream from endcell=1) and endcell array for each cell -- not parallelizable
  cellindex=array(0,dim=ncells)
  endcell=array(0,dim=ncells)
  for (c in 1:ncells){
    if (c%%signif(ncells/20,1)==0){print(paste0("indexing done ",round(c/ncells*100,1),"%"))}
    #write(c,stdout())
    route=showRoute(c,routing)
    cellindex[route]=seq(length(route),1,-1)
    endcell[c]=route[length(route)]
  }
  daysinflated=rep(rep(nDays,each=67420),30)
  dim(daysinflated)=c(ncells,12,30)
  discharge_monthly_1976_2005=readMonthly(inFile = paste0(iFol,"historic_noefr_ilim_landuse/mdischarge.bin"),startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1976,headersize = 0,getyearstop = 2005)*daysinflated/10^3 # from hm3/d to km3/month
  discharge_1990_monthly=apply(discharge_monthly_1976_2005,c(1,2),mean)
  discharge_1990=rowSums(discharge_1990_monthly)
  efrs_1990=calcEFRs(discharge = discharge_monthly_1976_2005,method = "PBpaper")
  save(efr_deficit_2000,efrs_1990,discharge_1990,discharge_monthly_1976_2005,discharge_1990_monthly,efrs_1685,lat, lon,cellarea,endcell,endcells,routing,drainage,nacells,cellindex,discharge_2000,discharge_1685,efr,efr_deficit_2000,file=outFile,version = 2)
}else{
  load(file=outFile)
}

load(file=paste0(oFol,"population2006_SSP2.RData")) # pop[67420]
popdens=pop/cellarea*10^4 #cap/km2
#plotGlobalWtoscreen(data = popdens,title = "",pow2min = 4,pow2max = 6,legendtitle = "",legYes = T)

remote=pop*0
remote[popdens<0.001]=1 #for a start define remoteness by popdens below 0.001 cap/km2
load(file=paste0(oFol,"intactForests.RData"))
remote[intactForests==1]=1 #add intact Forest landscapes to remote areas
upRemote=remote*0 #calculate upstream cells of all remote areas
remoteCells=which(remote==1)
for (c in remoteCells){
  upRemote[upstreamCells(c,routing)]=1
}
remote[upRemote==1]=1 #add the upStream cells of remote cells as well (makes computation easier and does not change much)
#par(oma=c(0,0,0,0),mar=c(0,0,0,0))
#plotGlobalWtoscreen(data = remote*10,title = "remote: popdens<0.001",pow2min = 0,pow2max = 4,legendtitle = "",legYes = T)

outflowRemoteCells=remote*0
remoteCells=which(remote==1)
for (c in remoteCells){
  route=showRoute(c,routing) #get the downstream route from this cell 
  inters=intersect(route,remoteCells) #get the portion of the downstream route, that is remote
  outflowRemoteCells[inters[which.min(cellindex[inters])]]=1 #set the last cell of that route (min cellindex) to 1; this is the outflowcell
}
#plotGlobalWtoscreen(data = outflowRemoteCells*10,title = "",pow2min = 0,pow2max = 4,legendtitle = "",legYes = T)
fullRemote=remote
fullRemote[outflowRemoteCells==1]=2 
#plotGlobalWtoscreen(data = fullRemote*5,title = "",pow2min = 0,pow2max = 4,legendtitle = "",legYes = T)

# =============== calc global PB from monthly ISIMIP data  ========================
basinPB = sumUpPBandEFRs(endcell = endcell,routing=routing,cellindex = cellindex,efr=rowSums(efrs_1990),discharge = discharge_1990,remotes=F)
basinPBremote = sumUpPBandEFRs(endcell = endcell,routing=routing,cellindex = cellindex,efr=rowSums(efrs_1990),discharge = discharge_1990,remotes=fullRemote)

basinPBremote_monthly=array(0,dim=c(ncells,2,12))
for (m in 1:12){
  basinPBremote_monthly[,,m]=sumUpPBandEFRs(endcell = endcell,routing=routing,cellindex = cellindex,efr=efrs_1990[,m],discharge = discharge_1990_monthly[,m],remotes=fullRemote)
}
pbEFRstats(basinPB,endcells)
pbEFRstats(basinPBremote,endcells)
pbEFRstats(basinPBremote_monthly,endcells)

basinPBall = sumUpPBandEFRs(endcell = endcell,routing=routing,cellindex = cellindex,efr=rowSums(efrs_1990),discharge = discharge_1990,remotes=F)
plotGlobalWlin(data = basinPBremote[,2]/discharge_1990*100,title = "PB/discharge*100 remote excluded",  file=paste0(fplotspath,"perc_pb_remote.png") ,max = 100,min=0,colPos = "Spectral",legendtitle = "",legYes = T,onlyPos = T,eps = F)
plotGlobalWlin(data = basinPBremote[,1]/discharge_1990*100,title = "EFR/discharge*100 remote excluded",file=paste0(fplotspath,"perc_efr_remote.png"),max = 100,min=0,colPos = "Spectral",legendtitle = "",legYes = T,onlyPos = T,eps = F)
plotGlobalWlin(data = basinPBall[,2]/discharge_1990*100,title = "PB/discharge*100 all cells",  file=paste0(fplotspath,"perc_pb_all.png") ,max = 100,min=0,colPos = "Spectral",legendtitle = "",legYes = T,onlyPos = T,eps = F)
plotGlobalWlin(data = basinPBall[,1]/discharge_1990*100,title = "EFR/discharge*100 all cells",file=paste0(fplotspath,"perc_efr_all.png"),max = 100,min=0,colPos = "Spectral",legendtitle = "",legYes = T,onlyPos = T,eps = F)

efrDiff=basinEFR-pureBasinEFR
sum(efrDiff[efrDiff<0])# rel. wenig
sum(efrDiff[efrDiff>0])# deutlich mehr --> pure>basin

# =============== calc global PB with CRU data and daily discharge  ========================

#discharge_1901_1930_daily=readDaily(inFile = "/p/projects/open/Fabian/runs/IRRHIST/cru_ts3.21/historic_noefr_ilim_landuse/d_discharge.bin",startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1901,headersize = 0,getyearstop = 1930)*24*3600 # m3/s to m3/day
#save(discharge_1901_1930_daily,file="/media/All/PB_water_revised/data.RData")
load(file="/media/All/PB_water_revised/data.RData")
discharge_1901_1930_monthly_full=aperm(apply(X = discharge_1901_1930_daily,MARGIN = c(1,3),FUN = daily2monthly,method="sum")/10^9,c(2,1,3)) # daily2monthly + mean 10y period + convert from m3/month to km3/month
discharge_1901_1930_monthly=apply(discharge_1901_1930_monthly_full,c(1,2),mean)
efrs_1901_1930_monthly=calcEFRs(discharge = discharge_1901_1930_monthly_full,method = "VMF")
basinPBremote_CRU = sumUpPBandEFRs(endcell = endcell,routing=routing,cellindex = cellindex,efr=rowSums(efrs_1901_1930_monthly),discharge = rowSums(discharge_1901_1930_monthly),remotes=fullRemote)
pbEFRstats(basinPBremote_CRU,endcells)
basinPBremote_CRU_remote = sumUpPBandEFRs(endcell = endcell,routing=routing,cellindex = cellindex,efr=rowSums(efrs_1901_1930_monthly),discharge = rowSums(discharge_1901_1930_monthly),remotes=F)
pbEFRstats(basinPBremote_CRU_remote,endcells)

#Q90
discharge_1901_1930_daily_noHF=removeHighFlows(discharge_1901_1930_daily,threshold=0.9,period = "raindays")
high_flows=aperm(apply(X = discharge_1901_1930_daily-discharge_1901_1930_daily_noHF,MARGIN = c(1,3),FUN = daily2monthly,method="sum")/10^9,c(2,1,3)) #from m3 to km3
remove(discharge_1901_1930_daily)
high_flows_sum=rowSums(apply(high_flows,c(1,2),mean)) #from 30 years monthly to average yearly 
discharge_sum=rowSums(apply(discharge_1901_1930_monthly_full,c(1,2),mean)) #from 30 years monthly to average yearly 
plotGlobalW(data = high_flows_sum,title = "high flows yearly sum",file=paste0(fplotspath,"high_flows_yearly_sum_1901_1930mean.png"),pow2max = 10,pow2min=1,colPos = "Blues",legendtitle = "km3/yr",legYes = T,onlyPos = T,eps = F)
plotGlobalWlin(data = high_flows_sum/discharge_sum*100,title = "high_flows_sum/discharge_sum*100",colrev=T,file=paste0(fplotspath,"high_flows_yearly_perc_discharge_1901_1930mean_raindays.png"),min = 0,max=100,colPos = "Spectral",legendtitle = "%",legYes = T,onlyPos = T,eps = F,colNeg = "Reds")

discharge_1901_1930_monthly_noHF_full=aperm(apply(X = discharge_1901_1930_daily_noHF,MARGIN = c(1,3),FUN = daily2monthly,method="sum")/10^9,c(2,1,3)) # daily2monthly + mean 10y period + convert from m3/month to km3/month
remove(discharge_1901_1930_daily_noHF)
efrs_1901_1930_monthly_noHF=calcEFRs(discharge = discharge_1901_1930_monthly_noHF_full,method = "VMF")
discharge_1901_1930_monthly_noHF=apply(discharge_1901_1930_monthly_noHF_full,c(1,2),mean)
basinPBremote_woHighFlow = sumUpPBandEFRs(endcell = endcell,routing=routing,cellindex = cellindex,efr=rowSums(efrs_1901_1930_monthly_noHF),discharge = rowSums(discharge_1901_1930_monthly_noHF),remotes=fullRemote)
pbEFRstats(basinPBremote_woHighFlow,endcells)

#Q80
load(file="/media/All/PB_water_revised/data.RData")
discharge_1901_1930_daily_noHF80=removeHighFlows(discharge_1901_1930_daily,threshold=0.8)
remove(discharge_1901_1930_daily)
discharge_1901_1930_monthly_noHF80_full=aperm(apply(X = discharge_1901_1930_daily_noHF80,MARGIN = c(1,3),FUN = daily2monthly,method="sum")/10^9,c(2,1,3)) # daily2monthly + mean 10y period + convert from m3/month to km3/month
remove(discharge_1901_1930_daily_noHF80)
efrs_1901_1930_monthly_noHF80=calcEFRs(discharge = discharge_1901_1930_monthly_noHF80_full,method = "VMF")
discharge_1901_1930_monthly_noHF80=apply(discharge_1901_1930_monthly_noHF80_full,c(1,2),mean)
basinPBremote_woHighFlow80 = sumUpPBandEFRs(endcell = endcell,routing=routing,cellindex = cellindex,efr=rowSums(efrs_1901_1930_monthly_noHF80),discharge = rowSums(discharge_1901_1930_monthly_noHF80),remotes=fullRemote)#,removeHighFlow = T)
pbEFRstats(basinPBremote_woHighFlow80,endcells)


# =============== current PB status and overuse ==============
iFol="/p/projects/open/Fabian/runs/IRRHIST/cru_ts3.21/"
discharge_2000=apply(readMonthly(inFile = paste0(iFol,"historic_noefr_ilim_landuse/mdischarge.bin"),startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1996,headersize = 0,getyearstop = 2005),c(1,2),mean)*rep(nDays,each=67420)/10^3 # 10 year mean from hm3/d to km3/year
days=rep(nDays,each=67420,times=30)
dim(days)=c(67420,12,30)
discharge_1915=readMonthly(inFile = paste0(iFol,"historic_noefr_ilim_landuse/mdischarge.bin"),startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1901,headersize = 0,getyearstop = 1930)*days/10^3 # 10 year mean from hm3/d to km3/month
efrs_1915=calcEFRs(discharge = discharge_1915,method = "VMF") # in km3/yr
waterwd_2000=rowMeans(readYearly(inFile = paste0(iFol,"historic_noefr_ilim_landuse/awaterwithdrawal_hil.bin"),startyear = 1901,stopyear = 2005,size = 4,getyearstart = 1996,headersize = 0,getyearstop = 2005))/1000/10^9 # 10 year mean from l/yr to km3/yr
waterwd_2000=rep(waterwd_2000/12,times=12)
dim(waterwd_2000)=c(ncells,12)

irrigation_2000=apply(readMonthly(paste0(iFol,"historic_noefr_ilim_landuse/mirrig.bin"),1901,2005,4,0,1996,2005),c(1,2),mean)/1000*cellarea/10^9 # from mm/month to km3/yr
convlossdrain_2000=apply(readMonthly(paste0(iFol,"historic_noefr_ilim_landuse/mconv_loss_drain.bin"),1901,2005,4,0,1996,2005),c(1,2),mean)/1000*cellarea/10^9 # from mm/month to km3/yr
convlossevap_2000=apply(readMonthly(paste0(iFol,"historic_noefr_ilim_landuse/mconv_loss_evap.bin"),1901,2005,4,0,1996,2005),c(1,2),mean)/1000*cellarea/10^9 # from mm/month to km3/yr
irrig_wd_2000=irrigation_2000+convlossdrain_2000+convlossevap_2000

#monthly
overuse=efrs_1915+irrig_wd_2000+waterwd_2000-discharge_2000
overuse[overuse<0]=0
overuse=rowSums(overuse)
sum(irrig_wd_2000+waterwd_2000)
sum(overuse)
plotGlobalW(data = overuse,title = "",file=paste0(fplotspath,"overuse_PBwater2000_monthly.png"),pow2max = 4,pow2min=-7,colPos = "Reds",legendtitle = "km3/yr",legYes = T,onlyPos = T,eps = F)

#yearly
overuse_yearly=rowSums(efrs_1915)+rowSums(irrig_wd_2000)+rowSums(waterwd_2000)-rowSums(discharge_2000)
overuse_yearly[overuse_yearly<0]=0
sum(overuse_yearly)
plotGlobalW(data = overuse_yearly,title = "",file=paste0(fplotspath,"overuse_PBwater2000_yearly.png"),pow2max = 4,pow2min=-7,colPos = "Reds",legendtitle = "km3/yr",legYes = T,onlyPos = T,eps = F)

plotGlobalW(data=discharge_2000,file=paste0(fplotspath,"discharge2000.png"),title = "",pow2min = 0,pow2max = 12,legendtitle = "km3/yr",legYes = T,eps = F,onlyPos = T)

# ===== Jahresgang plot =====
load(file="/media/All/PB_water_revised/data.RData")
daily_display=discharge_1901_1930_daily[17048,,15]/10^9
quantiles= quantile(daily_display,probs=c(0.8,0.9))
#remove(discharge_1901_1930_daily)

expFormat="pdf"
  file=paste0(fplotspath,"daily_discharge_Q90.png")
  if (expFormat=="eps"){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=18, height=18,paper="special")
  }else if (expFormat=="pdf"){  
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"pdf"),collapse=".")
    pdf(file,width=6,height=3,family = c("Helvetica"),pointsize = 12,paper='special',version = "1.5")
  }else{
    png(file, width=3, height=2, units="in", res=400, pointsize=6,type="cairo")
  }
  par(mar=c(4.5,4.5,0.1,0.1),oma=c(0,0,0,0))
  plot(x=1:365,y=daily_display,type="h",lwd=0.3,lend=1,ylab="discharge [km3/day]",xlab="Day of the year")
  #axis(side = 1,at = seq(1,351,50),labels = seq(1,351,50))
  par(xpd=T)
  xp=5
  ce=1
  offs=0.003
  lines(x=c(1,365),y=c(quantiles[1],quantiles[1]),lty="dashed",col="gray40")
  text(x=xp,y=(quantiles[1]-offs),"Q80",col="gray40",cex=ce,adj=0)
  lines(x=c(1,365),y=c(quantiles[2],quantiles[2]),lty="dashed",col="gray40")
  text(x=xp,y=(quantiles[2]+offs),"Q90",col="gray40",cex=ce,adj=0)
  dev.off()
# ===== showNetwork =====
cells=17048#lon=-62.75,lat=10.25
for (c in cells){
  expFormat="pdf"
  file=paste0(fplotspath,"efr_",c,".png")
  if (expFormat=="eps"){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=18, height=18,paper="special")
  }else if (expFormat=="pdf"){  
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"pdf"),collapse=".")
    pdf(file,width=6,height=4,family = c("Helvetica"),pointsize = 12,paper='special',version = "1.5")
  }else{
    png(file, width=3, height=2, units="in", res=400, pointsize=6,type="cairo")
  }
  months=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  ms=c("J","F","M","A","M","J","J","A","S","O","N","D")
  co=c("lightblue",alpha(c("blue","darkgreen"),alpha = 0.4),"red")
  par(mar=c(4.5,4.5,0.1,0.1))#,oma=c(0,0,0,0))
  ma=max(c(discharge_1685[c,],discharge_2000[c,]))
  plot(NA, xlim=c(1,12),ylim=c(-1*max(efr_deficit[c,]),ma),xaxt="n",ylab="water [km3/month]",xlab="Month")
  axis(side = 1,at = 1:12,labels = ms)
  polygon(x=c(1,1:12,12),y=c(0,discharge_1685[c,],0),col=co[1],border = F)
  polygon(x=c(1,1:12,12),y=c(0,discharge_2000[c,],0),col=co[2],border = F)
  polygon(x=c(1,1:12,12),y=c(0,efr[c,],0),col=co[3],border = F)
  polygon(x=c(1,1:12,12),y=c(0,efr_deficit[c,]*-1,0),col=co[4],border = F)
  maf=mean(discharge_1685[c,])
  offs=(1/30*ma)
  xp=1.5
  ce=0.6
  lines(x=c(1,12),y=c(maf,maf),lty="dashed",col="gray40")
  text(x=xp,y=(maf+offs),"MAF",col="gray40",cex=ce,adj=0)
  lines(x=c(1,12.5),y=c(maf*0.8,maf*0.8),lty="dashed",col="gray40")
  text(x=xp,y=(maf*0.8+offs),"80% MAF",col="gray40",cex=ce,adj=0)
  lines(x=c(1,12.5),y=c(maf*0.4,maf*0.4),lty="dashed",col="gray40")
  text(x=xp,y=(maf*0.4+offs),"40% MAF",col="gray40",cex=ce,adj=0)
  text(x=rep(12.2,3),y=c(maf,0.6*maf,0.2*maf),labels = c("high flow","interm. flow","low flow"),col="gray40",cex=ce,srt = -90)
  legend("topleft",legend = c("pristine discharge (MMF)","recent discharge","EFRs","EFR deficit"),fill=co,cex=0.9,border=NA,bg = NA,bty = "n")
  dev.off()
}

plotBasin(cell=9544,endcellList = endcell,zoom = T,off=25) #6777
# lv and i are global variables to be filled by showNetwork function
lv=array(list(),dim=c(max(cellindex[which(endcell==endcell[cell])])))
i=1
cell=endcell[28025]
showNetwork(ind=cell,routingTable = routing)

if (basinExportCSV){
  lv_flat=unlist(lv)
  fileW=paste0(fplotspath,"efr_data_norway.csv")
  dat=data.frame(lv_flat,routing[lv_flat],round(efr[lv_flat],2),round(discharge_1685[lv_flat],2),round(discharge_1685[lv_flat]-efr[lv_flat],2),round(cellEFR[lv_flat],2),round(cellPB[lv_flat],2))
  colnames(dat)=c("ID","nextID","EFR","Discharge","purePB (Discharge-EFR)","calcEFR","calcPB")
  write.csv(x = dat,file=fileW)
}


# ===== Jahresgang plot =====
cells=17048#lon=-62.75,lat=10.25
for (c in cells){
  expFormat="pdf"
  file=paste0(fplotspath,"efr_",c,".png")
  if (expFormat=="eps"){
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width=18, height=18,paper="special")
  }else if (expFormat=="pdf"){  
    file=strsplit(file,".",fixed=TRUE)[[1]]
    file=paste(c(file[1:(length(file)-1)],"pdf"),collapse=".")
    pdf(file,width=6,height=4,family = c("Helvetica"),pointsize = 12,paper='special',version = "1.5")
  }else{
    png(file, width=3, height=2, units="in", res=400, pointsize=6,type="cairo")
  }
  months=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  ms=c("J","F","M","A","M","J","J","A","S","O","N","D")
  co=c("lightblue",alpha(c("blue","darkgreen"),alpha = 0.4),"red")
  par(mar=c(4.5,4.5,0.1,0.1))#,oma=c(0,0,0,0))
  ma=max(c(discharge_1685[c,],discharge_2000[c,]))
  plot(NA, xlim=c(1,12),ylim=c(-1*max(efr_deficit[c,]),ma),xaxt="n",ylab="water [km3/month]",xlab="Month")
  axis(side = 1,at = 1:12,labels = ms)
  polygon(x=c(1,1:12,12),y=c(0,discharge_1685[c,],0),col=co[1],border = F)
  polygon(x=c(1,1:12,12),y=c(0,discharge_2000[c,],0),col=co[2],border = F)
  polygon(x=c(1,1:12,12),y=c(0,efr[c,],0),col=co[3],border = F)
  polygon(x=c(1,1:12,12),y=c(0,efr_deficit[c,]*-1,0),col=co[4],border = F)
  maf=mean(discharge_1685[c,])
  offs=(1/30*ma)
  xp=1.5
  ce=0.6
  lines(x=c(1,12),y=c(maf,maf),lty="dashed",col="gray40")
  text(x=xp,y=(maf+offs),"MAF",col="gray40",cex=ce,adj=0)
  lines(x=c(1,12.5),y=c(maf*0.8,maf*0.8),lty="dashed",col="gray40")
  text(x=xp,y=(maf*0.8+offs),"80% MAF",col="gray40",cex=ce,adj=0)
  lines(x=c(1,12.5),y=c(maf*0.4,maf*0.4),lty="dashed",col="gray40")
  text(x=xp,y=(maf*0.4+offs),"40% MAF",col="gray40",cex=ce,adj=0)
  text(x=rep(12.2,3),y=c(maf,0.6*maf,0.2*maf),labels = c("high flow","interm. flow","low flow"),col="gray40",cex=ce,srt = -90)
  legend("topleft",legend = c("pristine discharge (MMF)","recent discharge","EFRs","EFR deficit"),fill=co,cex=0.9,border=NA,bg = NA,bty = "n")
  dev.off()
}
# ===== discharge cells algorithm =====
expFormat="png"
file=paste0(fplotspath,"cell_efr_pb.png")
if (expFormat=="eps"){
  file=strsplit(file,".",fixed=TRUE)[[1]]
  file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
  ps.options(family = c("Helvetica"), pointsize = 18)
  postscript(file,horizontal = FALSE, onefile = FALSE, width=18, height=18,paper="special")
}else if (expFormat=="pdf"){  
  file=strsplit(file,".",fixed=TRUE)[[1]]
  file=paste(c(file[1:(length(file)-1)],"pdf"),collapse=".")
  pdf(file,width=4,height=6,family = c("Helvetica"),pointsize = 12,paper='special',version = "1.5")
}else{
  png(file, width=3, height=3, units="in", res=400, pointsize=6,type="cairo")
}
colz=c("white","lightblue","lightgreen")
par(mar=c(0.5,2,0.1,0.1),oma=c(0,0,0,0),fig=c(0,0.9,0,1))
wd=0.5
dataD=matrix(c(0,0,1.4,1.4,2.3,0.45,1.6,2.7,0.55),nrow = 3,byrow = T )
a=barplot(height = dataD,ylim = c(0,5),col=colz,space=2,border = colz,xlim=c(0,5),width=wd,cex.names=0.8,names.arg = c("","",""),ylab="",axes=F,yaxs="r")
arrows(x0 = c(a[1]+wd,a[3]-wd),x1 = c(a[2]-wd,a[2]+wd),y0 = c(1.4,1.8),y1 = c(2.3,2.3),code=2,length = 0.15)
arrows(x0=0.5,x1=0.5,y0=0,y1=3,angle = 90,code = 3,length = 0.07)
text(x = 0.2,y=1.5,"total discharge",srt=90)
legend(0,5,legend = c("EFR","available discharge"),fill = colz[2:3],border = colz[2:3])
par(mar=c(0.5,0,0.1,0.1),oma=c(0,0,0,0),fig=c(0.9,1,0,1),new=T)
plot(NA,xlim=c(0,1),ylim=c(0,5),axes=F,yaxs="r")
arrows(x0=c(0.3),x1=c(0.3),y0=c(0),y1=c(1.85),angle = 90,code = 3,length = 0.07,col="red")
text(x = 0.7,y=1,c("upstream EFR sum"),srt=90)
dev.off()

# ===== wateruse scenarios isimip2a ####
require(ncdf4) # package for netcdf manipulation
types=c("airrww","adomww","aliveww","aindww")
ce=1.5

expFormat="pdf"
file=paste0(fplotspath,"future_wu.png")
if (expFormat=="eps"){
  file=strsplit(file,".",fixed=TRUE)[[1]]
  file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
  ps.options(family = c("Helvetica"), pointsize = 18)
  postscript(file,horizontal = FALSE, onefile = FALSE, width=18, height=18,paper="special")
}else if (expFormat=="pdf"){  
  file=strsplit(file,".",fixed=TRUE)[[1]]
  file=paste(c(file[1:(length(file)-1)],"pdf"),collapse=".")
  pdf(file,width=12,height=6,family = c("Helvetica"),pointsize = 12,paper='special',version = "1.5")
}else{
  png(file, width=6, height=3, units="in", res=400, pointsize=6,type="cairo")
}
par(mar=c(5,5,2,2))
plot(0,xlim=c(2006,2099),ylim=c(0,5000),cex.axis=ce,cex.lab=ce,xlab="Year",ylab=expression(paste('Freshwater withdrawals [',km^3,'/yr]',sep='')),xaxs="i")
colz=RColorBrewer::brewer.pal(8,"Set2")[c(3,2,5,6)]
for (i in 1:4){
  files=list.files(path="/p/projects/open/Fabian/dissData/", pattern=paste0("*",types[i],"*"))
  dat=array(0,dim=c(94,length(files)))
  for (k in 1:length(files)){
    nc_data=nc_open(file=paste0("/p/projects/open/Fabian/dissData/",files[k]))
    dat[,k]=rowSums(matrix(ncvar_get(nc_data, types[i])/10^9/1000*30*24*3600,ncol = 12,byrow = T)) # from kg/s to km3/month
  }
  polygon(x=c(seq(2006,2099,1),seq(2099,2006,-1)),y=c(apply(dat,c(1),max),rev(apply(dat,c(1),min))),col=alpha(colz[i],alpha=0.4),border=NA)
  lines(x=seq(2006,2099,1),y=apply(dat,c(1),mean),col=alpha(colz[i],alpha=0.8))
}
legend(2006,4100,legend = c("Irrigation","Domestic","Livestock","Industry"),fill= colz,cex=1.5,border = colz,bg = "white")
dev.off()


# ===== wateruse scenarios isimip2a 2005soc aXXXww ####
require(ncdf4) # package for netcdf manipulation
types=c("airrww","adomww","aliveww","aindww")
ce=1.5

expFormat="pdf"
file=paste0(fplotspath,"future_2005soc_aXXXww_ISIMIP2b_water_global.png")
if (expFormat=="eps"){
  file=strsplit(file,".",fixed=TRUE)[[1]]
  file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
  ps.options(family = c("Helvetica"), pointsize = 18)
  postscript(file,horizontal = FALSE, onefile = FALSE, width=18, height=18,paper="special")
}else if (expFormat=="pdf"){  
  file=strsplit(file,".",fixed=TRUE)[[1]]
  file=paste(c(file[1:(length(file)-1)],"pdf"),collapse=".")
  pdf(file,width=12,height=6,family = c("Helvetica"),pointsize = 12,paper='special',version = "1.5")
}else{
  png(file, width=6, height=3, units="in", res=400, pointsize=6,type="cairo")
}
par(mar=c(5,5,2,2))
plot(0,xlim=c(2006,2099),ylim=c(0,5000),cex.axis=ce,cex.lab=ce,xlab="Year",ylab=expression(paste('Freshwater withdrawals [',km^3,'/yr]',sep='')),xaxs="i")
colz=RColorBrewer::brewer.pal(8,"Set2")[c(3,2,5,6,8,4)]
for (i in 1:4){
  files=list.files(path="/p/projects/open/Fabian/dissData/act/2005soc/", pattern=paste0("*",types[i],"*"))
  dat=array(0,dim=c(94,length(files)))
  for (k in 1:length(files)){
    nc_data=nc_open(file=paste0("/p/projects/open/Fabian/dissData/act/2005soc/",files[k]))
    dat[,k]=rowSums(matrix(ncvar_get(nc_data, types[i])/10^9/1000*30*24*3600,ncol = 12,byrow = T)) # from kg/s to km3/month
    lines(x=seq(2006,2099,1),y=dat[,k],col=alpha(colz[i],alpha=0.8))
    #text(x=2040,y=mean(dat[,k]),labels=files[k],cex=0.2)
  }
  #polygon(x=c(seq(2006,2099,1),seq(2099,2006,-1)),y=c(apply(dat,c(1),max),rev(apply(dat,c(1),min))),col=alpha(colz[i],alpha=0.4),border=NA)
  #lines(x=seq(2006,2099,1),y=apply(dat,c(1),mean),col=alpha(colz[i],alpha=0.8))
  if (F){ #i==2
    edat=colMeans(dat)
    dim(edat)=c(length(edat),1)
    rownames(edat)=files
  }
}

text(x=2040,y=c(4100,2200,1700,1400,750,300),labels=c("matsiro","lpjml","cwatm","h08","pcr-globwb","mpi-hm"),cex=1)
text(x=2020,y=c(600,450,50,1050),labels=c("pcr-globwb","h08","pcr-globwb","pcr-globwb"),cex=1)

legend(2006,3700,legend = c("Irrigation","Domestic","Livestock","Industry",""),fill= colz,cex=1.0,border = colz,bg = "white")
dev.off()

# ===== wateruse scenarios isimip2a rcpXXsoc pXXXuse ####
require(ncdf4) # package for netcdf manipulation
types=c("pdomuse","pmanuse")
ce=1.5

expFormat="pdf"
file=paste0(fplotspath,"future_rcpXXsoc_pXXXuse_ISIMIP2b_water_global.png")
if (expFormat=="eps"){
  file=strsplit(file,".",fixed=TRUE)[[1]]
  file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
  ps.options(family = c("Helvetica"), pointsize = 18)
  postscript(file,horizontal = FALSE, onefile = FALSE, width=18, height=18,paper="special")
}else if (expFormat=="pdf"){  
  file=strsplit(file,".",fixed=TRUE)[[1]]
  file=paste(c(file[1:(length(file)-1)],"pdf"),collapse=".")
  pdf(file,width=12,height=6,family = c("Helvetica"),pointsize = 12,paper='special',version = "1.5")
}else{
  png(file, width=6, height=3, units="in", res=400, pointsize=6,type="cairo")
}
par(mar=c(5,5,2,2))
plot(0,xlim=c(2006,2099),ylim=c(0,200),cex.axis=ce,cex.lab=ce,xlab="Year",ylab=expression(paste('Freshwater consumption [',km^3,'/yr]',sep='')),xaxs="i")
colz=RColorBrewer::brewer.pal(8,"Set2")[c(8,4)]
for (i in 1:2){
  files=list.files(path="/p/projects/open/Fabian/dissData/pot/", pattern=paste0("*",types[i],"*"))
  dat=array(0,dim=c(94,length(files)))
  for (k in 1:length(files)){
    nc_data=nc_open(file=paste0("/p/projects/open/Fabian/dissData/pot/",files[k]))
    dat[,k]=rowSums(matrix(ncvar_get(nc_data, types[i])/10^9/1000*30*24*3600,ncol = 12,byrow = T)) # from kg/s to km3/month
    lines(x=seq(2006,2099,1),y=dat[,k],col=alpha(colz[i],alpha=0.8))
    #text(x=2040,y=mean(dat[,k]),labels=files[k],cex=0.5)
  }
  if (F){ #i==2
    edat=colMeans(dat)
    dim(edat)=c(length(edat),1)
    rownames(edat)=files
  }
}

legend(2006,200,legend = c("pdomuse","pmanuse"),fill= colz,cex=1.0,border = colz,bg = "white")
dev.off()
# ===== wateruse scenarios isimip2a rcpXXsoc aXXXww ####
require(ncdf4) # package for netcdf manipulation
types=c("airrww","adomww","amanww")
ce=1.5

expFormat="pdf"
file=paste0(fplotspath,"future_rcpXXsoc_aXXXww_ISIMIP2b_water_global.png")
if (expFormat=="eps"){
  file=strsplit(file,".",fixed=TRUE)[[1]]
  file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
  ps.options(family = c("Helvetica"), pointsize = 18)
  postscript(file,horizontal = FALSE, onefile = FALSE, width=18, height=18,paper="special")
}else if (expFormat=="pdf"){  
  file=strsplit(file,".",fixed=TRUE)[[1]]
  file=paste(c(file[1:(length(file)-1)],"pdf"),collapse=".")
  pdf(file,width=12,height=6,family = c("Helvetica"),pointsize = 12,paper='special',version = "1.5")
}else{
  png(file, width=6, height=3, units="in", res=400, pointsize=6,type="cairo")
}
par(mar=c(5,5,2,2))
plot(0,xlim=c(2006,2099),ylim=c(0,8000),cex.axis=ce,cex.lab=ce,xlab="Year",ylab=expression(paste('Freshwater withdrawals [',km^3,'/yr]',sep='')),xaxs="i")
colz=RColorBrewer::brewer.pal(8,"Set2")[c(3,4,2)]
for (i in 1:3){
  files=list.files(path="/p/projects/open/Fabian/dissData/act/rcpXXsoc_aXXXww/", pattern=paste0("*",types[i],"*"))
  dat=array(0,dim=c(94,length(files)))
  for (k in 1:length(files)){
    nc_data=nc_open(file=paste0("/p/projects/open/Fabian/dissData/act/rcpXXsoc_aXXXww/",files[k]))
    dat[,k]=rowSums(matrix(ncvar_get(nc_data, types[i])/10^9/1000*30*24*3600,ncol = 12,byrow = T)) # from kg/s to km3/month
    lines(x=seq(2006,2099,1),y=dat[,k],col=alpha(colz[i],alpha=0.8))
    #text(x=2040,y=mean(dat[,k]),labels=files[k],cex=0.2)
  }
  #polygon(x=c(seq(2006,2099,1),seq(2099,2006,-1)),y=c(apply(dat,c(1),max),rev(apply(dat,c(1),min))),col=alpha(colz[i],alpha=0.4),border=NA)
  #lines(x=seq(2006,2099,1),y=apply(dat,c(1),mean),col=alpha(colz[i],alpha=0.8))
  if (F){ #i==2
    edat=colMeans(dat)
    dim(edat)=c(length(edat),1)
    rownames(edat)=files
  }
}

#text(x=2040,y=c(4100,2200,1700,1400,750,300),labels=c("matsiro","lpjml","cwatm","h08","pcr-globwb","mpi-hm"),cex=1)
#text(x=2020,y=c(600,450,50,1050),labels=c("pcr-globwb","h08","pcr-globwb","pcr-globwb"),cex=1)

legend(2006,8000,legend = c("Irrigation","Domestic","Manufacturing/Industry"),fill= colz,cex=1.0,border = colz,bg = "white")
dev.off()


# ===== reanalysis ISIMIP2b data - 1st gen BE ####

expFormat="pdf"
file=paste0(fplotspath,"future_rcpXXsoc_aXXXww_ISIMIP2b_water_global.png")
if (expFormat=="eps"){
  file=strsplit(file,".",fixed=TRUE)[[1]]
  file=paste(c(file[1:(length(file)-1)],"eps"),collapse=".")
  ps.options(family = c("Helvetica"), pointsize = 18)
  postscript(file,horizontal = FALSE, onefile = FALSE, width=18, height=18,paper="special")
}else if (expFormat=="pdf"){  
  file=strsplit(file,".",fixed=TRUE)[[1]]
  file=paste(c(file[1:(length(file)-1)],"pdf"),collapse=".")
  pdf(file,width=12,height=6,family = c("Helvetica"),pointsize = 12,paper='special',version = "1.5")
}else{
  png(file, width=6, height=3, units="in", res=400, pointsize=6,type="cairo")
}
par(mar=c(5,5,2,2))
plot(0,xlim=c(2006,2099),ylim=c(0,8000),cex.axis=ce,cex.lab=ce,xlab="Year",ylab=expression(paste('Freshwater withdrawals [',km^3,'/yr]',sep='')),xaxs="i")
colz=RColorBrewer::brewer.pal(8,"Set2")[c(3,4,2)]
for (i in 1:3){
  files=list.files(path="/p/projects/open/Fabian/dissData/act/rcpXXsoc_aXXXww/", pattern=paste0("*",types[i],"*"))
  dat=array(0,dim=c(94,length(files)))
  for (k in 1:length(files)){
    nc_data=nc_open(file=paste0("/p/projects/open/Fabian/dissData/act/rcpXXsoc_aXXXww/",files[k]))
    dat[,k]=rowSums(matrix(ncvar_get(nc_data, types[i])/10^9/1000*30*24*3600,ncol = 12,byrow = T)) # from kg/s to km3/month
    lines(x=seq(2006,2099,1),y=dat[,k],col=alpha(colz[i],alpha=0.8))
    #text(x=2040,y=mean(dat[,k]),labels=files[k],cex=0.2)
  }
  #polygon(x=c(seq(2006,2099,1),seq(2099,2006,-1)),y=c(apply(dat,c(1),max),rev(apply(dat,c(1),min))),col=alpha(colz[i],alpha=0.4),border=NA)
  #lines(x=seq(2006,2099,1),y=apply(dat,c(1),mean),col=alpha(colz[i],alpha=0.8))
  if (F){ #i==2
    edat=colMeans(dat)
    dim(edat)=c(length(edat),1)
    rownames(edat)=files
  }
}

#text(x=2040,y=c(4100,2200,1700,1400,750,300),labels=c("matsiro","lpjml","cwatm","h08","pcr-globwb","mpi-hm"),cex=1)
#text(x=2020,y=c(600,450,50,1050),labels=c("pcr-globwb","h08","pcr-globwb","pcr-globwb"),cex=1)

legend(2006,8000,legend = c("Irrigation","Domestic","Manufacturing/Industry"),fill = colz,cex = 1.0,border = colz,bg = "white")
dev.off()





################## quick check for PB 3.0 ###################
iFol <- "/p/projects/open/Fabian/runs/IRRHIST/cru_ts3.21/"
discharge_PI <- readMonthly(inFile = paste0(iFol,"efrcalc_efr_ilim_nolu/mdischarge.bin"),startyear = 1670,stopyear = 1699,
                         size = 4,getyearstart = 1670,headersize = 0,getyearstop = 1699)*rep(nDays,each = 67420)/10^3 #from hm3/d to km3/month
quantiles <- apply(discharge_PI,c(1,2),quantile,probs = c(0.5,0.9,0.05,0.95))
