#' Aggregate discharge composition as EFRs and PBs from grid cell to global
#'
#' Aggregate discharge composition as EFRs and PBs from grid cell to global
#'
#' @param endcell array containing the index of each cells final drainage cell
#'
#' @param routing array containing the index of each cells next drainage cell
#'
#' @param cellindex array containing the rank of each cell in the drainage
#'        network (1: final drainage cell, 2: draining directly into final
#'        drainage cell, ...)
#'
#' @param efr list
#'
#' @param discharge list
#'
#' @param remotes list
#'
#' @param flowseason list
#'
#' @return global array of aggregated EFR and PB
#'
#' @examples
#' \dontrun{
#' }
efr_aggregation <- function(endcell,
                            routing,
                            cellindex,
                            efr,
                            discharge,
                            remotes=F,
                            flowseason=F){
  basins=which(routing==0) # indices of ocean draining cells (due to adding 1, -1 becomes 0)
  totalEFR=routing*0
  totalPB=routing*0
  a=1:length(basins)
  if (length(remotes>1)){
    efr[remotes==2]=discharge[remotes==2] #outflowcells of remote areas: remotes==2
    efr[remotes==1]=discharge[remotes==1] #remote areas: remotes==1
  }
  #if (length(flowseason>1)){
  #  #fs=flowSeason(discharge=discharge,method = "VMF")
  #  #efr[fs=="h"]=discharge[fs=="h"] #outflowcells of remote areas: remotes==2
  #}
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

upstreamCells <- function(ind,routingTable){
  return(which(routingTable==ind))
}

#' @param endcell array containing the index of each cells final drainage cell
#' @param spatial_scale aggregate to "global" or "basin"
pbEFRstats <- function(efrPBarray, endcell, spatial_scale = "global"){
  final_drainage_cell_indices <- sort(unique(endcell))
  if (spatial_scale == "global"){
    return(colSums(efrPBarray[final_drainage_cell_indices,]))
  }else if (spatial_scale == "basin"){
    return(efrPBarray[final_drainage_cell_indices,])
  }else{
    stop(paste("Unknown setting '",spatial_scale,"' for parameter spatial_scale. Use 'global' or 'basin'"))
  }
}

#' Apply data for each basin to all cells of that basin
#'
#' Apply data for each basin to all cells of that basin
#'
#' @param data array with information for each basin
#' @param endcell lpjml array giving the index of the final drainage cell
#'                for each cell
#'
#' @return array with data for each grid cell, based on the final drainage cell
#'
#' @examples
#'
#' @export
applyBasinData <- function(data, endcell){
  ncells <- length(endcell)
  basins <- sort(unique(endcell))
  out <- array(NA,dim=ncells)
  for (i in 1:length(basins)){
    basin <- basins[i]
    out[which(endcell==basin)] <- data[i]
  }
  return(out)
}
