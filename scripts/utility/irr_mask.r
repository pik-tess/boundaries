
ncells = 67420
showRoute <- function(ind,routingTable){
  if (routingTable[ind]<1){return(ind)} #can be 0 or -8 -> endcell or nacell
  else{return(c(ind,showRoute(routingTable[ind],routingTable)))}
}
drainage = autoReadInput(
    inFile = "/p/projects/lpjml/input/historical/input_VERSION2/drainage.bin",
    manu = TRUE,
    msize = 4,
    mheadersize = 43
)
routing=drainage[1,]+1 #add 1 since in C indexing starts at 0 but in R at 1
endcell_indices=which(routing==0) # indices of ocean draining cells (due to adding 1, -1 becomes 0)
nacells=which(routing==-8) # cells that were not part of STN, LPJmL treats all values below 0 as outflow cells
endcell=array(0,dim=ncells)
cellindex = endcell
for (c in 1:ncells){
  #write(c,stdout())
  route=showRoute(c,routing)
  cellindex[route]=seq(length(route),1,-1)
  endcell[c]=route[length(route)]
}

  basin.ids <- sort(unique(endcell))
  irrmask_basin=array(0,ncells)
  for(i in 1:length(basin.ids)){
    basincells <- which(endcell==basin.ids[i])
    if(sum(air[basincells])>0) {
      irrmask_basin[basincells]=1
    }
  }

repl <- rep(0,67420)
repl[endcell] <- 1


test2 <- LPJmLTools::lpjmlinfo$cells_raster

test[LPJmLTools::lpjmlinfo$cellnumbers] <- repl

test2[LPJmLTools::lpjmlinfo$cellnumbers] <- irrmask_basin
plot(LPJmLTools::lpjmlinfo$cells_raster)