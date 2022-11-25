# written by Fabian Stenzel, based on work by Sebastian Ostberg
# 2022 - stenzel@pik-potsdam.de

#requirements:
require(lpjmliotools)

################# gamma calc functions  ###################

#' Calculate the ecosystem change metric gamma between 2 simulations/timesteps
#'
#' Function to calculate the ecosystem change metric gamma according
#' to Sykes (1999), Heyder (2011), and Ostberg (2015,2018).
#' This is a reformulated version in R, not producing 100% similar values
#' than the C/bash version from Ostberg 2018, but following their methodology.
#'
#' @param folderRef folder of reference run
#' @param folderRef2 2nd folder of reference run, if split (default NULL -- no split)
#' @param folderScen folder of scenario run
#' @param saveFileData file to save read in data to (default NULL)
#' @param saveFileGamma file to save gamma data to (default NULL)
#' @param exportGamma flag whether to export gamma data to global environment (default F)
#' @param combined read in soilc+litc, evap+interc, rh+harvestc separate or as combined outputs (default F -- separate)
#' @param nitrogen include nitrogen outputs for pools and fluxes into gamma calculation (default F)
#' @param weighting apply "old" (Ostberg-like), "new", or "equal" weighting of deltaV weights (default "old")
#' @param varnames data.frame with names of output files -- can be specified to account for variable file names (default NULL -- standard names are used)
#' @param headerout headersize of the output files (default 0)
#' @param timespan_full_ref full timespan of output files in reference period
#' @param timespan_full_scen full timespan of output files in scenario period
#' @param timespan_focus_ref specific years to use as reference period
#' @param timespan_focus_ref2 specific years to use as 2nd part of reference period in case this is split (default NULL)
#' @param timespan_focus_scen specific years to use as scenario period
#' @param npfts number of natural plant functional types (== bands in fpc - 1)
#' @param ncfts number of crop functional types
#' @param nbfts number of bfts (bioenergy functional types)
#' @param ncells number of cells in lpjml grid
#'
#' @return list data object containing arrays of deltaV, local_change,
#'         global_change, ecosystem_balance, carbon_fluxes, carbon_stocks,
#'         water_fluxes, gamma_total (+ nitrogen_fluxes and nitrogen_stocks if used)
#'
#' @examples
#' \dontrun{
#' gamma_today <- calcGamma(folderRef = "/p/projects/elis/ostberg/biome_shift/archive_Ostberg2018/simulations/scenario_PNV_trunkupdate_noCC/HadGEM2-ES/r1i1p1/",
#'                folderScen = "/p/projects/elis/ostberg/biome_shift/archive_Ostberg2018/simulations/scenario_LUC_trunkupdate/rcp26/HadGEM2-ES/r1i1p1/luc_rcp26_50bg50btflexi/",
#'                readPreviouslySavedData = TRUE, saveFileData = "/media/All/Gamma/Gamma_today_data.RData",
#'                saveFileGamma = "/media/All/Gamma/Gamma_today_gamma.RData",
#'                timespan_full_ref = c(1991,2099), timespan_full_scen = c(1991,2099),
#'                timespan_focus_ref = c(2000,2029), timespan_focus_scen = c(2000,2029),
#'                npfts = 9, ncfts = 32, nbfts = 12, ncells = 64240,
#'                combined = TRUE, weighting = "old")
#' }
#' @export
calcGamma <- function(folderRef, folderRef2 = NULL, folderScen, readPreviouslySavedData = TRUE,
                      saveFileData = NULL, saveFileGamma = NULL, combined = FALSE, nitrogen = FALSE,
                      weighting = "old", varnames = NULL, headerout = 0, timespan_full_ref,
                      timespan_full_scen, timespan_focus_ref, timespan_focus_ref2 = NULL,
                      timespan_focus_scen, npfts, ncfts, nbfts, ncells,dimensionsOnlyLocal=T) {
   require(lpjmliotools)

   if (is.null(varnames)) {
     print("variable name list not provided, using standard list, which might not be applicable for this case ...")
     varnames <- data.frame(row.names = c("grid",    "fpc",    "fpc_bft",    "cftfrac",    "firec",    "rh_harvest",    "npp",    "evapinterc",    "runoff",    "transp",    "soillitc",    "vegc",    "swcsum",    "firef",    "rh",    "harvestc",    "evap",    "interc",    "soilc",    "litc",    "swc",    "vegn",    "soilnh4",    "soilno3",    "leaching",    "n2o_denit",    "n2o_nit",    "n2o_denit",    "n2_emis",    "bnf",    "n_volatilization"),
                              outname = c("grid.bin","fpc.bin","fpc_bft.bin","cftfrac.bin","firec.bin","rh_harvest.bin","npp.bin","evapinterc.bin","runoff.bin","transp.bin","soillitc.bin","vegc.bin","swcsum.bin","firef.bin","rh.bin","harvestc.bin","evap.bin","interc.bin","soilc.bin","litc.bin","swc.bin","vegn.bin","soilnh4.bin","soilno3.bin","leaching.bin","n2o_denit.bin","n2o_nit.bin","n2o_denit.bin","n2_emis.bin","bnf.bin","n_volatilization.bin"),
                             timestep = c("Y",       "Y",      "Y",          "Y",          "Y",        "Y",             "Y",      "Y",             "Y",         "Y",         "Y",           "Y",       "Y",         "Y",        "Y",     "Y",           "Y",       "Y",         "Y",        "Y",       "Y",      ,"Y",      "Y",          "Y",          "Y",           "Y",            "Y",          "Y",            "Y",          "Y",      "Y"))
   }
   if (is.null(folderRef2)) {
      regular <-  TRUE
      nyears <- timespan_focus_ref[2] - timespan_focus_ref[1] + 1
   }else{
      regular <-  FALSE
      nyears1 <- timespan_focus_ref[2] - timespan_focus_ref[1] + 1
      nyears2 <- timespan_focus_ref2[2] - timespan_focus_ref2[1] + 1
      nyears <- nyears1 + nyears2
   }

   if (readPreviouslySavedData) {
     if (!is.null(saveFileData)) {
       print(paste("Loading saved data from:",saveFileData))
       load(file = saveFileData)
     }else{
       stop("saveFileData is not specified as parameter, nothing to load ... exiting")
     }

   }else{
     #first read in all lpjml output files required for computing gamma from reference and scenario folder
     returned_vars <- readGammaData(folderRef = folderRef, folderRef2 = folderRef2, folderScen = folderScen, saveFile = saveFileData, export = F,
                   combined = combined, nitrogen = nitrogen, varnames = varnames, headerout = headerout,
                   timespan_full_ref = timespan_full_ref, timespan_full_scen = timespan_full_scen,
                   timespan_focus_ref = timespan_focus_ref, timespan_focus_ref2 = timespan_focus_ref2, timespan_focus_scen = timespan_focus_scen,
                   npfts = npfts, ncfts = ncfts, nbfts = nbfts, ncells = ncells)
     state_ref <- returned_vars$state_ref
     mean_state_ref <- returned_vars$mean_state_ref
     state_scen <- returned_vars$state_scen
     mean_state_scen <- returned_vars$mean_state_scen
     fpc_ref <- returned_vars$fpc_ref
     fpc_scen <- returned_vars$fpc_scen
     bft_ref <- returned_vars$bft_ref
     bft_scen <- returned_vars$bft_scen
     cft_ref <- returned_vars$cft_ref
     cft_scen <- returned_vars$cft_scen
     lat <- returned_vars$lat
     lon <- returned_vars$lon
     cellArea <- returned_vars$cellArea
     rm(returned_vars)
   }


   ############## calc deltaV and variability of deltaV within reference period S(deltaV,sigma_deltaV) ##########
   fpc_ref_mean <- apply(fpc_ref, c(1,2), mean)
   bft_ref_mean <- apply(bft_ref, c(1,2), mean)
   cft_ref_mean <- apply(cft_ref, c(1,2), mean)
   sigma_deltaV_ref_list <- array(0, dim = c(ncells, nyears))
   # calculate for every year of the reference period, deltaV between that year and the average reference period year
   # this gives the variability of deltaV within the reference period
   for (y in 1:nyears) {
      sigma_deltaV_ref_list[,y] <- calcDeltaV(fpcRef = fpc_ref_mean,fpcScen = fpc_ref[,,y],
                                              bftRef = bft_ref_mean,bftScen = bft_ref[,,y],cftRef = cft_ref_mean,
                                              cftScen = cft_ref[,,y], weighting = weighting)
   }
   deltaVsd <- apply(sigma_deltaV_ref_list,c(1),sd)# calculate the standard deviation over the reference period for each gridcell
   # calculate deltaV between average reference and average scenario period
   deltaV <- calcDeltaV(fpcRef = fpc_ref_mean, fpcScen = apply(fpc_scen,c(1,2),mean),
                        bftRef = bft_ref_mean, bftScen = apply(bft_scen,c(1,2),mean), cftRef = cft_ref_mean,
                        cftScen = apply(cft_scen,c(1,2),mean), weighting = weighting)
   #
   ####

   ############## calc gamma components ################
   # variable names for the state vector
   # 1:3 carbon fluxes
   # 4:6 water fluxes
   # 7:8 carbon pools/stocks,
   # 9:10 additional variables for global/local difference, but not included in stocks/fluxes
   # 11:12 nitrogen pools/stocks
   # 13:15 nitrogen fluxes
   #                 1         2         3         4         5        6          7        8      9       10      11      12       13      14        15
   var_names <- c("firec","rh_harvest","npp","evapinterc","runoff","transp","soillitc","vegc","swcsum","firef","soiln","vegn","leaching","bnf","aggregated_n_emissions")

   delta <- deltaV*S_change_to_var_ratio(deltaV, deltaVsd)# deltaV
   lc <- calcComponent(ref = state_ref, scen = state_scen, local = TRUE, cellArea = cellArea)         # local change
   gc <- calcComponent(ref = state_ref, scen = state_scen, local = FALSE, cellArea = cellArea)         # global importance
   eb <- calcEcosystemBalance(ref = state_ref, scen = state_scen)                                # ecosystem balance
   if (dimensionsOnlyLocal == T){
     cf <- calcComponent(ref = state_ref[,,1:3], scen = state_scen[,,1:3], local = TRUE, cellArea = cellArea)     # carbon fluxes (local change)
     cs <- calcComponent(ref = state_ref[,,7:8], scen = state_scen[,,7:8], local = TRUE, cellArea = cellArea)     # carbon stocks (local change)
     wf <- calcComponent(ref = state_ref[,,4:6], scen = state_scen[,,4:6], local = TRUE, cellArea = cellArea)     # water fluxes (local change)
     if (nitrogen) {
       ns <- calcComponent(ref = state_ref[,,11:12], scen = state_scen[,,11:12], local = TRUE, cellArea = cellArea) # nitrogen stocks (local change)
       nf <- calcComponent(ref = state_ref[,,13:15], scen = state_scen[,,13:15], local = TRUE, cellArea = cellArea) # nitrogen fluxes (local change)
     }
   }else{
     cf <- (   calcComponent(ref = state_ref[,,1:3], scen = state_scen[,,1:3], local = TRUE, cellArea = cellArea)     # carbon fluxes (local change)
               + calcComponent(ref = state_ref[,,1:3], scen = state_scen[,,1:3], local = FALSE, cellArea = cellArea)
               + calcEcosystemBalance(ref = state_ref[,,1:3], scen = state_scen[,,1:3]) )/3
     cs <- (   calcComponent(ref = state_ref[,,7:8], scen = state_scen[,,7:8], local = TRUE, cellArea = cellArea)     # carbon stocks (local change)
               + calcComponent(ref = state_ref[,,7:8], scen = state_scen[,,7:8], local = FALSE, cellArea = cellArea)
               + calcEcosystemBalance(ref = state_ref[,,7:8], scen = state_scen[,,7:8]) )/3
     wf <- (   calcComponent(ref = state_ref[,,4:6], scen = state_scen[,,4:6], local = TRUE, cellArea = cellArea)     # water fluxes (local change)
               + calcComponent(ref = state_ref[,,4:6], scen = state_scen[,,4:6], local = FALSE, cellArea = cellArea)
               + calcEcosystemBalance(ref = state_ref[,,4:6], scen = state_scen[,,4:6]) )/3
     if (nitrogen) {
       ns <- (   calcComponent(ref = state_ref[,,11:12], scen = state_scen[,,11:12], local = TRUE, cellArea = cellArea)  # nitrogen stocks (local change)
                 + calcComponent(ref = state_ref[,,11:12], scen = state_scen[,,11:12], local = FALSE, cellArea = cellArea)
                 + calcEcosystemBalance(ref = state_ref[,,11:12], scen = state_scen[,,11:12]) )/3
       nf <- (   calcComponent(ref = state_ref[,,13:15], scen = state_scen[,,13:15], local = TRUE, cellArea = cellArea) # nitrogen fluxes (local change)
                 + calcComponent(ref = state_ref[,,13:15], scen = state_scen[,,13:15], local = FALSE, cellArea = cellArea)
                 + calcEcosystemBalance(ref = state_ref[,,13:15], scen = state_scen[,,13:15]) )/3
    }
   }

   # calc total gamma as the average of the 4 components
   gammaFull <- (delta + lc + gc + eb)/4 #check for NAs
   #
   ### FS: 2022-11-09 removed
   #gamma = expression(gamma, deltav, local, global, balance,
   #                  cstocks, cfluxes, wfluxes, nstocks, nfluxes)
   if (nitrogen) {
     gamma <- list(gamma_total = gammaFull, deltaV = delta, local_change = lc, global_change = gc, ecosystem_balance = eb,
                   carbon_stocks = cs, carbon_fluxes = cf, water_fluxes = wf, nitrogen_stocks = ns, nitrogen_fluxes = nf)
   }else {
     gamma <- list(gamma_total = gammaFull, deltaV = delta, local_change = lc, global_change = gc, ecosystem_balance = eb,
                   carbon_stocks = cs, carbon_fluxes = cf, water_fluxes = wf)
   }


   ############## export and save data if requested #############
   if (!(is.null(saveFileGamma))) {
      print(paste0("Saving gamma data to: ",saveFileGamma))
      save( deltaV, deltaVsd, gamma, file = saveFileGamma)
   }
   #
   ###
   return(gamma)
}

#' Read in output data from LPJmL to calculate the ecosystem change metric gamma
#'
#' Utility function to read in output data from LPJmL for calculation of gamma
#'
#' @param folderRef folder of reference run
#' @param folderRef2 2nd folder of reference run, if split (default NULL -- no split)
#' @param folderScen folder of scenario run
#' @param saveFile file to save read in data to (default NULL)
#' @param export flag whether to export réad in data to global environment (default F)
#' @param combined read in soilc+litc, evap+interc, rh+harvestc separate or as combined outputs (default F -- separate)
#' @param nitrogen include nitrogen outputs for pools and fluxes into gamma calculation (default F)
#' @param varnames data.frame with names of output files -- can be specified to account for variable file names (default NULL -- standard names are used)
#' @param headerout headersize of the output files (default 0)
#' @param timespan_full_ref full timespan of output files in reference period
#' @param timespan_full_scen full timespan of output files in scenario period
#' @param timespan_focus_ref specific years to use as reference period
#' @param timespan_focus_ref2 specific years to use as 2nd part of reference period in case this is split (default NULL)
#' @param timespan_focus_scen specific years to use as scenario period
#' @param npfts number of natural plant functional types (== bands in fpc - 1)
#' @param ncfts number of crop functional types
#' @param nbfts number of bfts (bioenergy functional types)
#' @param ncells number of cells in lpjml grid
#'
#' @return list data object containing arrays of state_ref, mean_state_ref,
#'         state_scen, mean_state_scen, fpc_ref, fpc_scen, bft_ref, bft_scen,
#'         cft_ref, cft_scen, lat, lon, cellArea
#'
#' @examples
#' \dontrun{
#' returned_vars <- readGammaData(folderRef = folderRef, folderRef2 = folderRef2,
#'                  folderScen = folderScen, saveFile = saveFileData, export = F,
#'                  combined = combined, nitrogen = nitrogen, varnames = varnames,
#'                  headerout = headerout, timespan_full_ref = timespan_full_ref,
#'                  timespan_full_scen = timespan_full_scen,
#'                  timespan_focus_ref = timespan_focus_ref,
#'                  timespan_focus_ref2 = timespan_focus_ref2,
#'                  timespan_focus_scen = timespan_focus_scen,
#'                  npfts = npfts, ncfts = ncfts, nbfts = nbfts, ncells = ncells)
#' }
#' @export
readGammaData <- function(folderRef, folderRef2 = NULL, folderScen, saveFile = NULL,
                          export = FALSE, combined = FALSE, varnames, headerout = 0,
                          timespan_full_ref, timespan_full_scen, timespan_focus_ref,
                          timespan_focus_ref2 = NULL, timespan_focus_scen, nitrogen,
                          npfts, ncfts, nbfts, ncells) {
   require(lpjmliotools)

   if (!(is.null(folderRef2) == is.null(timespan_focus_ref2))) {
      stop("Setting for second reference folder and timespan differ - both need to be set or unset (NULL)")
   }
   if (is.null(folderRef2)) {
      regular <- T
      nyears <- timespan_focus_ref[2] - timespan_focus_ref[1] + 1
   }else{
      regular <- F
      nyears1 <- timespan_focus_ref[2] - timespan_focus_ref[1] + 1
      nyears2 <- timespan_focus_ref2[2] - timespan_focus_ref2[1] + 1
      nyears <- nyears1 + nyears2
   }

   readGridOutputBin(inFile = paste0(folderRef,varnames["grid","outname"]), headersize = headerout, ncells = ncells)
   cellArea <- calcCellarea(lat = lat)

   ### read in lpjml output
   # for deltaV (fpc,fpc_bft,cftfrac)
   print("Reading in fpc,fpc_bft,cftfrac")
   if (regular) {# reference period not split in two folders (folderRef2==NULL)
      fpc_ref <- readCFToutput(inFile = paste0(folderRef,varnames["fpc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,bands = (npfts + 1),
                               headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      bft_ref <- readCFToutput(inFile = paste0(folderRef,varnames["fpc_bft","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,bands = nbfts,
                               headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      cft_ref <- readCFToutput(inFile = paste0(folderRef,varnames["cftfrac","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,bands = ncfts,
                               headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
   }else {# reference period split in two folders
      fpc_ref <- array(0, dim = c(ncells,(npfts + 1), nyears))
      fpc_ref[,,1:nyears1] <- readCFToutput(inFile = paste0(folderRef,varnames["fpc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,bands = (npfts + 1),
                                            headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      fpc_ref[,,(nyears1 + 1):nyears] <- readCFToutput(inFile = paste0(folderRef2,varnames["fpc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,bands = (npfts + 1),
                                                     headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
      bft_ref <- array(0,dim = c(ncells, nbfts, nyears))
      bft_ref[,,1:nyears1] <- readCFToutput(inFile = paste0(folderRef,varnames["fpc_bft","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,bands = nbfts,
                                            headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      bft_ref[,,(nyears1 + 1):nyears] <- readCFToutput(inFile = paste0(folderRef2,varnames["fpc_bft","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,bands = nbfts,
                                                     headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
      cft_ref <- array(0,dim = c(ncells, ncfts, nyears))
      cft_ref[,,1:nyears1] <- readCFToutput(inFile = paste0(folderRef,varnames["cftfrac","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,bands = ncfts,
                                            headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      cft_ref[,,(nyears1 + 1):nyears] <- readCFToutput(inFile = paste0(folderRef2,varnames["cftfrac","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,bands = ncfts,
                                                     headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
   }
   fpc_scen <- readCFToutput(inFile = paste0(folderScen,varnames["fpc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,bands = (npfts + 1),
                             headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
   bft_scen <- readCFToutput(inFile = paste0(folderScen,varnames["fpc_bft","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,bands = nbfts,
                             headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
   cft_scen <- readCFToutput(inFile = paste0(folderScen,varnames["cftfrac","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,bands = ncfts,
                             headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)

   # cffiles = ( firec rh_harvest npp ) yearly carbon fluxes
   print("Reading in firec, rh_harvest, npp")
   if (regular) {
     if (varnames["firec","timestep"] == "Y") {
       firec_ref <- readYearly(inFile = paste0(folderRef,varnames["firec","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["firec","timestep"] == "M") {
       firec_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["firec","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }
      if (combined) {
         rh_harvest_ref <- readYearly(inFile = paste0(folderRef,varnames["rh_harvest","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                      headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      }else{
        if (varnames["rh","timestep"] == "Y") {
          rh_ref <- readYearly(inFile = paste0(folderRef,varnames["rh","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
        }else if (varnames["rh","timestep"] == "M") {
          rh_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["rh","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
        }
        if (varnames["harvestc","timestep"] == "Y") {
          harvest_ref <- readYearly(inFile = paste0(folderRef,varnames["harvestc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                    headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
        }else if (varnames["harvestc","timestep"] == "M") {
          harvest_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["harvestc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                     headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
        }
        rh_harvest_ref <- rh_ref + harvest_ref
      }
     if (varnames["npp","timestep"] == "Y") {
       npp_ref <- readYearly(inFile = paste0(folderRef,varnames["npp","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                             headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["npp","timestep"] == "M") {
       npp_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["npp","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                             headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }

   }else{ # regular == F
      firec_ref <- array(0,dim = c(ncells, nyears))
      firec_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["firec","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                          headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      firec_ref[,(nyears1 + 1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["firec","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                   headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
      if (combined) {
         rh_harvest_ref <- array(0,dim = c(ncells, nyears))
         rh_harvest_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["rh_harvest","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                                  headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
         rh_harvest_ref[,(nyears1 + 1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["rh_harvest","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                           headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)

      }else{
         rh_ref <- array(0, dim = c(ncells, nyears))
         rh_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["rh","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                          headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
         rh_ref[,(nyears1 + 1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["rh","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                   headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
         harvest_ref <- array(0, dim = c(ncells, nyears))
         harvest_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["harvestc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                               headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
         harvest_ref[,(nyears1 + 1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["harvestc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                        headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
         rh_harvest_ref  <-  rh_ref + harvest_ref
      }

      npp_ref <- array(0, dim = c(ncells, nyears))
      npp_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["npp","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                        headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      npp_ref[,(nyears1 + 1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["npp","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                 headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
   }

   if (combined) {
      rh_harvest_scen <- readYearly(inFile = paste0(folderScen,varnames["rh_harvest","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                    headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
   }else{
     if (varnames["rh","timestep"] == "Y") {
       rh_scen <- readYearly(inFile = paste0(folderScen,varnames["rh","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                             headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["rh","timestep"] == "M") {
       rh_scen <- apply(readMonthly(inFile = paste0(folderRef,varnames["rh","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                   headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["harvestc","timestep"] == "Y") {
       harvest_scen <- readYearly(inFile = paste0(folderScen,varnames["harvestc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                  headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["harvestc","timestep"] == "M") {
       harvest_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["harvestc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                  headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
     }
     rh_harvest_scen  <-  rh_scen + harvest_scen
   }
   if (varnames["firec","timestep"] == "Y") {
     firec_scen <- readYearly(inFile = paste0(folderScen,varnames["firec","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                              headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
   }else if (varnames["firec","timestep"] == "M") {
     firec_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["firec","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                              headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
   }
   if (varnames["npp","timestep"] == "Y") {
     npp_scen <- readYearly(inFile = paste0(folderScen,varnames["npp","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                            headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
   }else if (varnames["npp","timestep"] == "M") {
     npp_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["npp","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                            headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
   }

   # wffiles = (evapinterc runoff transp) - yearly water fluxes
   print("Reading in evapinterc, runoff, transp")
   if (regular) {
      if (combined) {
         evapinterc_ref <- readYearly(inFile = paste0(folderRef,varnames["evapinterc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                      headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      }else{
        if (varnames["evap","timestep"] == "Y") {
          evap_ref <- readYearly(inFile = paste0(folderRef,varnames["evap","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                 headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
        }else if (varnames["evap","timestep"] == "M") {
          evap_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["evap","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                 headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
        }
        if (varnames["interc","timestep"] == "Y") {
          interc_ref <- readYearly(inFile = paste0(folderRef,varnames["interc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                   headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
        }else if (varnames["interc","timestep"] == "M") {
          interc_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["interc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                   headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
        }
        evapinterc_ref <- evap_ref + interc_ref
      }

     if (varnames["runoff","timestep"] == "Y") {
       runoff_ref <- readYearly(inFile = paste0(folderRef,varnames["runoff","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["runoff","timestep"] == "M") {
       runoff_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["runoff","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["transp","timestep"] == "Y") {
       transp_ref <- readYearly(inFile = paste0(folderRef,varnames["transp","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["transp","timestep"] == "M") {
       transp_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["transp","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }
    }else{ # regular == F
      if (combined) {
         evapinterc_ref <- array(0, dim = c(ncells, nyears))
         evapinterc_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["evapinterc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                                  headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
         evapinterc_ref[,(nyears1 + 1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["evapinterc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                           headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
      }else{
         evap_ref <- array(0, dim = c(ncells, nyears))
         evap_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["evap","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                            headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
         evap_ref[,(nyears1 + 1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["evap","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                     headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
         interc_ref <- array(0, dim = c(ncells, nyears))
         interc_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["interc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                              headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
         interc_ref[,(nyears1 + 1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["interc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                       headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
         evapinterc_ref <- evap_ref + interc_ref
      }
      runoff_ref <- array(0, dim = c(ncells, nyears))
      runoff_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["runoff","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                           headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      runoff_ref[,(nyears1 + 1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["runoff","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                    headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
      transp_ref <- array(0, dim = c(ncells, nyears))
      transp_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["transp","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                           headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      transp_ref[,(nyears1 + 1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["transp","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                    headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
   }
   if (combined) {
      evapinterc_scen <- readYearly(inFile = paste0(folderScen,varnames["evapinterc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                    headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
   }else{
     if (varnames["evap","timestep"] == "Y") {
       evap_scen <- readYearly(inFile = paste0(folderScen,varnames["evap","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["evap","timestep"] == "M") {
       evap_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["evap","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["interc","timestep"] == "Y") {
       interc_scen <- readYearly(inFile = paste0(folderScen,varnames["interc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                 headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["interc","timestep"] == "M") {
       interc_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["interc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                 headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
     }
     evapinterc_scen <- evap_scen + interc_scen

   }

   if (varnames["runoff","timestep"] == "Y") {
     runoff_scen <- readYearly(inFile = paste0(folderScen,varnames["runoff","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
   }else if (varnames["runoff","timestep"] == "M") {
     runoff_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["runoff","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
   }
   if (varnames["transp","timestep"] == "Y") {
     transp_scen <- readYearly(inFile = paste0(folderScen,varnames["transp","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
   }else if (varnames["transp","timestep"] == "M") {
     transp_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["transp","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
   }

   # csfiles = ( soillitc vegc_avg ) #carbon pools
   print("Reading in soillitc, vegc_avg")
   if (regular) {
      if (combined) {
         soillitc_ref <- readYearly(inFile = paste0(folderRef,varnames["soillitc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                    headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      }else{
        if (varnames["soilc","timestep"] == "Y") {
          soil_ref <- readYearly(inFile = paste0(folderRef,varnames["soilc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                 headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
        }else if (varnames["soilc","timestep"] == "M") {
          soil_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["soilc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                 headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
        }
        if (varnames["litc","timestep"] == "Y") {
          litc_ref <- readYearly(inFile = paste0(folderRef,varnames["litc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                 headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
        }else if (varnames["litc","timestep"] == "M") {
          litc_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["litc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                 headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
        }
        soillitc_ref <- soil_ref + litc_ref
      }
     if (varnames["vegc","timestep"] == "Y") {
       vegc_ref <- readYearly(inFile = paste0(folderRef,varnames["vegc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                              headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["vegc","timestep"] == "M") {
       vegc_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["vegc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                              headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }
    }else{ # regular == F
      if (combined) {
         soillitc_ref <- array(0, dim = c(ncells, nyears))
         soillitc_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["soillitc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                                headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
         soillitc_ref[,(nyears1 + 1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["soillitc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                         headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
      }else{
         soil_ref <- array(0, dim = c(ncells, nyears))
         soil_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["soilc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                            headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
         soil_ref[,(nyears1 + 1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["soilc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                     headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
         litc_ref <- array(0, dim = c(ncells, nyears))
         litc_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["litc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                            headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
         litc_ref[,(nyears1 + 1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["litc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                     headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
         soillitc_ref <- soil_ref + litc_ref

      }
      vegc_ref <- array(0, dim = c(ncells, nyears))
      vegc_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["vegc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                         headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      vegc_ref[,(nyears1 + 1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["vegc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                  headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
   }
   if (combined) {
      soillitc_scen <- readYearly(inFile = paste0(folderScen,varnames["soillitc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                  headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
   }else{
     if (varnames["soilc","timestep"] == "Y") {
       soil_scen <- readYearly(inFile = paste0(folderScen,varnames["soilc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["soilc","timestep"] == "M") {
       soil_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["soilc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["litc","timestep"] == "Y") {
       litc_scen <- readYearly(inFile = paste0(folderScen,varnames["litc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["litc","timestep"] == "M") {
       litc_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["litc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
     }
     soillitc_scen <- soil_scen + litc_scen
   }
   if (varnames["litc","timestep"] == "Y") {
     vegc_scen <- readYearly(inFile = paste0(folderScen,varnames["vegc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                             headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
   }else if (varnames["litc","timestep"] == "M") {
     vegc_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["vegc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                             headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
   }

   # allfiles = (swcsum firef) # other system-internal processes
   print("Reading in swcsum, firef")
   if (regular) {
      if (combined) {
         swcsum_ref <- readYearly(inFile = paste0(folderRef,varnames["swcsum","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                  headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      }else{
        if (varnames["swc","timestep"] == "Y") {
          swcsum_ref <- apply(readCFToutput(inFile = paste0(folderRef,varnames["swc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4, bands = 6,
                                            headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
        }else if (varnames["swc","timestep"] == "M") {
          swcsum_ref <- apply(readMonthlyCFToutput(inFile = paste0(folderRef,varnames["swc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4, bands = 6,
                                            headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,4),sum)
        }
               }
     if (varnames["firef","timestep"] == "Y") {
       firef_ref <- readYearly(inFile = paste0(folderRef,varnames["firef","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["firef","timestep"] == "M") {
       firef_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["firef","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }
   }else{ # regular == F
      if (combined) {
         swcsum_ref <- array(0, dim = c(ncells, nyears))
         swcsum_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["swc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                              headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
         swcsum_ref[,(nyears1+1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["swc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                       headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
      }else{
         swcsum_ref <- array(0, dim = c(ncells, nyears))
         swcsum_ref[,1:nyears1] <- apply(readCFToutput(inFile = paste0(folderRef,varnames["swc","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4, bands = 6,
                                              headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
         swcsum_ref[,(nyears1+1):nyears] <-  apply(readCFToutput(inFile = paste0(folderRef2,varnames["swc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4, bands = 6,
                                                       headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells),c(1,3),sum)

      }
      firef_ref <- array(0, dim = c(ncells, nyears))
      firef_ref[,1:nyears1] <- readYearly(inFile = paste0(folderRef,varnames["firef","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                          headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
      firef_ref[,(nyears1+1):nyears] <- readYearly(inFile = paste0(folderRef2,varnames["firef","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                                   headersize = headerout,getyearstart = timespan_focus_ref2[1],getyearstop = timespan_focus_ref2[2],ncells = ncells)
   }
   if (combined) {
      swcsum_scen <- readYearly(inFile = paste0(folderScen,varnames["swcsum","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
   }else{
     if (varnames["swc","timestep"] == "Y") {
       swcsum_scen <- apply(readCFToutput(inFile = paste0(folderScen,varnames["swc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                       headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells, bands = 6),c(1,3),sum)
     }else if (varnames["swc","timestep"] == "M") {
       swcsum_scen <- apply(readMonthlyCFToutput(inFile = paste0(folderScen,varnames["swc","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                       headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells, bands = 6),c(1,4),sum)
     }
   }
   if (varnames["firef","timestep"] == "Y") {
     firef_scen <- readYearly(inFile = paste0(folderScen,varnames["firef","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                              headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
   }else if (varnames["firef","timestep"] == "M") {
     firef_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["firef","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                              headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
   }

   # nitrogen variables



   if (nitrogen) {
     print("Reading in n-pools: soilnh4, soilno3 + fluxes: leaching, bnf, n_volatilization, n2o_nit, n2o_denit n2_emis")

     # reference state
     # pools: soilnh4, soilno3
     if (varnames["soilnh4","timestep"] == "Y") {
       soilnh4_ref <- readYearly(inFile = paste0(folderRef,varnames["soilnh4","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                              headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["soilnh4","timestep"] == "M") {
       soilnh4_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["soilnh4","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                     headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["soilno3","timestep"] == "Y") {
       soilno3_ref <- readYearly(inFile = paste0(folderRef,varnames["soilno3","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                              headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["soilno3","timestep"] == "M") {
       soilno3_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["soilno3","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                     headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["vegn","timestep"] == "Y") {
       vegn_ref <- readYearly(inFile = paste0(folderRef,varnames["vegn","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                 headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["vegn","timestep"] == "M") {
       vegn_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["vegn","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                        headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }
     # fluxes: leaching, n2o_nit, n2o_denit n2_emis, bnf, n_volatilization
     if (varnames["leaching","timestep"] == "Y") {
       leaching_ref <- readYearly(inFile = paste0(folderRef,varnames["leaching","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                              headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["leaching","timestep"] == "M") {
       leaching_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["leaching","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                     headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["n2o_denit","timestep"] == "Y") {
       n2o_denit_ref <- readYearly(inFile = paste0(folderRef,varnames["n2o_denit","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                              headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["n2o_denit","timestep"] == "M") {
       n2o_denit_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["n2o_denit","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                     headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["n2o_nit","timestep"] == "Y") {
       n2o_nit_ref <- readYearly(inFile = paste0(folderRef,varnames["n2o_nit","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                              headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["n2o_nit","timestep"] == "M") {
       n2o_nit_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["n2o_nit","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                     headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["n2_emis","timestep"] == "Y") {
       n2_emis_ref <- readYearly(inFile = paste0(folderRef,varnames["n2_emis","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                 headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["n2_emis","timestep"] == "M") {
       n2_emis_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["n2_emis","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                        headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["bnf","timestep"] == "Y") {
       bnf_ref <- readYearly(inFile = paste0(folderRef,varnames["bnf","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                 headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["bnf","timestep"] == "M") {
       bnf_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["bnf","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                        headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["n_volatilization","timestep"] == "Y") {
       n_volatilization_ref <- readYearly(inFile = paste0(folderRef,varnames["n_volatilization","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                 headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells)
     }else if (varnames["n_volatilization","timestep"] == "M") {
       n_volatilization_ref <- apply(readMonthly(inFile = paste0(folderRef,varnames["n_volatilization","outname"]),startyear = timespan_full_ref[1],stopyear = timespan_full_ref[2],size = 4,
                                        headersize = headerout,getyearstart = timespan_focus_ref[1],getyearstop = timespan_focus_ref[2],ncells = ncells),c(1,3),sum)
     }#

     # scenario state
     # pools
     if (varnames["soilnh4","timestep"] == "Y") {
       soilnh4_scen <- readYearly(inFile = paste0(folderScen,varnames["soilnh4","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["soilnh4","timestep"] == "M") {
       soilnh4_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["soilnh4","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                      headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["soilno3","timestep"] == "Y") {
       soilno3_scen <- readYearly(inFile = paste0(folderScen,varnames["soilno3","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["soilno3","timestep"] == "M") {
       soilno3_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["soilno3","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                      headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["vegn","timestep"] == "Y") {
       vegn_scen <- readYearly(inFile = paste0(folderScen,varnames["vegn","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                  headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["vegn","timestep"] == "M") {
       vegn_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["vegn","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                         headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
     }
     # fluxes
     if (varnames["leaching","timestep"] == "Y") {
       leaching_scen <- readYearly(inFile = paste0(folderScen,varnames["leaching","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["leaching","timestep"] == "M") {
       leaching_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["leaching","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                      headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["n2o_denit","timestep"] == "Y") {
       n2o_denit_scen <- readYearly(inFile = paste0(folderScen,varnames["n2o_denit","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["n2o_denit","timestep"] == "M") {
       n2o_denit_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["n2o_denit","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                      headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["n2o_nit","timestep"] == "Y") {
       n2o_nit_scen <- readYearly(inFile = paste0(folderScen,varnames["n2o_nit","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                               headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["n2o_nit","timestep"] == "M") {
       n2o_nit_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["n2o_nit","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                      headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["n2_emis","timestep"] == "Y") {
       n2_emis_scen <- readYearly(inFile = paste0(folderScen,varnames["n2_emis","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                   headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["n2_emis","timestep"] == "M") {
       n2_emis_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["n2_emis","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                          headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["bnf","timestep"] == "Y") {
       bnf_scen <- readYearly(inFile = paste0(folderScen,varnames["bnf","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                    headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["bnf","timestep"] == "M") {
       bnf_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["bnf","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                           headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
     }
     if (varnames["n_volatilization","timestep"] == "Y") {
       n_volatilization_scen <- readYearly(inFile = paste0(folderScen,varnames["n_volatilization","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                  headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells)
     }else if (varnames["n_volatilization","timestep"] == "M") {
       n_volatilization_scen <- apply(readMonthly(inFile = paste0(folderScen,varnames["n_volatilization","outname"]),startyear = timespan_full_scen[1],stopyear = timespan_full_scen[2],size = 4,
                                         headersize = headerout,getyearstart = timespan_focus_scen[1],getyearstop = timespan_focus_scen[2],ncells = ncells),c(1,3),sum)
     }
     # Calculating compound n emissions vector
     aggregated_n_emissions_scen <- n_volatilization_scen + n2o_nit_scen + n2o_denit_scen + n2_emis_scen
     aggregated_n_emissions_ref <- n_volatilization_ref + n2o_nit_ref + n2o_denit_ref + n2_emis_ref
     soiln_scen <- soilno3_scen + soilnh4_scen
     soiln_ref <- soilno3_ref + soilnh4_ref

   }
   if (nitrogen) {
     mean_state_ref <- cbind(rowMeans(firec_ref),rowMeans(rh_harvest_ref),rowMeans(npp_ref),rowMeans(evapinterc_ref),rowMeans(runoff_ref),
                             rowMeans(transp_ref),rowMeans(soillitc_ref),rowMeans(vegc_ref),rowMeans(swcsum_ref),rowMeans(firef_ref),
                             rowMeans(soiln_ref), rowMeans(vegn_ref), rowMeans(leaching_ref), rowMeans(bnf_ref), rowMeans(aggregated_n_emissions_ref))
     mean_state_scen <- cbind(rowMeans(firec_scen),rowMeans(rh_harvest_scen),rowMeans(npp_scen),rowMeans(evapinterc_scen),rowMeans(runoff_scen),
                              rowMeans(transp_scen),rowMeans(soillitc_scen),rowMeans(vegc_scen),rowMeans(swcsum_scen),rowMeans(firef_scen),
                              rowMeans(soiln_scen), rowMeans(vegn_scen), rowMeans(leaching_scen), rowMeans(bnf_scen), rowMeans(aggregated_n_emissions_scen))
     require(abind) #       1           2             3          4                5          6          7             8          9           10        11             12            13            14           15
     state_ref <- abind( firec_ref, rh_harvest_ref, npp_ref, evapinterc_ref, runoff_ref, transp_ref, soillitc_ref, vegc_ref, swcsum_ref, firef_ref, soiln_ref,  vegn_ref,  leaching_ref,  bnf_ref,  aggregated_n_emissions_ref,along = 3)
     state_scen <- abind(firec_scen,rh_harvest_scen,npp_scen,evapinterc_scen,runoff_scen,transp_scen,soillitc_scen,vegc_scen,swcsum_scen,firef_scen,soiln_scen, vegn_scen, leaching_scen, bnf_scen, aggregated_n_emissions_scen,along = 3)
   }else{
     mean_state_ref <- cbind(rowMeans(firec_ref),rowMeans(rh_harvest_ref),rowMeans(npp_ref),rowMeans(evapinterc_ref),rowMeans(runoff_ref),
                             rowMeans(transp_ref),rowMeans(soillitc_ref),rowMeans(vegc_ref),rowMeans(swcsum_ref),rowMeans(firef_ref))
     mean_state_scen <- cbind(rowMeans(firec_scen),rowMeans(rh_harvest_scen),rowMeans(npp_scen),rowMeans(evapinterc_scen),rowMeans(runoff_scen),
                              rowMeans(transp_scen),rowMeans(soillitc_scen),rowMeans(vegc_scen),rowMeans(swcsum_scen),rowMeans(firef_scen))
     require(abind)
     state_ref <- abind(firec_ref,rh_harvest_ref,npp_ref,evapinterc_ref,runoff_ref,transp_ref,soillitc_ref,vegc_ref,swcsum_ref,firef_ref,along = 3)
     state_scen <- abind(firec_scen,rh_harvest_scen,npp_scen,evapinterc_scen,runoff_scen,transp_scen,soillitc_scen,vegc_scen,swcsum_scen,firef_scen,along = 3)
   }

   if (!(is.null(saveFile))) {
     print(paste0("Saving data to: ",saveFile))
     save(state_ref,mean_state_ref,state_scen,mean_state_scen,fpc_ref,fpc_scen,bft_ref,bft_scen,cft_ref,cft_scen,lat,lon,cellArea,file = saveFile)
   }
   if (export) {
      print(paste0("Creating return data list object"))
      exportVars2GlobalEnv(state_ref,mean_state_ref,state_scen,mean_state_scen,fpc_ref,fpc_scen,bft_ref,bft_scen,cft_ref,cft_scen,lat,lon,cellArea)
   }
  return(list(state_ref = state_ref, mean_state_ref = mean_state_ref,
                 state_scen = state_scen, mean_state_scen = mean_state_scen,
                 fpc_ref = fpc_ref, fpc_scen = fpc_scen, bft_ref = bft_ref,
                 bft_scen = bft_scen, cft_ref = cft_ref, cft_scen = cft_scen,lat = lat, lon = lon, cellArea = cellArea))
}

#' Calculates changes in vegetation structure (deltaV)
#'
#' Utility function to calculate changes in vegetation structure (deltaV)
#' for calculation of gamma
#'
#' @param fpcRef reference fpc array (dim: [ncells,npfts+1])
#' @param fpcScen scenario fpc array (dim: [ncells,npfts+1])
#' @param bftRef reference bft array (dim: [ncells,nbfts])
#' @param bftScen scenario bft array (dim: [ncells,nbfts])
#' @param cftRef reference cft array (dim: [ncells,ncfts])
#' @param cftScen scenario cft array (dim: [ncells,ncfts])
#' @param weighting apply "old" (Ostberg-like), "new", or "equal" weighting of deltaV weights (default "old")
#'
#' @return deltaV array of size ncells with the deltaV value [0,1] for each cell
#'
#' @examples
#' \dontrun{
#' deltaV <- calcDeltaV(fpcRef = fpc_ref_mean,fpcScen = apply(fpc_scen,c(1,2),mean),
#'           bftRef = bft_ref_mean,bftScen = apply(bft_scen,c(1,2),mean),
#'           cftRef = cft_ref_mean, cftScen = apply(cft_scen,c(1,2),mean),
#'           weighting = "new")
#' }
#' @export
calcDeltaV <- function(fpcRef, fpcScen, bftRef, bftScen, cftRef, cftScen,
                       weighting) {
   require(lpjmliotools)
   di <- dim(fpcRef)
   ncells <- di[1]
   npfts <- di[2] - 1

   fpcRef[fpcRef < 0] = 0
   fpcScen[fpcScen < 0] = 0
   bftRef[bftRef < 0] = 0
   bftScen[bftScen < 0] = 0
   cftRef[cftRef < 0] = 0
   cftScen[cftScen < 0] = 0
   if (npfts == 9) {
      # barren     =      crop area             natural vegetation area               barren under bioenergy trees
      barrenAreaRef <- 1 - rowSums(cftRef) - rowSums(fpcRef[,2:10])*fpcRef[,1] + rowSums(cftRef[,c(16,32)])*(1 - rowSums(bftRef[,c(1:4,7:10)]))
      barrenAreaRef[barrenAreaRef < 0] <- 0
      treeAreaRef <- array(0,dim = c(ncells,11))
      treeAreaRef[,1:7] <- fpcRef[,2:8]*fpcRef[,1] # natural tree area fractions scaled by total natural frac
      treeAreaRef[,8:9] <- cftRef[,16]*bftRef[,1:2]/rowSums(bftRef[,c(1,2,4)]) # # fraction of rainfed tropical and temperate BE trees scaled by total rainfed bioenergy tree area and relative fpc of bioenergy trees and grass under bioenergy trees
      treeAreaRef[,10:11] <- cftRef[,32]*bftRef[,7:8]/rowSums(bftRef[,c(7,8,10)]) # fraction of irrigated tropical and temperate BE trees scaled by total irrigated bioenergy tree area and relative fpc of bioenergy trees and grass under bioenergy trees
      grassAreaRef <- array(0,dim = c(ncells,20))
      grassAreaRef[,1:2] <- fpcRef[,9:10]*fpcRef[,1] # natural grass
      grassAreaRef[,3:15] <- cftRef[,1:13] + cftRef[,17:29] # crops
      grassAreaRef[,16] <- cftRef[,14] # managed grass rf
      grassAreaRef[,17] <- cftRef[,30] # managed grass irr
      grassAreaRef[,18] <- cftRef[,15] + cftRef[,31] # bioenergy grass
      grassAreaRef[,19] <- cftRef[,16]*bftRef[,4]/rowSums(bftRef[,c(1,2,4)]) # fraction of rainfed grass under bioenergy trees
      grassAreaRef[,20] <- cftRef[,32]*bftRef[,10]/rowSums(bftRef[,c(7,8,10)]) # fraction of irrigated grass under bioenergy trees

      barrenAreaScen <- 1 - rowSums(cftScen) - rowSums(fpcScen[,2:10])*fpcScen[,1] + rowSums(cftScen[,c(16,32)])*(1 - rowSums(bftScen[,c(1:4,7:10)])) # barren
      barrenAreaScen[barrenAreaScen < 0] = 0
      treeAreaScen <- array(0,dim = c(ncells,11))
      treeAreaScen[,1:7] <- fpcScen[,2:8]*fpcScen[,1] # natural tree area fractions scaled by total natural frac
      treeAreaScen[,8:9] <- cftScen[,16]*bftScen[,1:2]/rowSums(bftScen[,c(1,2,4)]) # fraction of rainfed tropical and temperate BE trees scaled by total rainfed bioenergy tree area and relative fpc of bioenergy trees and grass under bioenergy trees
      treeAreaScen[,10:11] <- cftScen[,32]*bftScen[,7:8]/rowSums(bftScen[,c(7,8,10)]) # fraction of irrigated tropical and temperate BE trees scaled by total irrigated bioenergy tree area and relative fpc of bioenergy trees and grass under bioenergy trees
      grassAreaScen <- array(0,dim = c(ncells,20))
      grassAreaScen[,1:2] <- fpcScen[,9:10]*fpcScen[,1] # natural grass
      grassAreaScen[,3:15] <- cftScen[,1:13] + cftScen[,17:29] # crops
      grassAreaScen[,16] <- cftScen[,14] # managed grass rf
      grassAreaScen[,17] <- cftScen[,30] # managed grass irr
      grassAreaScen[,18] <- cftScen[,15] + cftScen[,31] # bioenergy grass
      grassAreaScen[,19] <- cftScen[,16]*bftScen[,4]/rowSums(bftScen[,c(1,2,4)]) # fraction of rainfed grass under bioenergy trees
      grassAreaScen[,20] <- cftScen[,32]*bftScen[,10]/rowSums(bftScen[,c(7,8,10)]) # fraction of irrigated grass under bioenergy trees

      treeAttributes <- matrix(c( #evergreenness,needleleavedness,tropicalness,borealness,naturalness
         c(  1,    0,  1,  0,  1),#1 TrBE
         c(  0,    0,  1,  0,  1),#2 TrBR
         c(  1,    1,  0,  0,  1),#3 TeNE
         c(  1,    0,  0,  0,  1),#4 TeBE
         c(  0,    0,  0,  0,  1),#5 TeBS
         c(  1,    1,  0,  1,  1),#6 BoNE
         c(  0, 0.25,  0,  1,  1),#7 BoS (including larchs)
         c(  1,    0,  1,  0,  0),#8 TrBi tropical bioenergy rainfed
         c(  0,    0,  0,  0,  0),#9 TeBi temperate bioenergy rainfed
         c(  1,    0,  1,  0,  0),#10 TrBi tropical bioenergy irrigated
         c(  0,    0,  0,  0,  0) #11 TeBi temperate bioenergy irrigated
      ),nrow = 11,byrow = T)

      if (weighting == "equal") {
        treeWeights <- c(0.2, 0.2 ,0.2, 0.2, 0.2)
      }else if (weighting == "new") { # changed compared to Sebastian
        treeWeights <- c(0.2, 0.2 ,0.3, 0.3, 0.3)/1.3
      }else if (weighting == "old") { # Sebastian's method (no downscaling to weightsum 1)
        treeWeights <- c(0.2, 0.2 ,0.3, 0.3, 0.3)
      }else{
        stop("Unknown method of weighting.")
      }

      grassAttributes <- array(0,dim = c(ncells,20,2))
      #1 C3grass
      #2 C4grass
      #3 TemperateCereals
      #4 Rice
      #5 Maize
      #6 TropicalCereals
      #7 Pulses
      #8 TemperateRoots
      #9 TropicalRoots
      #10 Sunflower
      #11 Soybean
      #12 Groundnut
      #13 Rapeseed
      #14 Sugarcane
      #15 Others
      #16 Managed grass rainfed
      #17 Managed grass irrigated
      #18 Bioenergy grass
      #19 Grass under rainfed Bioenergy trees
      #20 Grass under irrigated Bioenergy trees
      grassAttributes[,,1] <- rep(c(0,1,0,1,1,1,0.5,0,1,0.5,1,1,0.5,1,0.5,NA,NA,1,NA,NA),each = ncells) #tropicalness
      grassAttributes[,,2] <- rep(c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),each = ncells) #naturalness
      dynGrassAttributes <- cbind(bftScen[,6]/rowSums(bftScen[,5:6]),bftScen[,12]/rowSums(bftScen[,11:12])) # dynamic share of tropicalness for rf/irr grasslands taken from ratio of bioenergy grasses
      dynGrassAttributes[!is.finite(dynGrassAttributes)] <- 0
      grassAttributes[,16:17,1] <- dynGrassAttributes # managed grass rf/irr
      grassAttributes[,19:20,1] <- dynGrassAttributes # grass under biotrees rf/irr (taken from managed grass)

      if (weighting == "equal") {
        grassWeights <- c(0.2, 0.2)
      }else if (weighting == "new") { # changed compared to Sebastian
        grassWeights <- c(0.5,0.5)
      }else if (weighting == "old") { # Sebastian's method (no downscaling to weightsum 1)
        grassWeights <- c(0.3,0.3)
      }else{
        stop("Unknown method of weighting.")
      }

   }else if (npfts == 11) {
      # barren     =      crop area             natural vegetation area               barren under bioenergy trees
      barrenAreaRef <- 1 - rowSums(cftRef) - rowSums(fpcRef[,2:12])*fpcRef[,1] + rowSums(cftRef[,c(16,32)])*(1 - rowSums(bftRef[,c(4:9,13:18)]))
      barrenAreaRef[barrenAreaRef < 0] <- 0
      treeAreaRef <- array(0,dim = c(ncells,12))
      treeAreaRef[,1:8] <- fpcRef[,2:9]*fpcRef[,1] # natural tree area fractions scaled by total natural frac
      treeAreaRef[,9:10] <- cftRef[,16]*bftRef[,7:8]/rowSums(bftRef[,4:8]) # fraction of rainfed tropical and temperate BE trees scaled by total rainfed bioenergy tree area and relative fpc of bioenergy trees and grass under bioenergy trees
      treeAreaRef[,11:12] <- cftRef[,32]*bftRef[,16:17]/rowSums(bftRef[,13:17]) # fraction of irrigated tropical and temperate BE trees scaled by total irrigated bioenergy tree area and relative fpc of bioenergy trees and grass under bioenergy trees
      grassAreaRef <- array(0,dim = c(ncells,21))
      grassAreaRef[,1:3] <- fpcRef[,10:12]*fpcRef[,1] # natural grass
      grassAreaRef[,4:16] <- cftRef[,1:13] + cftRef[,17:29] # crops
      grassAreaRef[,17] <- cftRef[,14] # managed grass rf
      grassAreaRef[,18] <- cftRef[,30] # managed grass irr
      grassAreaRef[,19] <- cftRef[,15] + cftRef[,31] # bioenergy grass
      grassAreaRef[,20] <- cftRef[,16]*rowSums(bftRef[,4:6])/rowSums(bftRef[,4:8]) # fraction of rainfed grass under bioenergy trees
      grassAreaRef[,21] <- cftRef[,32]*rowSums(bftRef[,13:15])/rowSums(bftRef[,13:17]) # fraction of irrigated grass under bioenergy trees

      # barren     =      crop area             natural vegetation area               barren under bioenergy trees
      barrenAreaScen <- 1 - rowSums(cftScen) - rowSums(fpcScen[,2:12])*fpcScen[,1] + rowSums(cftScen[,c(16,32)])*(1 - rowSums(bftScen[,c(4:9,13:18)]))
      barrenAreaScen[barrenAreaScen < 0] <- 0
      treeAreaScen <- array(0,dim = c(ncells,12))
      treeAreaScen[,1:8] <- fpcScen[,2:9]*fpcScen[,1] # natural tree area fractions scaled by total natural frac
      treeAreaScen[,9:10] <- cftScen[,16]*bftScen[,7:8]/rowSums(bftScen[,4:8]) # fraction of rainfed tropical and temperate BE trees scaled by total rainfed bioenergy tree area and relative fpc of bioenergy trees and grass under bioenergy trees
      treeAreaScen[,11:12] <- cftScen[,32]*bftScen[,16:17]/rowSums(bftScen[,13:17]) # fraction of irrigated tropical and temperate BE trees scaled by total irrigated bioenergy tree area and relative fpc of bioenergy trees and grass under bioenergy trees
      grassAreaScen <- array(0,dim = c(ncells,21))
      grassAreaScen[,1:3] <- fpcScen[,10:12]*fpcScen[,1] # natural grass
      grassAreaScen[,4:16] <- cftScen[,1:13] + cftScen[,17:29] # crops
      grassAreaScen[,17] <- cftScen[,14] # managed grass rf
      grassAreaScen[,18] <- cftScen[,30] # managed grass irr
      grassAreaScen[,19] <- cftScen[,15] + cftScen[,31] # bioenergy grass
      grassAreaScen[,20] <- cftScen[,16]*rowSums(bftScen[,4:6])/rowSums(bftScen[,4:8]) # fraction of rainfed grass under bioenergy trees
      grassAreaScen[,21] <- cftScen[,32]*rowSums(bftScen[,13:15])/rowSums(bftScen[,13:17]) # fraction of irrigated grass under bioenergy trees

      treeAttributes <- matrix(c( #evergreenness,needleleavedness,tropicalness,borealness,naturalness
         c(  1,    0,  1,  0,  1),#1 TrBE
         c(  0,    0,  1,  0,  1),#2 TrBR
         c(  1,    1,  0,  0,  1),#3 TeNE
         c(  1,    0,  0,  0,  1),#4 TeBE
         c(  0,    0,  0,  0,  1),#5 TeBS
         c(  1,    1,  0,  1,  1),#6 BoNE
         c(  0,    0,  0,  1,  1),#7 BoBS
         c(  0,    1,  0,  1,  1),#8 BoNS
         c(  1,    0,  1,  0,  0),#9 TrBi tropical bioenergy rainfed
         c(  0,    0,  0,  0,  0),#10 TeBi temperate bioenergy rainfed
         c(  1,    0,  1,  0,  0),#11 TrBi tropical bioenergy irrigated
         c(  0,    0,  0,  0,  0) #12 TeBi temperate bioenergy irrigated
      ),nrow = 12,byrow = T)
      if (weighting == "equal") {
        treeWeights <- c(0.2, 0.2, 0.2, 0.2, 0.2)
      }else if (weighting == "new") { # changed compared to Sebastian
        treeWeights <- c(0.2, 0.2, 0.3, 0.3, 0.3)/1.3
      }else if (weighting == "old") { # Sebastian's method (no downscaling to weightsum 1)
        treeWeights <- c(0.2, 0.2, 0.3, 0.3, 0.3)
      }else{
        stop("Unknown method of weighting.")
      }

      grassAttributes <- array(0,dim = c(ncells,21,3))
      #1 C4grass tropic
      #2 C3grass temperate
      #3 C3grass polar
      #4 TemperateCereals
      #5 Rice
      #6 Maize
      #7 TropicalCereals
      #8 Pulses
      #9 TemperateRoots
      #10 TropicalRoots
      #11 Sunflower
      #12 Soybean
      #13 Groundnut
      #14 Rapeseed
      #15 Sugarcane
      #16 Others
      #17 Managed grass rainfed
      #18 Managed grass irrigated
      #19 Bioenergy grass
      #20 Grass under rainfed Bioenergy trees
      #21 Grass under irrigated Bioenergy trees
      grassAttributes[,,1] <- rep(c(1,0,0,0,1,1,1,0.5,0,1,0.5,1,1,0.5,1,0.5,NA,NA,1,NA,NA), each = ncells) #tropicalness
      grassAttributes[,,2] <- rep(c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,NA,NA,0,NA,NA), each = ncells) #borealness
      grassAttributes[,,3] <- rep(c(1,1,1,0,0,0,0,0  ,0,0,0  ,0,0,0  ,0,0  ,0 ,0 ,0,0 ,0), each = ncells) #naturalness
      dynTropGrassAttributes <- cbind(bftScen[,1]/rowSums(bftScen[,1:3]),    # dynamic share of tropicalness for rf grasslands
                                   bftScen[,10]/rowSums(bftScen[,10:12]), # dynamic share of tropicalness for irr grasslands
                                   bftScen[,4]/rowSums(bftScen[,4:6]),    # dynamic share of tropicalness for grass under rf biotrees
                                   bftScen[,13]/rowSums(bftScen[,13:15])) # dynamic share of tropicalness for grass under irr biotrees
      dynTropGrassAttributes[!is.finite(dynTropGrassAttributes)] <- 0
      grassAttributes[,c(17,18,20,21),1] <- dynTropGrassAttributes # managed grass rf/irr, grass under biotrees rf/irr
      dynBorealGrassAttributes <- cbind(bftScen[,3]/rowSums(bftScen[,1:3]),    # dynamic share of borealness for rf grasslands
                                     bftScen[,12]/rowSums(bftScen[,10:12]), # dynamic share of borealness for irr grasslands
                                     bftScen[,6]/rowSums(bftScen[,4:6]),    # dynamic share of borealness for grass under rf biotrees
                                     bftScen[,15]/rowSums(bftScen[,13:15])) # dynamic share of borealness for grass under irr biotrees
      dynBorealGrassAttributes[!is.finite(dynBorealGrassAttributes)] <- 0
      grassAttributes[,c(17,18,20,21),2] <- dynBorealGrassAttributes # managed grass rf/irr, grass under biotrees rf/irr
      if (weighting == "equal"){
        grassWeights <- c(0.2, 0.2, 0.2)
      }else if (weighting == "old" || weighting == "new"){
        grassWeights <- c(0.3333333,0.3333333,0.3333333)
      }else{
        stop("Unknown method of weighting.")
      }
   }else{
      stop("Unknown number of pfts.")
   }

   # compute deltaV
   require(fBasics) # for rowMins
   barrenV <- rowMins(cbind(barrenAreaRef,barrenAreaScen))
   treesV <- rowMins(cbind(rowSums(treeAreaRef,na.rm = T),rowSums(treeAreaScen,na.rm = T)))
   grassV <- rowMins(cbind(rowSums(grassAreaRef,na.rm = T),rowSums(grassAreaScen,na.rm = T)))
   innerSumTrees <- abs(rowSums(treeAreaRef[,]*rep(treeAttributes[,1],each = ncells),na.rm = T) - rowSums(treeAreaScen[,]*rep(treeAttributes[,1],each = ncells),na.rm = T))*treeWeights[1] +       # evergreenness
                  abs(rowSums(treeAreaRef[,]*rep(treeAttributes[,2],each = ncells),na.rm = T) - rowSums(treeAreaScen[,]*rep(treeAttributes[,2],each = ncells),na.rm = T))*treeWeights[2] +       # needleleavedness
                  abs(rowSums(treeAreaRef[,]*rep(treeAttributes[,3],each = ncells),na.rm = T) - rowSums(treeAreaScen[,]*rep(treeAttributes[,3],each = ncells),na.rm = T))*treeWeights[3] +       # tropicalness
                  abs(rowSums(treeAreaRef[,]*rep(treeAttributes[,4],each = ncells),na.rm = T) - rowSums(treeAreaScen[,]*rep(treeAttributes[,4],each = ncells),na.rm = T))*treeWeights[4] +       # borealness
                  abs(rowSums(treeAreaRef[,]*rep(treeAttributes[,5],each = ncells),na.rm = T) - rowSums(treeAreaScen[,]*rep(treeAttributes[,5],each = ncells),na.rm = T))*treeWeights[5]         # naturalness
   if (npfts == 9) {
      innerSumGrasses <- abs(rowSums(grassAreaRef[,]*grassAttributes[,,1],na.rm = T) - rowSums(grassAreaScen[,]*grassAttributes[,,1],na.rm = T))*grassWeights[1] + # tropicalness
                        abs(rowSums(grassAreaRef[,]*grassAttributes[,,2],na.rm = T) - rowSums(grassAreaScen[,]*grassAttributes[,,2],na.rm = T))*grassWeights[2]   # naturalness
   }else if (npfts == 11) {
      innerSumGrasses <- abs(rowSums(grassAreaRef[,]*grassAttributes[,,1],na.rm = T) - rowSums(grassAreaScen[,]*grassAttributes[,,1],na.rm = T))*grassWeights[1] + # tropicalness
                        abs(rowSums(grassAreaRef[,]*grassAttributes[,,2],na.rm = T) - rowSums(grassAreaScen[,]*grassAttributes[,,2],na.rm = T))*grassWeights[2] + # borealness
                        abs(rowSums(grassAreaRef[,]*grassAttributes[,,3],na.rm = T) - rowSums(grassAreaScen[,]*grassAttributes[,,3],na.rm = T))*grassWeights[3]   # naturalness
   }else{
      stop("Unknown number of pfts.")
   }
   deltaV <- 1 - barrenV - treesV*(1 - innerSumTrees) - grassV*(1 - innerSumGrasses)
   deltaV[deltaV < 0] <- 0
   deltaV[!is.finite(deltaV)] <- 0
   #plotGammaToScreen(data=deltaV,title = "deltaV")
   return(deltaV)
}

################# further gamma utility functions ##################
T_sigmoid_Trafo <- function(x) {
   return( -1/exp(3) + (1 + 1/exp(3))/(1 + exp(-6*(x - 0.5))) )
}
balance <- function(v1, v2) {
   return( 1 - sum(v1*v2)/( sqrt(sum(v1*v1)) * sqrt(sum(v2*v2)) ))
}
std_cellwise <- function(a) {
   return(apply(a,1,sd))
}
global_yearly_weighted_mean <- function(a, cellArea) {
   # a is matrix with dim=c(cells,years)
   # cellArea the corresponding cellArea array with dim=c(cells)
   return( sum(a*cellArea, na.rm = T)/(length(a[1,])*sum(cellArea, na.rm = T)) )
}
globally_weighted_mean_foreach_var <- function(x, cellArea) {
   # x is matrix with dim=c(ncells,vars)
   # cellArea the corresponding cellArea array to x with dim=c(ncells)
   return( colSums(x*cellArea, na.rm = T)/sum(cellArea, na.rm = T) )
}
S_change_to_var_ratio <- function(x, s) {
   return( 1/(1 + exp(-4*(x/s - 2))) )
}
state_Diff_local <- function(ref, scen, epsilon = 10^-4) { # based on Heyder 2011 eq. 6-9; epsilon case handling from code by Ostberg (not documented in papers)
   # param ref mean reference state vector of dimension c(ncells,variables)
   # param scen mean scenario state vector of dimension c(ncells,variables)
   # param epsilon threshold for variables to be treated as 0

   # Ostberg code: case change_metric_lu_comparison_jun2013.c
   di <- dim(ref)
   s_scen <- scen/ref #generally normalize the scenario state vector by the reference state
   s_ref <- array(1, dim = di) #initialize

   # for variables in places, where ref is small (<epsilon), but scen larger (Ostberg code, line 798)
   # set back scenario and reference vector, to keep the unscaled values (Ostberg code, line 804)
   s_scen[abs(ref) < epsilon & abs(scen) > epsilon] = scen[abs(ref) < epsilon & abs(scen) > epsilon]
   s_ref[abs(ref) < epsilon & abs(scen) > epsilon] = ref[abs(ref) < epsilon & abs(scen) > epsilon]

   # for variables in places, where ref and scen are small (<epsilon), return 0 (both are 1, difference is 0) (Ostberg code, line 809)
   s_scen[abs(ref) < epsilon & abs(scen) < epsilon] <- 1 # no change

   # normalize both state vectors by the sqrt(amount of state variables) to ensure length(s_ref)==1
   # (this is part of the weighting in the Ostberg code)
   s_ref <- s_ref/sqrt(di[2])
   s_scen <- s_scen/sqrt(di[2])

   return( sqrt(rowSums((s_scen - s_ref)*(s_scen - s_ref))) ) #length of the local difference vector s_scen (sl2) - s_ref (sl1)
}
state_Diff_global <- function(ref, scen, cellArea) { #c based on Heyder 2011 eq. 10-13
   # param ref mean reference state vector of dimension c(ncells,variables)
   # param scen mean scenario state vector of dimension c(ncells,variables)
   # param cellArea area of each cell as a vector of dim=c(ncells)
   # returns the length of the difference vector for each cell
   di <- dim(ref)
   ncells <- di[1]
   global_mean_ref <- globally_weighted_mean_foreach_var(ref,cellArea)
   global_mean_scen <- globally_weighted_mean_foreach_var(scen,cellArea)
   global_mean_ref[global_mean_ref == 0 & !global_mean_scen == 0] <- global_mean_scen[global_mean_ref == 0 & !global_mean_scen == 0] # if global mean state in ref period is 0 (e.g. for landuse vars in pnv run?) take the mean scen state instead
   global_mean_ref[global_mean_ref == 0 & global_mean_scen == 0] <- 1 # if both are 0 take 1, then the division is defined but 0 - 0 leads to no change, which is what gamma should show
   norm <- rep(global_mean_ref, each = ncells)
   dim(norm) <- dim(ref)
   s_scen <- scen/norm
   s_ref <- ref/norm

   # normalize both state vectors by the sqrt(amount of state variables) to ensure length(s_ref)==1
   # (this is part of the weighting in the Ostberg code)
   s_ref <- s_ref/sqrt(di[2])
   s_scen <- s_scen/sqrt(di[2])

   return( sqrt(rowSums((s_scen - s_ref)*(s_scen - s_ref))) ) #length of the difference vector s_scen (sl2) - s_ref (sl1) for each cell
}
calcComponent <- function(ref, scen, local, cellArea, export = FALSE) {
   # calc mean ref and scen state
   ref_mean <- apply(ref, c(1,3), mean)
   scen_mean <- apply(scen, c(1,3), mean)
   di <- dim(ref)
   ncells <- di[1]
   nyears <- di[2]
   nvars <- di[3]

   if (local) {
      x <- T_sigmoid_Trafo(state_Diff_local(ref = ref_mean, scen = scen_mean))
   }else{
      x <- T_sigmoid_Trafo(state_Diff_global(ref = ref_mean, scen = scen_mean, cellArea = cellArea))
   }
   # calculation of the change-to-variability ratio in my view is mathematically not correctly described in Heyder and Ostberg
   # the way I understand it: recalculate the c/g/b value for each year of the ref period compared to the mean
   # of the ref period as "scenario" and then calc the variability (sd) of that
   sigma_x_ref_list <- array(0,dim = c(ncells, nyears))
   for (i in 1:nyears) {
      if (local) {
         sigma_x_ref_list[,i] <- T_sigmoid_Trafo(state_Diff_local(ref = ref_mean, scen = ref[,i,])) #todo think about
      }else{
         sigma_x_ref_list[,i] <- T_sigmoid_Trafo(state_Diff_global(ref = ref_mean, scen = ref[,i,], cellArea = cellArea))
      }
   }
   sigma_x_ref <- apply(sigma_x_ref_list, 1, sd)
   if (export) exportVars2GlobalEnv(x,sigma_x_ref)

   return(x*S_change_to_var_ratio(x, sigma_x_ref))
}
balanceShift <- function(ref, scen, epsilon = 10^-4) {
   # param ref with dimension c(ncells,nvars)
   # param scen with dimension c(ncells,nvars)

   # first normalize as for local change
   s_scen <- scen/ref
   s_ref <- array(1, dim = dim(ref))

   # for variables in places, where ref is small (<epsilon), but scen larger
   # (Ostberg code, line 798/vector length calc in line 837)
   # set back scenario vector, to keep the unscaled values (Ostberg code, line 805)
   s_scen[abs(ref) < epsilon & abs(scen) > epsilon] <- scen[abs(ref) < epsilon & abs(scen) > epsilon]

   # for variables in places, where ref and scen are small (<epsilon),
   # set scen to 1 (both are 1, difference is 0 -- no change) (Ostberg code, line 809)
   s_scen[abs(ref) < epsilon & abs(scen) < epsilon] <- 1 # results in no change
   abs_ref <- sqrt(rowSums(s_ref*s_ref)) # absa(_ts) in Sebastians code
   abs_scen <- sqrt(rowSums(s_scen*s_scen)) # absb(_ts) in Sebastians code
   b1 <- 1 - (rowSums(s_ref*s_scen)/abs_ref/abs_scen) #=1-angle_ts

   # these special cases are based on lines 601-614 in Sebastians code, not sure if they are needed
   # tests yield no difference -- disabled so far (FS: 2022-05-23)
   #b1[abs(abs_ref) < epsilon & abs(abs_scen) < epsilon] <- 0 # 1 - 1
   #b1[abs(abs_ref) < epsilon & abs(abs_scen) > epsilon] <- 1 # or 0.5? I do not see the b=b1*2 as described in E.11 in diss in the code
   #b1[abs(abs_ref) > epsilon & abs(abs_scen) < epsilon] <- 1 # similarly

   angle <- acos(1 - b1)*360/2/pi
   angle[b1 == 1] <- 0
   b <- b1*2
   b[angle > 60] <- 1
   return(b)
}
calcEcosystemBalance <- function(ref, scen, export=FALSE) {
  ref_mean <- apply(ref,c(1,3),mean)
  scen_mean <- apply(scen,c(1,3),mean)
  di <- dim(ref)
  ncells <- di[1]
  nyears <- di[2]
  nvars <- di[3]

  b <- balanceShift(ref = ref_mean, scen = scen_mean)
  # calculation of the change-to-variability ratio in my view is mathematically not correctly described in Heyder and Ostberg
  # the way I understand it: recalculate the c/g/b value for each year of the ref period compared to the mean
  # of the ref period as "scenario" and then calc the variability (sd) of that
  sigma_b_ref_list <- array(0,dim = c(ncells, nyears))
  for (i in 1:nyears) {
    sigma_b_ref_list[,i] <- balanceShift(ref = ref_mean, scen = ref[,i,])
  }
  sigma_b_ref <- apply(sigma_b_ref_list, 1, sd)
  if (export) exportVars2GlobalEnv(b, sigma_b_ref)

  return(b*S_change_to_var_ratio(b, sigma_b_ref))
}

#' Create modified gamma data
#'
#' Function to create a modified gamma data file where each cell is compared
#' to the average cell of the precribed reference biome refBiom
#'
#' @param dataFileIn path to input data
#' @param dataFileOut path to save modified data to
#' @param temperature yearly average temperature
#' @param lat
#' @param lon
#' @param refBiom
#'
#' @examples
#' \dontrun{
#' evaluateGamma(dataFileIn = "/media/All/Gamma/Gamma_today_data.RData",
#'       dataFileOut = "/media/All/Gamma/Gamma_preind_compared_to_average_TrBD_data.RData",
#'       temperature = temperature_1716, lat = lat64240, lon = lon64240, refBiom = 4)
#' }
#' @export
evaluateGamma <- function(dataFileIn, dataFileOut, temperature, pft_lai, elevation, lat, lon, refBiom){
  if (dataFileIn == dataFileOut) {
    stop("Same file for input and output of data, would overwrite original data. Aborting.")
  }
  load(dataFileIn)

  # modify vars to compare all to the average cell in ref biom
  data_bc <- list(lat = lat, lon = lon, fpc = apply(fpc_ref, c(1,2), mean),
                  vegc = mean_state_ref[,8],temp = temperature,
                  pft_lai = pft_lai, elevation = elevation)

  biome_classes <- NULL # todo: replace with new cb function classify_biomes(data = data_bc, readOutput = F)
  ref_cells <- which(biome_classes$biome_id == refBiom)

  # first set all scen vars to the ref vars
  state_scen <- state_ref # [1:64240, 1:30, 1:10]
  mean_state_scen <- mean_state_ref
  fpc_scen <- fpc_ref #
  bft_scen <- bft_ref
  cft_scen <- cft_ref

  di_state <- dim(state_scen)
  di_fpc <- dim(fpc_scen)
  di_bft <- dim(bft_scen)
  di_cft <- dim(cft_scen)
  # now replace all ref cells with that of the mean ref biom cell
  av_year_state <- apply(state_scen[ref_cells,,],c(2,3),mean) #FS 2022-08-10: keeping the year-to-year variation
  state_ref <- rep(av_year_state, each = di_state[1]) #FS 2022-08-10: keeping the year-to-year variation
  dim(state_ref) <- di_state
  mean_state_ref <- rep(colMeans(av_year_state), each = di_state[1]) # is the same for each year, thus for the mean just take one year
  fpc_ref <- rep(apply(fpc_scen[ref_cells,,],c(2,3),mean), each = di_fpc[1])
  dim(fpc_ref) <- di_fpc
  bft_ref <- rep(apply(bft_scen[ref_cells,,],c(2,3),mean), each = di_bft[1])
  dim(bft_ref) <- di_bft
  cft_ref <- rep(apply(cft_scen[ref_cells,,],c(2,3),mean), each = di_cft[1])
  dim(cft_ref) <- di_cft

  # and write out the modified data
  save(state_ref,mean_state_ref,state_scen,mean_state_scen,fpc_ref,fpc_scen,bft_ref,bft_scen,cft_ref,cft_scen,lat,lon,cellArea,file = dataFileOut)

}

#' Create modified gamma data for crosstable
#'
#' Function to create a modified gamma data file where for each present day biome
#' the average cell is compared to the average cell of all other biomes. This
#' can then be used to compute a crosstable with the average difference between
#' each of them as in the SI of Ostberg et al. 2013
#' (Critical impacts of global warming on land ecosystems)
#'
#' @param dataFileIn path to input data
#' @param dataFileOut path to save modified data to
#' @param temperature yearly average temperature
#' @param pft_lai yearly pft_lai values for the biome classification
#' @param elevation elevation input file for the biome classification
#' @param lat
#' @param lon
#'
#' @examples
#' \dontrun{
#' gammaCrossTable(dataFileIn = "/media/All/Gamma/Gamma_today_data.RData",
#'       dataFileOut = "/media/All/Gamma/Gamma_preind_compared_to_average_TrBD_data.RData",
#'       temperature = temperature_1716, pft_lai = pft_lai_1996, elevation = elevation, lat = lat64240, lon = lon64240, refBiom = 4)
#' }
#' @export
gammaCrossTable <- function(dataFileIn, dataFileOut, temperature, pft_lai, elevation, pickCells = NULL, lat, lon){
  if (dataFileIn == dataFileOut) {
    stop("Same file for input and output of data, would overwrite original data. Aborting.")
  }
  load(dataFileIn)

  # modify vars to compare all to the average cell in ref biom
  data_bc <- list(lat = lat, lon = lon, fpc = apply(fpc_ref, c(1,2), mean),
                  vegc = mean_state_ref[,8],temp = temperature,
                  pft_lai = pft_lai, elevation = elevation)
  biome_classes <- NULL # todo: replace with new cb function classify_biomes(data = data_bc, readOutput = F)

  #save scenario state vectors, they contain relevant data (ref can go)
  state_scen_sav <- state_scen
  fpc_scen_sav <- fpc_scen
  bft_scen_sav <- bft_scen
  cft_scen_sav <- cft_scen


  nbiomes <- max(biome_classes$biome_id) # by default 19
  state_ref <- array(0,dim = c(nbiomes,nbiomes,dim(state_scen_sav)[2:3]))
  state_scen <- state_ref
  fpc_ref <- array(0,dim = c(nbiomes,nbiomes,dim(fpc_scen_sav)[2:3]))
  fpc_scen <- fpc_ref
  bft_ref <- array(0,dim = c(nbiomes,nbiomes,dim(bft_scen_sav)[2:3]))
  bft_scen <- bft_ref
  cft_ref <- array(0,dim = c(nbiomes,nbiomes,dim(cft_scen_sav)[2:3]))
  cft_scen <- cft_ref

  # now replace all ref cells with that of the mean ref biome cell
  for (b in sort(unique(biome_classes$biome_id))) {
    ref_cells <- which(biome_classes$biome_id == b)

    if (is.null(pickCells)) {
      av_state <- apply(state_scen_sav[ref_cells,,],c(2,3),mean) # average over cells, keeping the average year-to-year variation
      av_fpc <- apply(fpc_scen_sav[ref_cells,,],c(2,3),mean)
      av_bft <- apply(bft_scen_sav[ref_cells,,],c(2,3),mean)
      av_cft <- apply(cft_scen_sav[ref_cells,,],c(2,3),mean)
    }else{
      av_state <- state_scen_sav[pickCells[b],,] # average over cells, keeping the average year-to-year variation
      av_fpc <- fpc_scen_sav[pickCells[b],,]
      av_bft <- bft_scen_sav[pickCells[b],,]
      av_cft <- cft_scen_sav[pickCells[b],,]
    }

    state_ref[b,,,] <- rep(av_state, each = nbiomes)
    state_scen[,b,,] <- rep(av_state, each = nbiomes)

    mean_state_ref <- apply(state_ref,c(1,3),mean)
    mean_state_scen <- apply(state_scen,c(1,3),mean)

    fpc_ref[b,,,] <- rep(av_fpc, each = nbiomes)
    fpc_scen[,b,,] <- rep(av_fpc, each = nbiomes)

    bft_ref[b,,,] <- rep(av_bft, each = nbiomes)
    bft_scen[,b,,] <- rep(av_bft, each = nbiomes)

    cft_ref[b,,,] <- rep(av_cft, each = nbiomes)
    cft_scen[,b,,] <- rep(av_cft, each = nbiomes)

  }
  dim(state_ref) <- c(nbiomes*nbiomes,dim(state_scen_sav)[2:3])
  dim(state_scen) <- c(nbiomes*nbiomes,dim(state_scen_sav)[2:3])
  dim(fpc_ref) <- c(nbiomes*nbiomes,dim(fpc_scen_sav)[2:3])
  dim(fpc_scen) <- c(nbiomes*nbiomes,dim(fpc_scen_sav)[2:3])
  dim(bft_ref) <- c(nbiomes*nbiomes,dim(bft_scen_sav)[2:3])
  dim(bft_scen) <- c(nbiomes*nbiomes,dim(bft_scen_sav)[2:3])
  dim(cft_ref) <- c(nbiomes*nbiomes,dim(cft_scen_sav)[2:3])
  dim(cft_scen) <- c(nbiomes*nbiomes,dim(cft_scen_sav)[2:3])

  lat <- rep(0,nbiomes*nbiomes)
  lon <- rep(1,nbiomes*nbiomes)
  cellArea <- calcCellarea(lat = lat)
  # and write out the modified data
  save(state_ref,mean_state_ref,state_scen,mean_state_scen,fpc_ref,fpc_scen,bft_ref,bft_scen,cft_ref,cft_scen,lat,lon,cellArea,file = dataFileOut)

}
################# biome (dis-)aggregation functions ##################
#' Returns the average gamma values across regions
#'
#' Returns the average value across either 4 regions or 19 biomes for gamma
#' and each of the subcomponents for each
#'
#' @param gamma gamma data object with dimension c(regions/biomes[4/16],gamma_components[10/8],minmeanmax/quantiles[3])
#' @param biome_class path to save modified data to
#' @param type return minimum,mean,maximum or Q10,Q50,Q90
#' @param classes into which regions should be disaggregated? "4biomes" or "allbiomes" default "4biomes"
#'
#' @examples
#' \dontrun{
#' disaggregateGammaIntoBiomes(gamma = gamma_historic_nitrogen_classic_bnf,
#'                    biome_class = gamma_historic_withNitrogen_biome_classes,
#'                    type = "quantile",classes = "4biomes")
#' }
#' @export
disaggregateGammaIntoBiomes <- function(gamma, biome_class, type = "quantile", classes = "4biomes") {
  if (classes == "4biomes") {
    tropics <- c(1,2,9,10,11)
    temperate <- c(3,4,5,6,12,13,14)
    boreal <- c(7,8)
    arctic <- c(15,16)
    cell_list <- list(tropical_cells = which(biome_class$biome_id %in% tropics),
                      temperate_cells = which(biome_class$biome_id %in% temperate),
                      boreal_cells = which(biome_class$biome_id %in% boreal),
                      arctic_cells = which(biome_class$biome_id %in% arctic) )
    nclasses <- 4
  }else if (classes == "allbiomes") {
    nclasses <- max(unique(biome_class$biome_id))
  }else{
    stop(paste0("Unknown parameter classes: ",classes,", should be either '4biomes' or 'allbiomes'"))
  }


  gamma_dims <- length(gamma)
  gamma_biomes <- array(0,dim = c(nclasses,gamma_dims,3)) # c(biome,gamma_components,min/median/max)
  if (classes == "4biomes") { # aggregate to trop/temp/boreal/arctic
    for (b in 1:nclasses) {
      for (c in 1:gamma_dims) {
        if (type == "minmeanmax") {
          gamma_biomes[b,c,] <- c( min(gamma[[c]][cell_list[[b]]],na.rm = T),
                                   mean(gamma[[c]][cell_list[[b]]],na.rm = T),
                                   max(gamma[[c]][cell_list[[b]]],na.rm = T) )
        }else if (type == "quantile") {
          gamma_biomes[b,c,] <- c( quantile(gamma[[c]][cell_list[[b]]], probs = c(0.1,0.5,0.9), na.rm = T) )
        }else{stop(paste("type",type,"unknown. please choose either 'quantile' or 'minmeanmax'"))}
      }
    }
  }else if (classes == "allbiomes") { #calculate all biomes separately
    for (b in 1:nclasses) {
      for (c in 1:gamma_dims) {
        if (type == "minmeanmax") {
          gamma_biomes[b,c,] <- c( min(gamma[[c]][which(biome_class$biome_id == b)],na.rm = T),
                                   mean(gamma[[c]][which(biome_class$biome_id == b)],na.rm = T),
                                   max(gamma[[c]][which(biome_class$biome_id == b)],na.rm = T) )
        }else if (type == "quantile") {
          gamma_biomes[b,c,] <- c( quantile(gamma[[c]][which(biome_class$biome_id == b)], probs = c(0.1,0.5,0.9), na.rm = T) )
        }else{stop(paste("type",type,"unknown. please choose either 'quantile' or 'minmeanmax'"))}
      }
    }
  }else{
    stop(paste0("Unknown parameter classes: ",classes,", should be either '4biomes' or 'allbiomes'"))
  }

  return(gamma_biomes)
}
################# gamma plotting functions ##################
#' Plot gamma map to screen
#'
#' Function to plot a global map of gamma values [0-1] per grid cell to screen
#'
#' @param data folder of reference run
#' @param focusBiome highlight the biome with this id and desaturate all other (default NULL -- no highlight)
#' @param biome_classes to mask the focusBiome from
#' @param title character string title for plot
#' @param legendtitle character string legend title
#' @param legYes logical. whether to plot legend or not. defaults to TRUE
#' @param legScale scaling factor for legend. defaults to 1
#'
#' @return None
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
plotGammaMapToScreen <- function(data,focusBiome = NULL, biome_classes = NULL,
                                 title, legendtitle, titleSize = 1, legYes = T) {
   brks <- seq(0,1,0.1)
   data[data < brks[1]] <- brks[1]
   data[data > brks[length(brks)]] <- brks[length(brks)]
   palette <- c("white","steelblue1","royalblue",RColorBrewer::brewer.pal(7,"YlOrRd"))
   if (!is.null(focusBiome)){
     focus <- data
     focus[!(biome_classes == focusBiome)] <- NA
     palette_lowSat <- adjustcolor(palette,alpha.f = 0.25)
     ra_f <- raster::raster(ncols = 720, nrows = 360)
     ra_f[raster::cellFromXY(ra_f,cbind(lon,lat))] <-  focus
   }
   ra <- raster::raster(ncols = 720, nrows = 360)
   ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  data
   range <- range(data)
   extent <- raster::extent(c(-180, 180, -60, 90))
   par(mar = c(0,0,1,3), oma = c(0,0,0,0),bty = "n")
   if (is.null(focusBiome)) {
     raster::plot(ra, ext = extent, breaks = brks, col = palette, main = "", legend = FALSE, axes = FALSE)
   }else{
     raster::plot(ra, ext = extent, breaks = brks, col = palette_lowSat, main = "", legend = FALSE, axes = FALSE)
     raster::plot(ra_f, ext = extent, breaks = brks, col = palette, main = "", legend = FALSE, axes = FALSE,add = T)
   }
   title(main = title, line = -2, cex.main = titleSize)
   maps::map('world', add = TRUE, res = 0.4, lwd = 0.25, ylim = c(-60,90))
   if (legYes) {
      fields::image.plot(legend.only = TRUE, zlim = range(brks), col = palette, breaks = brks, lab.breaks = brks, legend.shrink = 0.7,
                         legend.args = list(legendtitle, side = 3, font = 2, line = 1))
   }
}
#' Plot gamma map to file
#'
#' Function to plot a global map of gamma values [0-1] per grid cell to file
#'
#' @param data folder of reference run
#' @param file to write into
#' @param focusBiome highlight the biome with this id and desaturate all other (default NULL -- no highlight)
#' @param biome_classes to mask the focusBiome from
#' @param title character string title for plot
#' @param legendtitle character string legend title
#' @param eps write as eps or png
#' @param legYes logical. whether to plot legend or not. defaults to TRUE
#' @param legScale scaling factor for legend. defaults to 1
#'
#' @return None
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
plotGammaMap <- function(data, file, focusBiome = NULL, biome_classes = NULL,
                         title, legendtitle, eps = FALSE, titleSize = 1, legYes = T) {
   if (eps) {
      file <- strsplit(file,".",fixed = TRUE)[[1]]
      file <- paste(c(file[1:(length(file) - 1)],"eps"),collapse = ".")
      ps.options(family = c("Helvetica"), pointsize = 18)
      postscript(file, horizontal = FALSE, onefile = FALSE, width = 22, height = 8.5, paper = "special")
   }else{
      png(file, width = 7.25, height = 3.5, units = "in", res = 300, pointsize = 6,type = "cairo")
   }
   plotGammaMapToScreen(data = data, focusBiome = focusBiome, biome_classes = biome_classes, title = title, legendtitle = legendtitle, titleSize = titleSize, legYes = legYes)
   dev.off()
}
#' Plot radial gamma plot to screen
#'
#' Function to plot an aggregated radial status of gamma values [0-1]
#' for the different sub-categories to screen
#'
#' @param data gamma data array c([ngammacomponents],3[min,median,max])
#' @param title character string title for plot
#' @param zoom scaling factor for circle plot. defaults to 1
#' @param type plot type, 'legend1' for variable and color legend,
#'             'legend2' for value legend, or 'regular' (default setting)
#'             for the regular gamma plot
#' @param titleSize scaling factor for tile. defaults to 1
#'
#' @return None
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
plotGammaRadialToScreen <- function(data, title, zoom = 1.0, type = "regular", titleSize = 2, titleline=-2,quantile=T) {
   suppressPackageStartupMessages(require(circlize))
   require(RColorBrewer)
   gamma_dims <- length(data[,1])
   if (gamma_dims == 10) {
     names <- c( gamma = "m-eco", deltav = "vegetation\nstructure",
                 local = "local\nchange", global = "global\nimportance",
                 balance =  "ecosystem\nbalance", cstocks = "carbon stocks", cfluxes = "carbon fluxes",
                 wfluxes = "water fluxes", nstocks = "nitrogen\nstocks", nfluxes = "nitrogen fluxes")
     set <- brewer.pal(12, "Set3") #c(blue-green, yellow,violet,red,blue,orange,green,pink,grey,purple,green-blue,yellow-orange)
     colz <- set[c(4,7,8,11,1,3,10,5,12,6)]
     #set <- brewer.pal(9, "Set1") #c(red,blue,green,purple,orange,yellow,brown,pink,grey)
     #colz <- c("limegreen", "darkgreen", "maroon","orchid4","bisque4",
     #         "orangered","sienna",brewer.pal(6, "PuBu")[6], "yellow" , "orange")
     #                  gamma    deltaV      lc       gc      eb          cs      cf       wf       ns      nf
     angles <- matrix(c(90,270, 216,252,  180,216, 144,180, 108,144,  -18,18, -54,-18, -90,-54,  54,90, 18,54 ),byrow = T,nrow = length(colz))
   }else if (gamma_dims == 8) {
     names <- c( gamma = expression(paste(Gamma, " metric")), deltav = "vegetation\nstructure",
                 local = "local\nchange", global = "global\nimportance",
                 balance =  "ecosystem\nbalance", cstocks = "carbon stocks",
                 cfluxes = "carbon fluxes", wfluxes = "water fluxes")
     colz <- c("darkgoldenrod", brewer.pal(5,"Greens")[5], brewer.pal(6, "Set1")[seq(2, 6, by = 2)],
               rev(brewer.pal(6, "Oranges")[c(4,5)]), brewer.pal(6, "PuBu")[6])
     angles <- matrix(c(234,270, 198,234, 162,198, 126,162, 90,126,     18,54, -18,18, -54,-18),byrow = T,nrow = length(colz))
   }else{
     stop(paste("Unknown number of dimensions for gamma data:",gamma_dims))
   }
   par(oma = c(0,0,0,0), mar = c(0,0,0,0))
   plot(c(-zoom, zoom), c(-zoom, zoom), type = "n", axes = FALSE, ann = FALSE, asp = 1, main = "")
   title(main = title, line = titleline, cex.main = titleSize)
   if (type == "legend1") {
      draw.sector(0, 360, rou1 = 1)
     ro = c(1,1.1,0.8,1.1,0.8, 1,1,1,1,1)
      for (i in 1:length(angles[,1])) {
         draw.sector(start.degree = angles[i,1] + 90, end.degree = angles[i,2] + 90, col = colz[i], clock.wise = F, rou1 = 0, rou2 = ro[i],border = "black")
      }
     if (gamma_dims == 10) {
       text(names,x = c(1.1,1.0,0.2,-0.8,-1.6, -0.4,0.7,1.05,   -1.7,-1.5),y = c(-0.15,-0.9,-1.3,-1.3,-0.9, 1.2,1,0.25,  0.3,1), adj = 0)
     }else if (gamma_dims == 8) {
       text(names,x = c(1.1,0.6,-0.2,-1.2,-1.7, -1.5,-0.4,0.7),y = c(-0.3,-1.1,-1.3,-1,-0.5, 1,1.2,1),adj = 0)
     }else{
       stop(paste("Unknown number of dimensions for gamma data:",gamma_dims))
     }
     draw.sector(start.degree = (angles[3,1]+angles[3,2])/2+90, end.degree = (angles[3,1]+angles[3,2])/2+90, rou1 = 0.7, rou2 = 1.1)# line lc
     draw.sector(start.degree = -9, end.degree = -9, rou1 = 0.9, rou2 = 1.05)# line m-eco
     draw.sector(start.degree = (angles[5,1]+angles[5,2])/2+90, end.degree = (angles[5,1]+angles[5,2])/2+90, rou1 = 0.7, rou2 = 1.1)# line eb
     draw.sector(start.degree = 180, end.degree = 180, clock.wise = F, rou1 = -1.2, rou2 = 1.2, border = "black",lwd = 2)
   }else if (type == "legend2") {
      text("+",x = 0, y = 0)
      draw.sector(0, 360, rou1 = 1)
      draw.sector(0, 360, rou1 = 0.65)
      draw.sector(0, 360, rou1 = 0.3)
      draw.sector(start.degree = -18 + 90, end.degree = 18 + 90, clock.wise = F, rou1 = 0.55, border = "black") #sector
      draw.sector(start.degree = 90, end.degree = 90, clock.wise = F, rou1 = 0.4, rou2 = 0.8, border = "black") #uncertainty arrow
      draw.sector(start.degree = -9 + 90, end.degree = 9 + 90, clock.wise = F, rou1 = 0.8, rou2 = 0.8, border = "black") #uncertainty lower
      draw.sector(start.degree = -9 + 90, end.degree = 9 + 90, clock.wise = F, rou1 = 0.4, rou2 = 0.4, border = "black") #uncertainty upper
      draw.sector(start.degree = 270 - 270, end.degree = 270 - 270, clock.wise = F, rou1 = 0.3, rou2 = 1.3, border = "black",lty = "dashed") #0.3
      draw.sector(start.degree = 280 - 270, end.degree = 280 - 270, clock.wise = F, rou1 = 0.65, rou2 = 1.3, border = "black",lty = "dashed") #0.65
      draw.sector(start.degree = 290 - 270, end.degree = 290 - 270, clock.wise = F, rou1 = 1, rou2 = 1.3, border = "black",lty = "dashed") #1.0
      text(c("0.3","0.65","1"),x = c(1.4,1.45,1.25), y = c(0,0.25,0.45))
      # plot how the whiskers are calculated
      if (quantile) text(c("Q90","Q50","Q10"),x = c(-0.3,-0.29,-0.26), y = c(0.8,0.48,0.35),cex=0.7) # quantile case
      else text(c("max","mean","min"),x = c(-0.3,-0.29,-0.26), y = c(0.8,0.48,0.35),cex=0.7) # minmeanmax case
   }else if (type == "regular") {
     draw.sector(180, 360, rou1 = 1,col = "gray80")
      for (i in 1:length(angles[,1])) {
         mangle <- mean(angles[i,])
         if (i==1) mangle <- -98
         dmin <- data[i,1]
         dmedian <- data[i,2]
         dmax <- data[i,3]
         draw.sector(start.degree = angles[i,1] + 90, end.degree = angles[i,2] + 90, col = colz[i], rou1 = dmedian,clock.wise = F, border = "black")
         draw.sector(start.degree = mangle + 90, end.degree = mangle + 90, clock.wise = F, rou1 = dmin, rou2 = dmax, border = "black") #uncertainty arrow
         draw.sector(start.degree = mangle - 9 + 90, end.degree = mangle + 9 + 90, clock.wise = F, rou1 = dmin, rou2 = dmin, border = "black") #uncertainty lower
         draw.sector(start.degree = mangle - 9 + 90, end.degree = mangle + 9 + 90, clock.wise = F, rou1 = dmax, rou2 = dmax, border = "black") #uncertainty upper
      }
      draw.sector(0, 360, rou1 = 1)
      draw.sector(0, 360, rou1 = 0.6)
      draw.sector(0, 360, rou1 = 0.3)
      draw.sector(start.degree = 180, end.degree = 180, clock.wise = F, rou1 = -1.2, rou2 = 1.2, border = "black",lwd = 2)

   }else {
      stop(paste0("Unknown type ",type,". Please use 'legend1' for variable and color legend,
                  'legend2' for value legend, or 'regular' (default setting) for the regular gamma plot."))
   }
}

#' Plot radial gamma plot to file
#'
#' Function to plot an aggregated radial status of gamma values [0-1]
#' for the different sub-categories to file
#'
#' @param data gamma data array c([ngammacomponents],3[min,median,max])
#' @param file to write into
#' @param title character string title for plot
#' @param type plot type, 'legend1' for variable and color legend,
#'             'legend2' for value legend, or 'regular' (default setting)
#'             for the regular gamma plot
#' @param eps write as eps or png
#' @param legYes logical. whether to plot legend or not. defaults to TRUE
#'
#' @return None
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
plotGammaRadial <- function(data, file, title, legYes = T, eps = FALSE,quantile=T) {
   #param data gamma data array c(8[ncomponents],3[min,median,max])
   #param title title for plot
   #param type plot type, 'legend1' for variable and color legend, 'legend2' for value legend, or 'regular' (default setting) for the regular gamma plot

   if (length(which(data < 0 | data > 1)) > 0) print("Warning: there are values in data outside the expected gamma range [0..1].")
   if (eps) {
      file <- strsplit(file,".",fixed = TRUE)[[1]]
      file <- paste(c(file[1:(length(file) - 1)],"eps"), collapse = ".")
      ps.options(family = c("Helvetica"), pointsize = 18)
      postscript(file, horizontal = FALSE, onefile = FALSE, width = 22, height = 8.5, paper = "special")
   }else {
      png(file, width = 7.25, height = 3.5, units = "in", res = 300, pointsize = 6,type = "cairo")
   }
   # adjust the margins, dependent on whether a legend should be plotted or not
   par(fig = c(0,0.7,0,1))#, oma=c(0,0,0,0),mar=c(0,0,0,0))

   #plot main gamma radial
   plotGammaRadialToScreen(data = data, title = title,zoom=1.0, type = "regular")

   if (legYes) {
      par(fig = c(0.7,1,0,0.5), new = TRUE)#, oma=c(0,0,0,0),mar=c(0,0,0,0))
      plotGammaRadialToScreen(data = data, title = "", zoom = 1.5, type = "legend1")
      par(fig = c(0.7,1,0.5,1), new = TRUE)#, oma=c(0,0,0,0),mar=c(0,0,0,0))
      plotGammaRadialToScreen(data = data, title = "",zoom = 1.5, type = "legend2", quantile=quantile)
   }
   dev.off()
}

#' Plot radial gamma plot to file with 4/16 biomes
#'
#' Function to plot an aggregated radial status of gamma values [0-1]
#' for the different sub-categories to file
#'
#' @param data gamma data array c(4[biomes],[ngammacomponents],3[min,median,max])
#' @param file to write into
#' @param title character string title for plot
#' @param eps write as eps or png
#'
#' @return None
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
plotGammaRadial4 <- function(data, biomeNames, file, eps = FALSE, quantile=T) {
  if (length(which(data < 0 | data > 1)) > 0) print("Warning: there are values in data outside the expected gamma range [0..1].")
  if (eps) {
    file <- strsplit(file,".",fixed = TRUE)[[1]]
    file <- paste(c(file[1:(length(file) - 1)],"eps"), collapse = ".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file, horizontal = FALSE, onefile = FALSE, width = 15, height = 10, paper = "special")
  }else {
    png(file, width = 5.25, height = 3.5, units = "in", res = 300, pointsize = 6,type = "cairo")
  }
  d <- length(data[,1,1])
  if (d==16 | d==4){
    k <- sqrt(d)
    xs <- seq(0,0.6,length.out = k+1)
    ys <- seq(0.98,0,length.out = k+1)
    for (x in 1:k){
      for (y in 1:k){
        if (x==1 & y==1) par(fig = c(xs[x],xs[x+1],ys[y+1],ys[y]), xpd = T, oma=c(0,0,0,0),mar=c(0,0,0,0))
        else             par(fig = c(xs[x],xs[x+1],ys[y+1],ys[y]), xpd = T, new=T)
        plotGammaRadialToScreen(data = data[(x-1)*k+y,,], title = "",zoom=1.0, type = "regular")
        mtext(text = biomeNames[(x-1)*k+y], side = 3, line = -0.5, cex = 1,font = 2)
      }
    }
  }else{
    stop(paste("Unknown number of biomes: ",length(data[,1,1])))
  }
  #legend
  par(fig = c(0.6,1,0.1,0.6), new = TRUE)#, oma=c(0,0,0,0),mar=c(0,0,0,0))
  plotGammaRadialToScreen(data = data[1,,], title = "", zoom = 1.5, type = "legend1")
  par(fig = c(0.6,1,0.5,1.0), new = TRUE)#, oma=c(0,0,0,0),mar=c(0,0,0,0))
  plotGammaRadialToScreen(data = data[1,,], title = "legend",zoom = 1.5,
                          type = "legend2",titleSize = 1,quantile = quantile)
  dev.off()
}

#' Plot radial gamma plot to file with 4/16 biomes
#'
#' Function to plot an aggregated radial status of gamma values [0-1]
#' for the different sub-categories to file
#'
#' @param data gamma data array c(4[biomes],[ngammacomponents],3[min,median,max])
#' @param file to write into
#' @param title character string title for plot
#' @param eps write as eps or png
#'
#' @return None
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
plotBiomesToScreen <- function(biome_ids, biome_class_names, title,
                               titleSize = 2, legYes = T, leg_scale = 0.5) {
   require(raster)
   require(RColorBrewer)
   #---- setting up colors and biome names ----------------------------------------------------------------#
   colz <- c(rev(brewer.pal(6, "YlOrBr")), # warm
                         rev(brewer.pal(9, "YlGn")[c(3,5,7,9)]),
                         rev(brewer.pal(9, "GnBu"))[c(2:4,6,8,9)], # cold below forest
                         "white", #"lightblue" # Water
                         "lightgrey", #Rocks & Ice
                         "magenta" #montane Tundra/Grassland
   )
   biome_class_cols <-  colz[c(1,2,7,8,9,10,13,12,3,4,5,14,15,16,19,11,6,17,18)]
   if (!(length(biome_class_names) == length(biome_class_cols))) stop("Size of biome class names and colors do not match -- should be 18.")
   #---- plotting ----------------------------------------------------------------#
   brks <- seq(min(biome_ids) - 0.5, max(biome_ids) + 0.5, 1)
   ra <- raster::raster(ncols = 720, nrows = 360)
   range <- range(biome_ids)
   ra[raster::cellFromXY(ra,cbind(lon,lat))] <-  biome_ids
   extent <- raster::extent(c(-180, 180, -60, 90))
   par(mar = c(0,0,0,0), oma = c(0,0,0,0),bty = "n")
   raster::plot(ra, ext = extent, breaks = brks, col = biome_class_cols, main = "", legend = FALSE, axes = FALSE)
   title(main = title, line = -2, cex.main = titleSize)
   if (legYes) {
      legend(x = -180, y = 27, legend = biome_class_names, fill = biome_class_cols, col = biome_class_cols, cex = leg_scale, bg = "white", bty = "o")
   }
   maps::map('world', add = TRUE, res = 0.4, lwd = 0.25, ylim = c(-60,90))
}

#' Plot radial gamma plot to file with 4/16 biomes
#'
#' Function to plot an aggregated radial status of gamma values [0-1]
#' for the different sub-categories to file
#'
#' @param data gamma data array c(4[biomes],[ngammacomponents],3[min,median,max])
#' @param file to write into
#' @param title character string title for plot
#' @param eps write as eps or png
#'
#' @return None
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
plotBiomes <- function(biome_ids, biome_class_names, file, title, titleSize = 2,
                       legYes = T, leg_scale = 1, eps = FALSE) {
   #---- plotting ----------------------------------------------------------------#
   if (eps) {
      file <- strsplit(file,".",fixed = TRUE)[[1]]
      file <- paste(c(file[1:(length(file) - 1)],"eps"),collapse = ".")
      ps.options(family = c("Helvetica"), pointsize = 18)
      postscript(file, horizontal = FALSE, onefile = FALSE, width = 22, height = 8.5, paper = "special")
   }else{
      png(file, width = 7.25, height = 3.5, units = "in", res = 300, pointsize = 6,type = "cairo")
   }
   plotBiomesToScreen(biome_ids = biome_ids, biome_class_names = biome_class_names, title = title, titleSize = titleSize, legYes = legYes, leg_scale = leg_scale)
   dev.off()
}

#' Plot radial gamma plot to file with 4/16 biomes
#'
#' Function to plot an aggregated radial status of gamma values [0-1]
#' for the different sub-categories to file
#'
#' @param data input data with dimension c(nbiome_classes,3) -- Q10,Q50,Q90 each
#' @param biome_class_names to write into
#' @param title character string title for plot
#' @param titleSize character string title for plot
#' @param leg_scale character string title for plot
#'
#' @return None
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
plotBiomesAveragesToScreen <- function(data, biome_class_names, title,
                               titleSize = 2, leg_scale = 0.5) {
  require(raster)
  require(RColorBrewer)
  #---- setting up colors and biome names ----------------------------------------------------------------#
  brks <- seq(0,1,0.1)
  data[data < brks[1]] <- brks[1]
  data[data > brks[length(brks)]] <- brks[length(brks)]
  palette <- c("white","steelblue1","royalblue",RColorBrewer::brewer.pal(7,"YlOrRd"))
  colIndex <- floor(data[,2]*10) + 1
  if (!(length(biome_class_names) == dim(data)[1])) stop("Size of biome class names and data input do not match.")
  #---- plotting ----------------------------------------------------------------#
  plot(NA, xlim = c(0,1), ylim = c(0,1), main = title, axes = F,cex.main = titleSize,xlab = "", ylab = "")
  legend(x = 0, y = 1, legend = biome_class_names, fill = palette[colIndex], col = palette[colIndex],border = palette[colIndex], cex = leg_scale, bg = "white", bty = "o")
}

#' Plot radial gamma plot to file with 4/16 biomes
#'
#' Function to plot an aggregated radial status of gamma values [0-1]
#' for the different sub-categories to file
#'
#' @param data gamma data array c(4[biomes],[ngammacomponents],3[min,median,max])
#' @param file to write into
#' @param biome_class_names to write into
#' @param title character string title for plot
#' @param titleSize character string title for plot
#' @param leg_scale character string title for plot
#' @param eps write as eps or png
#'
#' @return None
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
plotBiomesAverage <- function(data, file, biome_class_names, title, titleSize = 2,
                      leg_scale = 1, eps = FALSE) {
  #---- plotting ----------------------------------------------------------------#
  if (eps) {
    file <- strsplit(file,".",fixed = TRUE)[[1]]
    file <- paste(c(file[1:(length(file) - 1)],"eps"),collapse = ".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file, horizontal = FALSE, onefile = FALSE, width = 22, height = 8.5, paper = "special")
  }else{
    png(file, width = 4, height = 3, units = "in", res = 300, pointsize = 6,type = "cairo")
  }
  plotBiomesAveragesToScreen(data = data, biome_class_names = biome_class_names, title = title, titleSize = titleSize, leg_scale = leg_scale)
  dev.off()
}

#' Plot crosstable showing (dis-)similarity between average biome pixels
#'
#' Function to plot a crosstable showing (dis-)similarity between average
#' biome pixels based on M-ECO (former Gamma) metric from LPJmL simulations
#'
#' @param data crosstable data as array with [nbiomes,nbiomes] and row/colnames
#' @param lmar left margin for plot in lines (default: 3)
#'
#' @return None
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
plotGammaCrossTableToScreen <- function(data, lmar = 3) {
  #data prep
  data <- round(data,digits = 2)
  x = 1:ncol(data)
  y = 1:nrow(data)
  centers <- expand.grid(y,x)
  #coloring
  palette <- c("white","steelblue1","royalblue",RColorBrewer::brewer.pal(7,"YlOrRd"))
  brks <- seq(0,1,0.1)
  #color.picker <- function(z){
  #  if(is.na(z)){return("black")}
  #  else {return(palette[floor(z*10)+1])}
  #}

  # plot margins
  par(mar = c(0,lmar,2,0)) #bltr

  image(x, y, t(data),
        col = palette,
        breaks = brks,
        xaxt = 'n',
        yaxt = 'n',
        xlab = '',
        ylab = '',
        ylim = c(max(y) + 0.5, min(y) - 0.5)
  )
  text(centers[,2], centers[,1], c(data), col = "black")

  #add margin text
  mtext(attributes(data)$dimnames[[2]], at=1:ncol(data), padj = -1)
  mtext(attributes(data)$dimnames[[1]], at=1:nrow(data), side = 2, las = 1, adj = 1,line = 1)

  #add black lines
  abline(h=y + 0.5)
  abline(v=x + 0.5)
}

#' Plot crosstable to file showing (dis-)similarity between average biome pixels
#'
#' Function to plot to file a crosstable showing (dis-)similarity between average
#' biome pixels based on M-ECO (former Gamma) metric from LPJmL simulations
#'
#' @param data crosstable data as array with [nbiomes,nbiomes] and row/colnames
#' @param file to write into
#' @param lmar left margin for plot in lines (default: 3)
#' @param eps write as eps or png
#'
#' @return None
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
plotGammaCrossTable <- function(data, file, lmar=3, eps = FALSE) {
  #---- plotting ----------------------------------------------------------------#
  if (eps) {
    file <- strsplit(file,".",fixed = TRUE)[[1]]
    file <- paste(c(file[1:(length(file) - 1)],"eps"),collapse = ".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file, horizontal = FALSE, onefile = FALSE, width = 22, height = 8.5, paper = "special")
  }else{
    png(file, width = 6, height = 3, units = "in", res = 300, pointsize = 6,type = "cairo")
  }
  plotGammaCrossTableToScreen(data = data, lmar = lmar)
  dev.off()
}
