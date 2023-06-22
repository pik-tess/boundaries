#' Calculate biosphere status based on MCOL from a PNV run and LU run of LPJmL
#'
#' Calculate biosphere status based on MCOL from a PNV run and LU run of LPJmL
#' @param files_scenario list with variable names and corresponding file paths
#' (character string) of the scenario LPJmL run. All needed files are
#' provided in XXX. E.g.: list(leaching = "/temp/leaching.bin.json")
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the reference LPJmL run. All needed files are
#' provided in XXX. E.g.: list(leaching = "/temp/leaching.bin.json"). If not
#' needed for the applied method, set to NULL.
#' @param time_span_scenario time span to be used for the scenario run, defined
#' as a character string, e.g. `as.character(1982:2011)` (default)
#' @param time_span_reference time span to be used for the scenario run, defined
#' as an integer vector, e.g. `as.character(1901:1930)`. Can differ in offset
#' and length from `time_span_scenario`! If `NULL` value of `time_span_scenario`
#' is used
#' @param reference_npp_time_span time span to read reference npp from, using
#'     index years 10:39 from potential npp input if set to NULL (default: NULL)
#' @param reference_npp_file file to read reference npp from, using
#'        potential npp input if set to NULL (default: NULL)
#' @param gridbased logical are pft outputs gridbased or pft-based?
#' @param readPreviouslySavedData flag whether to read previously saved data
#'        instead of reading it in from output files (default FALSE)
#' @param saveDataFile whether to save input data to file (default FALSE)
#' @param dataFile file to save/read input data to/from (default NULL)
#' @param include_fire boolean include firec in calculation of MCOL? (default T)
#' @param external_fire instead of reading in firec for fire emissions, read in
#'        this external firec file from a separate spitfire run with disabled
#'        lighning. this will then include only human induced fires (default F)
#' @param external_wood_harvest include external wood harvest from LUH2_v2h
#'        (default F)
#' @param grass_scaling whether to scale pasture harvest according to
#'        data given via grass_harvest_file (default F)
#' @param npp_threshold lower threshold for npp (to mask out non-lu areas
#'        according to Haberl et al. 2007). Below MCOL will be set to 0.
#'        (default: 20 gC/m2)
#' @param grass_harvest_file file containing grazing data to rescale the
#'        grassland harvests according to Herrero et al. 2013. File contains:
#'        grazing_data list object with $name and $id of 29 world regions, and
#'        $Herrero_2000_kgDM_by_region containing for each of these regions and
#'        mapping_lpj67420_to_grazing_regions array with a mapping between 67420
#'        LPJmL cells and the 29 regions
#' @param external_fire_file path to external file with human induced fire
#'        fraction c(cell,month,year) since 1500
#' @param external_wood_harvest_file path to R-file containing processed
#'        timeline of maps for LUH2_v2h woodharvest
#'
#' @return pb_status list data object
#'
#' @examples
#' \dontrun{
#' }
#' @export
calc_biosphere_status <- function( files_scenario,
                            files_reference,
                            time_span_scenario,
                            time_span_reference = NULL,
                            reference_npp_time_span = NULL,
                            reference_npp_file = NULL,
                            gridbased = T,
                            readPreviouslySavedData = FALSE,
                            saveDataFile = FALSE,
                            dataFile = NULL,
                            include_fire = F,
                            external_fire = F,
                            external_wood_harvest = F,
                            grass_scaling = F,
                            npp_threshold = 20,
                            grass_harvest_file = "/p/projects/open/Fabian/LPJbox/grazing_data.RData",
                            external_fire_file = "/p/projects/open/Fabian/LPJbox/human_ignition_fraction.RData",
                            external_wood_harvest_file = "/p/projects/open/LanduseData/LUH2_v2h/wood_harvest_biomass_sum_1500-2014_67420.RData"
) {
  if (is.null(time_span_reference)) time_span_reference <- time_span_scenario
  if (grass_scaling && !file.exists(grass_harvest_file)) stop(paste0("Grass harvest scaling enabled, but grass_harvest_file does not exist in: ",grass_harvest_file))
  if (external_wood_harvest && !file.exists(external_wood_harvest_file)) stop(paste0("External wood harvest enabled, but external_wood_harvest_file does not exist in: ",external_wood_harvest_file))
  if (external_fire && !file.exists(external_fire_file)) stop(paste0("External fire fraction file enabled, but external_fire_file does not exist in: ",external_fire_file))
  # reading required data
  if (readPreviouslySavedData) {
    if (file.exists(dataFile)) {
      print(paste0("Reading in data from previously saved data file"))
      load(dataFile)
      wood_harvest[is.na(wood_harvest)] <- 0
    }else{
      stop(paste0("dataFile: '",dataFile,"' does not exist but is required since reading is set to FALSE."))
    }
    if (saveDataFile) {
      saveDataFile <- FALSE
      print("Both readPreviouslySavedData and saveDataFile have been set to TRUE. Overwriting with the same data does not make sense, saving disabled. ")
    }
  }else{
    print(paste0("Reading in data from outputs"))

    fileType <- tools::file_ext(files_reference$grid)

    if (fileType %in% c("json","clm")) {
      # read grid
      grid <- lpjmlkit::read_io(
        files_reference$grid,
      )$data %>% drop()
      # calculate cell area
      cell_area <- lpjmlkit::calc_cellarea(grid[, 2])
      ncells <- length(grid[, 1])

      npp <- lpjmlkit::read_io(
        files_scenario$npp,
        subset = list(year = as.character(time_span_scenario))) %>%
        lpjmlkit::transform(to = c("year_month_day")) %>%
        lpjmlkit::as_array(aggregate = list(month = sum)) %>% drop() # gC/m2

      if (!is.null(reference_npp_file)){
        npp_ref <- lpjmlkit::read_io(
          reference_npp_file,
          subset = list(year = as.character(reference_npp_time_span))) %>%
          lpjmlkit::transform(to = c("year_month_day")) %>%
          lpjmlkit::as_array(aggregate = list(month = sum)) %>% drop()#rem bands
      }

      pftnpp <- lpjmlkit::read_io(
        files_scenario$pft_npp,
        subset = list(year = as.character(time_span_scenario))) %>%
        lpjmlkit::transform(to = c("year_month_day")) %>%
        lpjmlkit::as_array(aggregate = list(month = sum))

      harvest <- lpjmlkit::read_io(
        files_scenario$pft_harvestc,
        subset = list(year = as.character(time_span_scenario))) %>%
        lpjmlkit::transform(to = c("year_month_day")) %>%
        lpjmlkit::as_array(aggregate = list(month = sum))

      rharvest <- lpjmlkit::read_io(
        files_scenario$pft_rharvestc,
        subset = list(year = as.character(time_span_scenario))) %>%
        lpjmlkit::transform(to = c("year_month_day")) %>%
        lpjmlkit::as_array(aggregate = list(month = sum))

      timber <- lpjmlkit::read_io(
        files_scenario$timber_harvestc,
        subset = list(year = as.character(time_span_scenario))) %>%
        lpjmlkit::transform(to = c("year_month_day")) %>%
        lpjmlkit::as_array(aggregate = list(month = sum)) %>% drop()#rem bands
      if (include_fire){
        # read fire in monthly res. if possible, then multiply with monthly
        # human/total ignition frac and aggregate to yearly. Otherwise aggregate
        # human/total ignition frac to yearly and multiply with yearly firec
        fire_raw <- lpjmlkit::read_io(
          files_scenario$firec,
          subset = list(year = as.character(time_span_scenario))) %>%
          lpjmlkit::transform(to = c("year_month_day")) %>%
          lpjmlkit::as_array(aggregate = list(band = sum)) #gC/m2
        if (external_fire) {
          load(external_fire_file) #frac = c(cell,month,year)
        }
        if ("month" %in% names(dim(fire_raw))) {
          if (external_fire) {
            fire <- apply(fire_raw*frac[,,year = time_span_scenario], c("cell","year"), sum, na.rm = T) #gC/m2
            rm(frac)
          }else{
            fire <- apply(fire_raw, c("cell","year"), sum, na.rm = T) #gC/m2
          }
          rm(fire_raw)
        }else{
          if (external_fire) {
            frac_yearly <- apply(frac[,,year = time_span_scenario],c("cell","year"),mean,na.rm = T)
            fire <- fire_raw*frac_yearly
            rm(frac_yearly,frac)
          }
        }
        gc()
      }else{
        fire <- timber*0
      }
      if (external_wood_harvest) {
        load(external_wood_harvest_file) #wh_lpj in kgC
        wh_years <- names(wh_lpj[1,])
        wood_harvest <- wh_lpj[,match(time_span_scenario,wh_years)]*10^3/cell_area #from kgC to gC/m2
        wood_harvest[is.na(wood_harvest)] <- 0 # the division can lead to NAs
        rm(wh_lpj,wh_years)
        gc()
      }else{
        wood_harvest <- fire*0
      }

      cftfrac <- lpjmlkit::read_io(
        files_scenario$cftfrac,
        subset = list(year = as.character(time_span_scenario))) %>%
        lpjmlkit::transform(to = c("year_month_day")) %>%
        lpjmlkit::as_array(aggregate = list(month = sum))

      npp_potential <- lpjmlkit::read_io(
        files_reference$npp,
        subset = list(year = as.character(time_span_reference))) %>%
        lpjmlkit::transform(to = c("year_month_day")) %>%
        lpjmlkit::as_array(aggregate = list(month = sum)) %>% drop() # gC/m2

      fpc <- lpjmlkit::read_io(
        files_scenario$fpc,
        subset = list(year = as.character(time_span_scenario))) %>%
        lpjmlkit::transform(to = c("year_month_day")) %>%
        lpjmlkit::as_array(aggregate = list(band = sum))

      cftbands <- lpjmlkit::read_meta(files_scenario$cftfrac)$nbands
      pftbands <- lpjmlkit::read_meta(files_scenario$fpc)$nbands - 1

    }else if (fileType == "nc") { # to be added
      stop("nc reading has not been updated to latest functionality. please contact Fabian")
    }else{
      stop(paste0("Unrecognized file type (",fileType,")"))
    }
    bp_bands <- c(15,16,31,32)
    grass_bands <- c(14,30)
    nat_bands <- 1:pftbands
    if (!gridbased) { # needs to be scaled with standfrac
      pftnpp[,,nat_bands] <- pftnpp[,,nat_bands]*fpc[,,2:(pftbands+1)]
      pftnpp[,,-c(nat_bands)] <- pftnpp[,,-c(nat_bands)]*cftfrac
      harvest <- harvest*cftfrac
    }
    pftnpp_grasslands <- apply(pftnpp[,,pftbands+grass_bands ],c(1,2),sum) #gC/m2 only from grassland bands
    pftnpp_cft <- apply(pftnpp[,,-c(nat_bands,pftbands+grass_bands,pftbands+bp_bands)], c(1,2), sum) #gC/m2 not from grassland and bioenergy bands
    pftnpp_bioenergy <- apply(pftnpp[,, pftbands+bp_bands ], c(1,2), sum) #gC/m2 only from bioenergy bands
    pftnpp_nat <- apply(pftnpp[,,nat_bands], c(1,2), sum) #gC/m2

    if (is.null(reference_npp_file)){
      pi_window <- 3:32
      npp_ref <- npp_potential[,pi_window]
    } # npp_ref

    harvest_grasslands <- apply(harvest[,,grass_bands],c(1,2),sum) #gC/m2 only from grassland bands
    harvest_bioenergy <- apply(harvest[,,bp_bands],c(1,2),sum) #gC/m2 only from bioenergy bands
    harvest_cft <- apply(harvest[,,-c(grass_bands,bp_bands)], c(1,2), sum) #gC/m2 not from grassland and bioenergy bands
    rharvest_cft <- apply(rharvest[,,-c(grass_bands,bp_bands)], c(1,2), sum) #gC/m2 not from grassland and bioenergy bands

    if (saveDataFile) {
      if (!file.exists(dataFile) ) {
        print(paste0("Writing data file: ",dataFile))
      }else{
        print(paste0("Data file (",dataFile,") already exists, old file renamed to: ",dataFile,"_sav"))
        file.rename(dataFile, paste0(dataFile,"_sav"))
      }
      save(npp_potential,npp,npp_ref,pftnpp_cft,pftnpp_nat,pftnpp_grasslands,
           pftnpp_bioenergy,harvest_cft,rharvest_cft,fire,timber,fpc,cftfrac,
           harvest_grasslands,harvest_bioenergy,wood_harvest,file = dataFile)
    }

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
      grassland_scaling_factor_cellwise[
          which(mapping_lpj67420_to_grazing_regions == r)] <- scaling_factor[r]
    }
    #plotGlobalMan(data = grassland_scaling_factor_cellwise,file = "~/scaling_factor_lpj_grazing.png",
    #         title = "grassland scaling factor",brks = c(seq(0,1,length.out = 11),seq(2,10,length.out = 9),60),
    #         palette = "RdBu",legendtitle = "",legYes = T,eps=F)
    harvest_grasslands <- harvest_grasslands*rep(grassland_scaling_factor_cellwise,
                                                 times = length(harvest_grasslands[1,]))
  }

  npp_act_overtime <- colSums(npp*cellarea)/10^15 # from gC/m2 to GtC
  npp_pot_overtime <- colSums(npp_potential*cellarea)/10^15 # from gC/m2 to GtC
  npp_eco_overtime <- colSums(pftnpp_nat*cellarea)/10^15 # from gC/m2 to GtC
  npp_luc_overtime <- npp_pot_overtime - npp_act_overtime

  harvest_cft_overtime <- colSums(harvest_cft*cellarea)/10^15 # from gC/m2 to GtC
  rharvest_cft_overtime <- colSums(rharvest_cft*cellarea)/10^15 # from gC/m2 to GtC
  harvest_grasslands_overtime <- colSums(harvest_grasslands*cellarea)/10^15 # from gC/m2 to GtC
  harvest_bioenergy_overtime <- colSums(harvest_bioenergy*cellarea)/10^15 # from gC/m2 to GtC

  timber_harvest_overtime <- colSums(timber*cellarea)/10^15 # from gC/m2 to GtC
  fire_overtime <- colSums(fire*cellarea)/10^15 # from gC/m2 to GtC
  wood_harvest_overtime <- colSums(wood_harvest*cellarea)/10^15 # from gC/m2 to GtC

  if (include_fire) {
    mcol_overtime <- harvest_cft_overtime + rharvest_cft_overtime +
      harvest_grasslands_overtime + harvest_bioenergy_overtime +
      timber_harvest_overtime + fire_overtime + npp_luc_overtime +
      wood_harvest_overtime
  }else{
    mcol_overtime <- harvest_cft_overtime + rharvest_cft_overtime +
      harvest_grasslands_overtime + harvest_bioenergy_overtime +
      timber_harvest_overtime + npp_luc_overtime +
      wood_harvest_overtime
  }

  mcol_overtime_perc_piref <- mcol_overtime/mean(colSums(npp_ref*cellarea)/10^15)*100
  mcol_luc <- npp_potential - npp
  # pick a PI window that excludes onset effects, but is reasonable early
  #mcol_luc_piref <- rep(rowMeans(npp_potential[,pi_window]),times = length(npp[1,])) - npp # not used, but handed over for checking comparison to pi_ref

  if (include_fire) {
    mcol_harvest <- harvest_cft + rharvest_cft + harvest_grasslands + harvest_bioenergy + timber + fire + wood_harvest
  }else{
    mcol_harvest <- harvest_cft + rharvest_cft + harvest_grasslands + harvest_bioenergy + timber + wood_harvest
  }
  mcol <- mcol_harvest + mcol_luc
  mcol[abs(npp_potential) < npp_threshold] <- 0 # set to 0 below lower threshold of NPP
  mcol_perc <- mcol/npp_potential*100 #actual NPPpot as ref


  mcol_perc_piref <- mcol/rowMeans(npp_ref)*100 # NPPpi as ref

  # todo: add param spatial_resolution "global"/"grid" switch
  # add avrg_args treatment
  pb_status <- list()
  pb_status$control_variable <- mcol_perc_piref
  # pb value normalized to 0-1, 1 is the
  # threshold between safe and increasing risk, >1 is transgressed
  pb_status$continuous_normalized <-
  pb_status$traffic_light <- ifelse(pb_status$control_variable > q95_area,
                                  2, 1) # only green and transgresssed
  pb_status$thresholds <- c(holocene = q50_area, pb = q95_area) # h: 0, pb: 1

  return(pb_status)

} # end of calc_biosphere_status

#' Calculate MCOL
#'
#' Wrapper function to calculate MCOL
#'
#' @param inFol_lu folder of landuse scenario run
#' @param inFol_pnv folder of pnv reference run
#' @param startyr first year of simulations
#' @param stopyr last year of simulations
#' @param reference_npp_time_span time span to read reference npp from, using
#'     index years 10:39 from potential npp input if set to NULL (default: NULL)
#' @param reference_npp_file file to read reference npp from, using
#'        potential npp input if set to NULL (default: NULL)
#' @param gridbased logical are pft outputs gridbased or pft-based?
#' @param readPreviouslySavedData flag whether to read previously saved data
#'        instead of reading it in from output files (default FALSE)
#' @param saveDataFile whether to save input data to file (default FALSE)
#' @param dataFile file to save/read input data to/from (default NULL)
#' @param include_fire boolean include firec in calculation of MCOL? (default T)
#' @param external_fire instead of reading in firec for fire emissions, read in
#'        this external firec file from a separate spitfire run with disabled
#'        lighning. this will then include only human induced fires (default F)
#' @param external_wood_harvest include external wood harvest from LUH2_v2h
#'        (default F)
#' @param grass_scaling whether to scale pasture harvest according to
#'        data given via grass_harvest_file (default F)
#' @param npp_threshold lower threshold for npp (to mask out non-lu areas
#'        according to Haberl et al. 2007). Below MCOL will be set to 0.
#'        (default: 20 gC/m2)
#' @param grass_harvest_file file containing grazing data to rescale the
#'        grassland harvests according to Herrero et al. 2013. File contains:
#'        grazing_data list object with $name and $id of 29 world regions, and
#'        $Herrero_2000_kgDM_by_region containing for each of these regions and
#'        mapping_lpj67420_to_grazing_regions array with a mapping between 67420
#'        LPJmL cells and the 29 regions
#' @param external_fire_file path to external file with human induced fire
#'        fraction c(cell,month,year) since 1500
#' @param external_wood_harvest_file path to R-file containing processed
#'        timeline of maps for LUH2_v2h woodharvest
#'
#' @return list data object containing MCOL and components as arrays: mcol,
#'         mcol_overtime, mcol_overtime_perc_piref, mcol_perc, mcol_perc_piref,
#'         npp_potential, npp_act_overtime, npp_pot_overtime, npp_eco_overtime,
#'         harvest_cft_overtime, npp_luc_overtime, rharvest_cft_overtime,
#'         fire_overtime, timber_harvest_overtime, harvest_cft, mcol_harvest,
#'         grassland_scaling_factor_cellwise,  mcol_luc, mcol_luc_piref
#'
#' @examples
#' \dontrun{
#' }
#' @export
calc_mcol <- function( inFol_lu,
                      inFol_pnv,
                      startyr,
                      stopyr,
                      reference_npp_time_span = NULL,
                      reference_npp_file = NULL,
                      varnames = NULL,
                      gridbased = T,
                      readPreviouslySavedData = FALSE,
                      saveDataFile = FALSE,
                      dataFile = NULL,
                      include_fire = F,
                      external_fire = F,
                      external_wood_harvest = F,
                      grass_scaling = F,
                      npp_threshold = 20,
                      grass_harvest_file = "/p/projects/open/Fabian/LPJbox/grazing_data.RData",
                      external_fire_file = "/p/projects/open/Fabian/LPJbox/human_ignition_fraction.RData",
                      external_wood_harvest_file = "/p/projects/open/LanduseData/LUH2_v2h/wood_harvest_biomass_sum_1500-2014_67420.RData"
) {
  if (is.null(varnames)) {
    print("Varnames not given, using standard values, which might not fit this specific configuration. Please check!")
    varnames <- data.frame(row.names = c("grid","npp","pft_npp","pft_harvest","pft_rharvest","firec","timber_harvest","cftfrac","fpc"),
                           outname = c("grid.bin.json",  "mnpp.bin.json","pft_npp.bin.json","pft_harvest.bin.json","pft_rharvest.bin.json","firec.bin.json","timber_harvestc.bin.json","cftfrac.bin.json","fpc.bin.json"),
                           timestep = c("Y","Y","Y","Y","Y","Y","Y","Y","Y"))
  }

  # translate varnames and folders to files_scenarios/reference lists
  files_scenario <- list(grid = paste0(inFol_lu,varnames["grid","outname"]),
                         npp = paste0(inFol_lu,varnames["npp","outname"]),
                         pft_npp = paste0(inFol_lu,varnames["pft_npp","outname"]),
                         pft_harvestc = paste0(inFol_lu,varnames["pft_harvest","outname"]),
                         pft_rharvestc = paste0(inFol_lu,varnames["pft_rharvest","outname"]),
                         firec = paste0(inFol_lu,varnames["firec","outname"]),
                         timber_harvestc = paste0(inFol_lu,varnames["timber_harvest","outname"]),
                         cftfrac = paste0(inFol_lu,varnames["cftfrac","outname"]),
                         fpc = paste0(inFol_lu,varnames["fpc","outname"])
  )
  files_reference <- list(grid = paste0(inFol_pnv,varnames["grid","outname"]),
                          npp = paste0(inFol_pnv,varnames["npp","outname"]),
                          pft_npp = paste0(inFol_pnv,varnames["pft_npp","outname"]),
                          pft_harvestc = paste0(inFol_pnv,varnames["pft_harvest","outname"]),
                          pft_rharvestc = paste0(inFol_pnv,varnames["pft_rharvest","outname"]),
                          firec = paste0(inFol_pnv,varnames["firec","outname"]),
                          timber_harvestc = paste0(inFol_pnv,varnames["timber_harvest","outname"]),
                          cftfrac = paste0(inFol_pnv,varnames["cftfrac","outname"]),
                          fpc = paste0(inFol_pnv,varnames["fpc","outname"])
  )
  # todo call calc_biosphere_status
  return(calc_biosphere_status(files_scenario = files_scenario,
                        files_reference = files_reference,
                        time_span_scenario = as.character(startyr:stopyr),
                        time_span_reference = as.character(startyr:stopyr),
                        reference_npp_time_span = reference_npp_time_span,
                        reference_npp_file = reference_npp_file,
                        gridbased = gridbased,
                        readPreviouslySavedData = readPreviouslySavedData,
                        saveDataFile = saveDataFile,
                        dataFile = dataFile,
                        include_fire = include_fire,
                        external_fire = external_fire,
                        external_wood_harvest = external_wood_harvest,
                        grass_scaling = grass_scaling,
                        npp_threshold = npp_threshold,
                        grass_harvest_file = grass_harvest_file,
                        external_fire_file = external_fire_file,
                        external_wood_harvest_file = external_wood_harvest_file
  ) )

}

#' Plot absolute MCOL, overtime, maps, and npp into given folder
#'
#' Wrapper function to plot absolute mcol, overtime, maps, and npp into given folder
#'
#' @param mcolData mcol data list object (returned from calc_mcol) containing
#'                 mcol, npp_eco_overtime, npp_act_overtime, npp_pot_overtime,
#'                 npp_bioenergy_overtime, mcol_overtime, npp_harv_overtime,
#'                 mcol_overtime_perc_piref, mcol_perc, mcol_perc_piref, npp
#'                 all in GtC
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
                     startyr, mapyear, mapyear_buffer = 5, highlightyear, eps = FALSE) {
  mapindex <- mapyear - startyr
  print(paste0("Plotting MCOL figures"))
  dir.create(file.path(outFol),showWarnings = F)
  lpjmliotools::plotGlobal(data = rowMeans(mcolData$mcol[,(mapindex-mapyear_buffer):(mapindex+mapyear_buffer)]),
                           file = paste0(outFol,"MCOL_absolute_",mapyear,".png"), type = "exp",
                           title = paste0("MCOL_abs in ",mapyear), pow2min = 0, pow2max = 12,
                           legendtitle = "GtC", legYes = T, onlyPos = F, eps = eps)
  lpjmliotools::plotGlobal(data = rowMeans(mcolData$mcol_luc[,(mapindex-mapyear_buffer):(mapindex+mapyear_buffer)]),
                           file = paste0(outFol,"MCOL_luc_",mapyear,".png"), type = "exp",
                           title = paste0("MCOL_luc in ",mapyear), pow2min = 0, pow2max = 12,
                           legendtitle = "GtC", legYes = T, onlyPos = F, eps = eps)
  #lpjmliotools::plotGlobal(data = rowMeans(mcolData$mcol_luc_piref[,(mapindex-mapyear_buffer):(mapindex+mapyear_buffer)]),
  #                         file = paste0(outFol,"MCOL_luc_piref_",mapyear,".png"), type = "exp",
  #                         title = paste0("MCOL_luc piref in ",mapyear), pow2min = 0, pow2max = 12,
  #                         legendtitle = "GtC", legYes = T, onlyPos = F, eps = eps)
  lpjmliotools::plotGlobal(data = rowMeans(mcolData$mcol_harvest[,(mapindex-mapyear_buffer):(mapindex+mapyear_buffer)]),
                           file = paste0(outFol,"MCOL_harv_",mapyear,".png"), type = "exp",
                           title = paste0("MCOL_harv in ",mapyear), pow2min = 0, pow2max = 12,
                           legendtitle = "GtC", legYes = T, onlyPos = F, eps = eps)
  plotMCOLovertime(mcolData = mcolData, file = paste0(outFol,"MCOL_overtime_LPJmL_",plotyears[1],"-",plotyears[2],".png"),
                   firstyr = startyr, plotyrs = plotyears, minVal = minVal, ref = "pi",
                   legendpos = legendpos, maxVal = maxVal, eps = eps, highlightyrs = highlightyear)
  plotMCOLmap(data = rowMeans(mcolData$mcol_perc[,(mapindex-mapyear_buffer):(mapindex+mapyear_buffer)]),
              file = paste0(outFol,"MCOL_LPJmL_",mapyear,".png"),legendtitle = "% of NPPpot", eps = eps,
              title = paste0("MCOL_perc ",mapyear-mapyear_buffer, " - ",mapyear+mapyear_buffer) )
  plotMCOLmap(data = rowMeans(mcolData$mcol_perc_piref[,(mapindex-mapyear_buffer):(mapindex+mapyear_buffer)]),
              file = paste0(outFol,"MCOL_piref_LPJmL_",mapyear,".png"),
              title = paste0("MCOL_perc ",mapyear-mapyear_buffer, " - ",mapyear+mapyear_buffer),legendtitle = "% of NPPpi", eps = eps)
  lpjmliotools::plotGlobalMan(data = rowMeans(mcolData$npp[,(mapindex-mapyear_buffer):(mapindex+mapyear_buffer)]),
                              file = paste0(outFol,"NPP_LPJmL_",mapyear,".png"), brks = seq(0,1000,100),
                              palette = c("orangered4","orangered","orange","gold","greenyellow","limegreen","green4","darkcyan","darkslategrey","navy"),
                              title = paste0("NPP average ",mapyear-mapyear_buffer, "-",mapyear+mapyear_buffer),
                              legendtitle = "gC/m2",legYes = T)

} # end of plotMCOL

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
plotMCOLmap <- function(data, file, title = "", legendtitle = "", zeroThreshold = 0.1, eps = FALSE) {
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
#' @param mcolData mcol data list object (returned from calc_mcol) containing
#'                  mcol, npp_eco_overtime, npp_act_overtime, npp_pot_overtime,
#'                  npp_bioenergy_overtime, mcol_overtime, npp_harv_overtime,
#'                  mcol_overtime_perc_piref, mcol_perc, mcol_perc_piref, npp
#'                 all in GtC
#' @param file character string for location/file to save plot to
#' @param firstyr first year of mcol object
#' @param plotyrs range of years to plot over time
#' @param highlightyrs year(s) that should be highlighted in overtime plot (default: 2000)
#' @param minVal y-axis minimum value for plot over time (default: 0)
#' @param maxVal y-axis maximum value for plot over time (default: 100)
#' @param legendpos position of legend (default: "topleft")
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
                             maxVal = 100, legendpos = "topleft", ext = FALSE, eps = FALSE, ref = "pi") {
  lastyr = firstyr + length(mcolData$npp_act_overtime) - 1
  colz = c("slateblue","gold","green3","darkorange","black","red3","green","brown","yellow","turquoise","darkgreen")
  if (eps) {
    file = strsplit(file,".",fixed=TRUE)[[1]]
    file = paste(c(file[1:(length(file) - 1)],"eps"),collapse=".")
    ps.options(family = c("Helvetica"), pointsize = 18)
    postscript(file,horizontal = FALSE, onefile = FALSE, width = 22,  height = 8.5,paper="special")
  }else{
    png(file, width=3.5,  height = 3, units = "in", res = 300, pointsize = 6,type="cairo")
  }
  par(bty="o",oma=c(0,0,0,0),mar=c(4,5,1,3))
  plot(NA,ylab="GtC/yr",xlab="Year",xlim=plotyrs,
       ylim=c(minVal, maxVal),xaxs="i",yaxs="i")
  grid()
  lines(x=seq(firstyr,lastyr,1),y=mcolData$npp_pot_overtime,type = "l",col=colz[1])
  lines(x=seq(firstyr,lastyr,1),y=mcolData$npp_act_overtime,type = "l",col=colz[2])
  lines(x=seq(firstyr,lastyr,1),y=mcolData$npp_eco_overtime,type = "l",col=colz[3])
  lines(x=seq(firstyr,lastyr,1),y=mcolData$npp_luc_overtime,type = "l",col=colz[4])
  lines(x=seq(firstyr,lastyr,1),y=mcolData$mcol_overtime,type = "l",col=colz[5])
  lines(x=seq(firstyr,lastyr,1),y=mcolData$harvest_cft_overtime,type = "l",col=colz[7])
  lines(x=seq(firstyr,lastyr,1),y=mcolData$rharvest_cft_overtime,type = "l",col=colz[8])
  lines(x=seq(firstyr,lastyr,1),y=mcolData$fire_overtime,type = "l",col=colz[9])
  lines(x=seq(firstyr,lastyr,1),y=mcolData$timber_harvest_overtime,type = "l",col=colz[10])
  lines(x=seq(firstyr,lastyr,1),y=mcolData$wood_harvest_overtime,type = "l",col=colz[11])


  par(bty="n",oma=c(0,0,0,0),mar=c(4,5,1,3), new = T)
  if (ref == "pi") {
    plot(x=seq(firstyr,lastyr,1),y=mcolData$mcol_overtime_perc_piref,ylab="",xlab="",xlim=plotyrs,
         ylim=c(0, 30),type = "l",col=colz[6],xaxs="i", yaxs="i", axes = F)
  } else if (ref == "act") {
    plot(x=seq(firstyr,lastyr,1),y=mcolData$mcol_overtime,ylab="",xlab="",xlim=plotyrs,
         ylim=c(0, 30),type = "l",col=colz[6],xaxs="i", yaxs="i", axes = F)
  }else stop(paste0("Unknown value for parameter ref: ",ref," - Aborting."))

  axis(side = 4, col = colz[6],col.axis = colz[6])
  mtext(text = "%", col=colz[6], side = 4,line = 2)

  if (!is.null(highlightyrs)) {
    for (y in highlightyrs) {
      lines(x=c(y,y),y=c(minVal,maxVal),col="grey40")
    }
  }
  legend(legendpos,legend = c("NPPpot (PNV)","NPPact (landuse)","NPPeco","NPPluc","M-COLabs","M-COL [% NPPpi]","harvestc","rharvest","firec","timber_harvest","wood_harvest"),col=colz ,lty=1,cex=1)
  dev.off()
}
