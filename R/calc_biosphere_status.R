#' Calculate biosphere status based on BioCol (HANPP) from a PNV run (reference) and
#' LU run (scenario) of LPJmL, both using time_span_scenario. Additionally
#' a separate reference NPP file (e.g. from a Holocene run) can be supplied as
#' reference_npp_file, which will use time_span_reference, or file index years
#' 3:32 if time_span_reference is not supplied
#'
#' Calculate biosphere status based on BioCol from a PNV run and LU run of LPJmL
#' @param files_scenario list with variable names and corresponding file paths
#' (character string) of the scenario LPJmL run. All needed files need to be
#' provided. E.g.: list(leaching = "/temp/leaching.bin.json",
#'                      npp = "/temp/leaching.bin.json")
#' @param files_reference list with variable names and corresponding file paths
#' (character string) of the PNV LPJmL run. All needed files need to be
#' provided. In this case only NPP. list(npp = "/temp/leaching.bin.json").
#' @param time_span_scenario time span to be used for the scenario run, defined
#' as a character string, e.g. `as.character(1982:2011)` (default)
#' @param time_span_reference time span to read reference_npp_file from, using
#'     index years 3:32 from potential npp input if set to NULL (default: NULL)
#' @param spatial_resolution character string indicating spatial resolution
#'        either "grid" or "global"
#' @param thresholds named character string with thresholds to be used to
#'        define the lower end of safe, increasing risk and high risk zone,
#'        e.g. c(holocene = 0.0, pb = 0.1, highrisk = 0.2).
#' @param reference_npp_file file to read reference npp from, using
#'        potential npp reference input if set to NULL (default: NULL)
#' @param gridbased logical are pft outputs gridbased or pft-based?
#' @param include_fire boolean include firec in calculation of BioCol? (default T)
#' @param external_fire instead of reading in firec for fire emissions, read in
#'        this external firec file from a separate spitfire run with disabled
#'        lighning. this will then include only human induced fires (default F)
#' @param external_wood_harvest include external wood harvest from LUH2_v2h
#'        (default F)
#' @param grass_scaling whether to scale pasture harvest according to
#'        data given via grass_harvest_file (default F)
#' @param npp_threshold lower threshold for npp (to mask out non-lu areas
#'        according to Haberl et al. 2007). Below BioCol will be set to 0.
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
calc_biosphere_status <- function(files_scenario,
                                  files_reference,
                                  time_span_scenario,
                                  time_span_reference = NULL,
                                  spatial_resolution,
                                  thresholds = NULL,
                                  reference_npp_file = NULL,
                                  gridbased = T,
                                  include_fire = F,
                                  external_fire = F,
                                  external_wood_harvest = F,
                                  grass_scaling = F,
                                  npp_threshold = 20,
                                  grass_harvest_file = "/p/projects/open/Fabian/LPJbox/grazing_data.RData",
                                  external_fire_file = "/p/projects/open/Fabian/LPJbox/human_ignition_fraction.RData",
                                  external_wood_harvest_file = "/p/projects/open/LanduseData/LUH2_v2h/wood_harvest_biomass_sum_1500-2014_67420.RData"
) {
  if (grass_scaling && !file.exists(grass_harvest_file)) stop(paste0("Grass harvest scaling enabled, but grass_harvest_file does not exist in: ",grass_harvest_file))
  if (external_wood_harvest && !file.exists(external_wood_harvest_file)) stop(paste0("External wood harvest enabled, but external_wood_harvest_file does not exist in: ",external_wood_harvest_file))
  if (external_fire && !file.exists(external_fire_file)) stop(paste0("External fire fraction file enabled, but external_fire_file does not exist in: ",external_fire_file))
  # reading required data
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

    # PNV NPP overtime (time_span_scenario!), required to calculate NPP_luc
    npp_potential <- lpjmlkit::read_io(
      files_reference$npp,
      subset = list(year = as.character(time_span_scenario))) %>%
      lpjmlkit::transform(to = c("year_month_day")) %>%
      lpjmlkit::as_array(aggregate = list(month = sum)) %>% drop() # gC/m2

    # reference NPP, only used to translate absolute HANPP to relative BioCol
    if (is.null(reference_npp_file)) { # use PNV NPP
      if (is.null(time_span_reference)) { # use first 30 years from PNV NPP
        npp_ref <- lpjmlkit::read_io(
          files_reference$npp,
          subset = list(year = 3:32)) %>%
          lpjmlkit::transform(to = c("year_month_day")) %>%
          lpjmlkit::as_array(aggregate = list(month = sum)) %>% drop()#rem bands
      } else{
        npp_ref <- lpjmlkit::read_io(
          files_reference$npp,
          subset = list(year = as.character(time_span_reference))) %>%
          lpjmlkit::transform(to = c("year_month_day")) %>%
          lpjmlkit::as_array(aggregate = list(month = sum)) %>% drop()#rem bands
      }
    } else { # separate reference NPP file
      if (is.null(time_span_reference)) { # use first 30 years from ref NPP
        npp_ref <- lpjmlkit::read_io(
          reference_npp_file,
          subset = list(year = 3:32)) %>%
          lpjmlkit::transform(to = c("year_month_day")) %>%
          lpjmlkit::as_array(aggregate = list(month = sum)) %>% drop()#rem bands
      } else{
        npp_ref <- lpjmlkit::read_io(
          reference_npp_file,
          subset = list(year = as.character(time_span_reference))) %>%
          lpjmlkit::transform(to = c("year_month_day")) %>%
          lpjmlkit::as_array(aggregate = list(month = sum)) %>% drop()#rem bands
      }
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

  harvest_grasslands <- apply(harvest[,,grass_bands],c(1,2),sum) #gC/m2 only from grassland bands
  harvest_bioenergy <- apply(harvest[,,bp_bands],c(1,2),sum) #gC/m2 only from bioenergy bands
  harvest_cft <- apply(harvest[,,-c(grass_bands,bp_bands)], c(1,2), sum) #gC/m2 not from grassland and bioenergy bands
  rharvest_cft <- apply(rharvest[,,-c(grass_bands,bp_bands)], c(1,2), sum) #gC/m2 not from grassland and bioenergy bands


  print(paste0("Calculating data"))
  if (grass_scaling) {
    load(grass_harvest_file)
    nregs <- length(grazing_data$name)
    lpj_grass_harvest_region <- array(0,dim = nregs)
    lpj_grass_harvest_2000 <- rowMeans(harvest_grasslands[,(1995-startyr+1):(2005-startyr+1)])*cell_area/1000*2 # from gC/m2 to kgDM
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

  npp_act_overtime <- colSums(npp*cell_area)/10^15 # from gC/m2 to GtC
  npp_pot_overtime <- colSums(npp_potential*cell_area)/10^15 # from gC/m2 to GtC
  npp_luc_overtime <- npp_pot_overtime - npp_act_overtime

  harvest_cft_overtime <- colSums(harvest_cft*cell_area)/10^15 # from gC/m2 to GtC
  rharvest_cft_overtime <- colSums(rharvest_cft*cell_area)/10^15 # from gC/m2 to GtC
  harvest_grasslands_overtime <- colSums(harvest_grasslands*cell_area)/10^15 # from gC/m2 to GtC
  harvest_bioenergy_overtime <- colSums(harvest_bioenergy*cell_area)/10^15 # from gC/m2 to GtC

  timber_harvest_overtime <- colSums(timber*cell_area)/10^15 # from gC/m2 to GtC
  fire_overtime <- colSums(fire*cell_area)/10^15 # from gC/m2 to GtC
  wood_harvest_overtime <- colSums(wood_harvest*cell_area)/10^15 # from gC/m2 to GtC

  if (include_fire) {
    biocol_overtime <- harvest_cft_overtime + rharvest_cft_overtime +
      harvest_grasslands_overtime + harvest_bioenergy_overtime +
      timber_harvest_overtime + fire_overtime + npp_luc_overtime +
      wood_harvest_overtime
  }else{
    biocol_overtime <- harvest_cft_overtime + rharvest_cft_overtime +
      harvest_grasslands_overtime + harvest_bioenergy_overtime +
      timber_harvest_overtime + npp_luc_overtime +
      wood_harvest_overtime
  }

  biocol_overtime_piref <- biocol_overtime/mean(colSums(npp_ref*cell_area)/10^15)
  biocol_luc <- npp_potential - npp

  if (include_fire) {
    biocol_harvest <- harvest_cft + rharvest_cft + harvest_grasslands + harvest_bioenergy + timber + fire + wood_harvest
  }else{
    biocol_harvest <- harvest_cft + rharvest_cft + harvest_grasslands + harvest_bioenergy + timber + wood_harvest
  }
  biocol <- biocol_harvest + biocol_luc
  biocol[abs(npp_potential) < npp_threshold] <- 0 # set to 0 below lower threshold of NPP

  ref_npp <- rowMeans(npp_ref)
  biocol_piref <- biocol/ref_npp # NPPpi as ref
  biocol_piref[ref_npp == 0] <- 0

  if (is.null(thresholds)) {
    if (spatial_resolution == "grid") {
      thresholds <- c(holocene = 0,
                      pb = 0.05,
                      highrisk = 0.75)
    }else if (spatial_resolution == "global") {
      thresholds <- c(holocene = 0,
                      pb = 0.2,
                      highrisk = 0.5)
    }
  }
  if (spatial_resolution == "grid") {
    control_variable <- abs(biocol_piref)
    attr(control_variable, "thresholds") <- thresholds
  } else if (spatial_resolution == "global") {
    control_variable <- abs(biocol_overtime_piref)
    attr(control_variable, "thresholds") <- thresholds
  }else{
    stop(paste("Unknown value for spatial_resolution: ", spatial_resolution))
  }
  return(control_variable)

} # end of calc_biosphere_status
