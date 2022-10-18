# written by Fabian Stenzel, based on work by Sebastian Ostberg, Johanna Braun, Jannes Breier
# 2022 - stenzel@pik-potsdam.de

#requirements:
require(lpjmliotools) # in at least version 0.2.17

#' Classify biomes
#'
#' Classify biomes based on foliage protected cover (FPC) LPJmL output
#' and temperature output plus either vegetation carbon or pft_lai depending on
#' the savannaProxy option and elevation if montaneArcticProxy requires this
#'
#' @param data list object. Contains all relevant data for computation of the
#'        biome classes as yearly averages of monthly values.
#'        fpc [dim=c(ncells,npfts+1)],pft_lai [dim=c(ncells,npfts+ncfts)],
#'        lat, lon, vegc, temp [all dim=c(ncells)] e.g. list(lat = lat67420,
#'        lon = lon67420, fpc = apply(fpc_scen, c(1,2), mean),
#'        vegc = mean_state_scen[,8],temp = average_annual_temperature_1916)
#' @param readOutput read output from folder or use supplied data
#'        (requires folder and timespan to be set)
#' @param folder to read outputs from
#' @param files to read as list list(grid="grid.bin",fpc="fpc.bin")
#' @param timespan as c(startyear,stopyear) to use for averaging outputs over
#' @param savannaProxy "vegc" or "natLAI". Use vegetation carbon or LAI in natural
#'        vegetation as a proxy threshold to distinguish forests and savannahs
#' @param montaneArcticProxy "elevation" or "latitude". Use elevation or latitude
#'        as a proxy threshold to distinguish arctic tundra and montane grassland
#' @param elevation_threshold threshold in m above which ArcticTundra is
#'        classified as Montane Grassland, if montaneArcticProxy is set to
#'        elevation - default: 1000
#' @param latitude_threshold threshold in degrees, south of which ArcticTundra is
#'        classified as Montane Grassland, if montaneArcticProxy is set to
#'        latitude - default: 55
#' @param lai_threshold threshold for "natLAI" proxy (default: 6 m2/m2)
#' @param vegc_threshold threshold for "vegc" proxy (default: 7500 gC/m2)
#' @param tree_cover_thresholds list with minimum tree cover thresholds for
#'        definition of forest, woodland, savanna and grassland. Only changes to
#'        the default have to be included in the list, for the rest the default
#'        is used. Default values:
#'        "boreal forest" = 0.6
#'        "temperate forest" = 0.6
#'        "temperate woodland" = 0.3
#'        "temperate savannah" = 0.1
#'        "tropical forest" = 0.6
#'        "tropical woodland" = 0.3
#'        "tropical savannah" = 0.1
#' @param lpjGridInput path to lpjml grid input to be used for conversion
#'        only required for nc input
#'        (default: "/p/projects/lpjml/input/historical/input_VERSION2/grid.bin")
#' @param elevationInput path to lpjml elevation input to be used as proxy
#'        (default: "/p/projects/lpjml/input/historical/input_VERSION2/elevation.bin")
#' @param lpjGridHeaderSize header size for lpjml grid input (default: 43)
#'        only required for nc input
#' @param lpjCells number of grid cells in lpjml grid input (default: 67420)
#'        only required for nc input
#'
#' @return list object containing biome_id (main biome per grid cell [dim=c(ncells)]),
#' and list of respective biome_names[dim=c(nbiomes)]
#'
#' @examples
#' \dontrun{
#' classify_biomes(readOutput = T, timespan = c(1982,2011)
#'      folder = "/p/projects/open/Fabian/runs/Gamma/output/historic_gamma/",
#'      files = list(grid = "grid.bin", fpc = "fpc.bin", vegc = "vegc.bin",
#'      pft_lai = "pft_lai.bin", temp = "temp.bin"))
#' }
#'
#' @export
classify_biomes <- function(data = NULL, readOutput = F, folder = NULL, files = NULL,
                   timespan = NULL, savannaProxy = "natLAI", lai_threshold = 6,
                   vegc_threshold = 7500, montaneArcticProxy = "elevation",
                   elevation_threshold = 1000, latitude_threshold = 55, tree_cover_thresholds = list(),
                   lpjGridInput = "/p/projects/lpjml/input/historical/input_VERSION2/grid.bin",
                   elevationInput = "/p/projects/lpjml/input/historical/input_VERSION2/elevation.bin",
                   lpjGridHeaderSize = 43, lpjCells = 67420) {

  require(lpjmliotools)

  if (readOutput) { #reading output from folder
    if (!is.null(data)) stop("In readOutput mode, data cannot be supplied. Aborting.")
    if ( is.null(folder)) stop("Missing required parameter folder. Aborting.")
    if ( is.null(files)) stop("Missing required parameter files Aborting.")
    if ( is.null(timespan)) stop("Missing required parameter timespan. Aborting.")

  }else{ # output is supplied via data list-object
    if ( is.null(data)) stop("In data supply mode, data cannot be NULL. Aborting.")
    if (!is.null(folder)) stop("In data supply mode, folder cannot be supplied. Aborting.")
    if (!is.null(files)) stop("In data supply mode, files cannot be supplied. Aborting.")
    if (!is.null(timespan)) stop("In data supply mode, timespan cannot be supplied. Aborting.")
  }

  # define default minimum tree cover for forest / woodland / savanna
  min_tree_cover <- list("boreal forest" = 0.6, "temperate forest" = 0.6,
                         "temperate woodland" = 0.3, "temperate savannah" = 0.1,
                         "tropical forest" = 0.6, "tropical woodland" = 0.3,
                         "tropical savannah" = 0.1)

  # replace default values by values defined in tree_cover_thresholds
  # parameter
  overwrite <- match(names(tree_cover_thresholds), names(min_tree_cover))
  if (any(is.na(overwite))) {
    stop(paste0(
      names(tree_cover_thresholds)[which(is.na(overwrite))],
      " is not valid. Please use a name of: ",
      paste0(names(min_tree_cover), collapse = ", ")
    ))
  }
  min_tree_cover[overwrite] <- tree_cover_thresholds

  # test if forest threshold is always > woodland threshold > savanna threshold
  if (min_tree_cover[["temperate forest"]] <=
        min_tree_cover[["temperate woodland"]] |
      min_tree_cover[["temperate woodland"]] <=
        min_tree_cover[["temperate savannah"]] |
      min_tree_cover[["tropical woodland"]] <=
        min_tree_cover[["tropical savannah"]] |
      min_tree_cover[["tropical forest"]] <=
        min_tree_cover[["tropical woodland"]]) {
    stop(paste0("Tree cover threshold for forest are not always higher than",
                "tree cover thresholds for woodland and savannah. Aborting."))
  }

  if (!readOutput) {
    #process grid
    lpjml_grid  <- rbind(data$lon,data$lat)

    # process foliage projected cover (fpc)
    fpc <- data$fpc
    di <- dim(fpc)
    npft <- di[2] - 1
    ncell <- di[1]

    if (savannaProxy == "vegc") {
      # process vegetation carbon output
      vegc <- data$vegc
    }else if (savannaProxy == "natLAI") {
      # process pft_lai input
      pft_lai <- data$pft_lai
    }else{
      stop(paste0("Unknown setting (",savannaProxy,") for savannaProxy. Enter either 'natLAI' or 'vegc'. Aborting."))
    }

    # process temperature input
    temp <- data$temp

    # process elevation input
    if (montaneArcticProxy == "elevation") {
      elevation <- data$elevation
    }
  }else{ # read in output

    grid_ending <- tail(strsplit(files$grid,".", fixed = T)[[1]], n = 1)
    if (grid_ending %in% c("bin","clm","raw")) {
      grid <- lpjmliotools::autoReadMetaOutput(metaFile = paste0(folder,"/",files$grid,".json"))
      ncell <- length(grid)/2
      lon   <- grid[c(1:ncell)*2 - 1]
      lat   <- grid[c(1:ncell)*2]
    }else if (grid_ending %in% c("nc","cdf")) {
      grid <- readGridInputBin(inFile = lpjGridInput, headersize = lpjGridHeaderSize, ncells = lpjCells)
      ncell <- lpjCells
      lon   <- grid$lon
      lat   <- grid$lat
    }else{
      stop(paste0("Unknown file ending (",grid_ending,"). Aborting."))
    }

    lpjml_grid <- rbind(lon,lat)

    fpc_ending <- tail(strsplit(files$fpc,".", fixed = T)[[1]], n = 1)
    if (fpc_ending %in% c("bin","clm","raw")) {
      fpc <- apply(lpjmliotools::autoReadMetaOutput(
                            metaFile = paste0(folder,"/",files$fpc,".json"),
                            getyearstart = timespan[1], getyearstop = timespan[2]),
                            c(1,2),mean)
    }else if (fpc_ending %in% c("nc","cdf")) {
      fpc <- lpjmliotools::netcdfCFT2lpjarray(ncInFile = files$fpc, var = "FPC", lon = lon, lat = lat)
    }else{
      stop(paste0("Unknown file ending (",fpc_ending,"). Aborting."))
    }
    di <- dim(fpc)
    npft <- di[2] - 1

    if (savannaProxy == "vegc") {
      vegc_ending <- tail(strsplit(files$vegc,".", fixed = T)[[1]], n = 1)
      if (vegc_ending %in% c("bin","clm","raw")) {
        vegc <- apply(lpjmliotools::autoReadMetaOutput(
                    metaFile = paste0(folder,"/",files$vegc,".json"),
                    getyearstart = timespan[1], getyearstop = timespan[2]),
                      c(1,2),mean)
      }else if (vegc_ending %in% c("nc","cdf")) {
        vegc <- lpjmliotools::netcdfCFT2lpjarray(ncInFile = files$vegc, var = "VegC", lon = lon, lat = lat)
      }else{
        stop(paste0("Unknown file ending (",vegc_ending,"). Aborting."))
      }
    }

    if (montaneArcticProxy == "elevation") {
        elevation <- lpjmliotools::autoReadInput(inFile = elevationInput)[1,]
        #plotGlobalWlin(data = elevation,file = "/home/stenzel/elevation.png",title = "",max = 6000,min=-100,legYes = T,legendtitle = "",eps = F)
    }

    temp_ending <- tail(strsplit(files$temp,".", fixed = T)[[1]], n = 1)
    if (temp_ending %in% c("bin","clm","raw")) {
      temp <- apply(lpjmliotools::autoReadMetaOutput(
                    metaFile = paste0(folder,"/",files$temp,".json"),
                    getyearstart = timespan[1], getyearstop = timespan[2]),
                    c(1,2),mean)
    }else if (temp_ending %in% c("nc","cdf")) {
      temp <- lpjmliotools::netcdfCFT2lpjarray(ncInFile = files$temp, var = "temp", lon = lon, lat = lat)
    }else{
      stop(paste0("Unknown file ending (",temp_ending,"). Aborting."))
    }

    if (savannaProxy == "natLAI") {
      pft_lai_ending <- tail(strsplit(files$pft_lai,".", fixed = T)[[1]], n = 1)
      if (pft_lai_ending %in% c("bin","clm","raw")) {
        pft_lai <- apply(lpjmliotools::autoReadMetaOutput(
                         metaFile = paste0(folder,"/",files$pft_lai,".json"),
                         getyearstart = timespan[1], getyearstop = timespan[2]),
                         c(1,2),mean)
      }else if (pft_lai_ending %in% c("nc","cdf")) {
        pft_lai <- lpjmliotools::netcdfCFT2lpjarray(ncInFile = files$pft_lai, var = "LAI", lon = lon, lat = lat)
      }else{
        stop(paste0("Unknown file ending (",pft_lai_ending,"). Aborting."))
      }
    }

  }

  # biome_names after biome classification in Ostberg et al. 2013
  # (https://doi.org/10.5194/esd-4-347-2013), Ostberg et al 2015
  # (https://doi.org/10.1088/1748-9326/10/4/044011) and Gerten et al. 2020
  # (https://doi.org/10.1038/s41893-019-0465-1)
  # biome names

  biome_mapping <- system.file("extdata", "biomes.csv",
                              package = "pbfunctions") %>%
                   readr::read_csv2()
  biome_names <- biome_mapping$id
  names(biome_names) <- biome_mapping$name

  if (npft == 9) {
    fpc_names <- c("natvegfrac", 												    #1
                   "Tropical Broadleaved Evergreen Tree",		#2
                   "Tropical Broadleaved Raingreen Tree",		#3
                   "Temperate Needleleaved Evergreen Tree",	#4
                   "Temperate Broadleaved Evergreen Tree",	#5
                   "Temperate Broadleaved Summergreen Tree",#6
                   "Boreal Needleleaved Evergreen Tree",		#7
                   "Boreal Broadleaved Summergreen Tree",   #8 actually Boreal Summergreen Tree (broadleaved and needleleaved, includes larchs)
                   "Temperate C3 Grass",										#9
                   "Tropical C4 Grass"											#10
    )
    fpc_temperate <- c(4,5,6)
    fpc_tropical  <- c(2,3)
    fpc_boreal    <- c(7,8)
    fpc_needle    <- c(4,7,8)
    fpc_grass     <- c(9,10)
  }else if (npft == 11) {
    fpc_names <- c("natvegfrac", 														#1
                   "Tropical Broadleaved Evergreen Tree",		#2
                   "Tropical Broadleaved Raingreen Tree",		#3
                   "Temperate Needleleaved Evergreen Tree",	#4
                   "Temperate Broadleaved Evergreen Tree",	#5
                   "Temperate Broadleaved Summergreen Tree",#6
                   "Boreal Needleleaved Evergreen Tree",		#7
                   "Boreal Broadleaved Summergreen Tree",   #8
                   "Boreal Needleleaved Summergreen Tree",  #9
                   "Tropical C4 Grass",										  #10
                   "Temperate C3 Grass",										#11
                   "Polar C3 Grass"											    #12
    )
    fpc_temperate <- c(4,5,6)
    fpc_tropical  <- c(2,3)
    fpc_boreal    <- c(7,8,9)
    fpc_needle    <- c(4,7,9)
    fpc_grass     <- c(10,11,12)
  }else{stop(paste("Unknown number of pfts:",npft))}
  fpc_trees     <- sort(unique(c(fpc_temperate, fpc_tropical, fpc_boreal, fpc_needle)))

  # indices for pft subsets
  fpc_tropical_trees <- c("Tropical Broadleaved Evergreen Tree",
                          "Tropical Broadleaved Raingreen Tree")
  fpc_temperate_trees <- c("Temperate Needleleaved Evergreen Tree",
                           "Temperate Broadleaved Evergreen Tree",
                           "Temperate Broadleaved Summergreen Tree")
  fpc_evergreen_trees <- c("Tropical Broadleaved Evergreen Tree",
                           "Temperate Needleleaved Evergreen Tree",
                           "Temperate Broadleaved Evergreen Tree",
                           "Boreal Needleleaved Evergreen Tree") # was "Temperate Broadleaved Summergreen Tree") #FS?
  if (npft == 9) {
    fpc_boreal_trees <- c("Boreal Needleleaved Evergreen Tree",
                          "Boreal Broadleaved Summergreen Tree")
    fpc_needle_trees <- c("Temperate Needleleaved Evergreen Tree",
                          "Boreal Needleleaved Evergreen Tree")
    fpc_grass <- c("Tropical C4 Grass",
                   "Temperate C3 Grass")
  }else if (npft == 11) {
    fpc_boreal_trees <- c("Boreal Needleleaved Evergreen Tree",
                          "Boreal Broadleaved Summergreen Tree",
                          "Boreal Needleleaved Summergreen Tree")
    fpc_needle_trees <- c("Temperate Needleleaved Evergreen Tree",
                          "Boreal Needleleaved Evergreen Tree",
                          "Boreal Needleleaved Summergreen Tree")
    fpc_grass <- c("Tropical C4 Grass",
                   "Temperate C3 Grass",
                   "Polar C3 Grass")
  }else{stop(paste("Unknown number of pfts:",npft))}

  fpc_trees <- c(fpc_tropical_trees,
                 fpc_temperate_trees,
                 fpc_boreal_trees)

  #process grid
  dim(lpjml_grid) <- c(coordinate = 2, cell = ncell)
  dimnames(lpjml_grid) <- list(coordinate = c("lon", "lat"),
                               cell = seq_len(ncell))

  # process foliage projected cover (fpc)
  dim(fpc) <- c(cell = ncell,
                band = npft + 1)
  dimnames(fpc) <- list(cell = seq_len(ncell),
                        band = fpc_names)

  if (savannaProxy == "vegc") {
    # process vegetation carbon output
    dim(vegc) <- c(cell = ncell)
    dimnames(vegc) <- list(cell = seq_len(ncell))
  }

  # process temperature input
  dim(temp) <- c(cell = ncell)
  dimnames(temp) <- list(cell = seq_len(ncell))

  if (savannaProxy == "natLAI") {
    # process pft_lai input
    di2 <- dim(pft_lai)
    dim(pft_lai) <- c(cell = di2[1],
                      band = di2[2])
    dimnames(pft_lai) <- list(cell = seq_len(di2[1]),
                              band = c(fpc_names[2:(npft + 1)],(npft + 1):di2[2]))
  }

  # latitudes (same dimension for vectorized biome classification)
  latitudes <- array(
    subset_array(lpjml_grid, list(coordinate = "lat")),
    dim = c(ncell)
  )

  fpc_tree_total <- apply(
    subset_array(fpc, list(band = fpc_trees)),
    c("cell"),
    sum,
    na.rm = TRUE
  )

  fpc_tree_tropical <- apply(
    subset_array(fpc, list(band = fpc_tropical_trees)),
    c("cell"),
    sum,
    na.rm = TRUE
  )
  fpc_tree_temperate <- apply(
    subset_array(fpc, list(band = fpc_temperate_trees)),
    c("cell"),
    sum,
    na.rm = TRUE
  )
  fpc_tree_boreal <- apply(
    subset_array(fpc, list(band = fpc_boreal_trees)),
    c("cell"),
    sum,
    na.rm = TRUE
  )
  fpc_tree_needle <- apply(
    subset_array(fpc, list(band = fpc_needle_trees)),
    c("cell"),
    sum,
    na.rm = TRUE
  )
  fpc_tree_evergreen <- apply(
    subset_array(fpc, list(band = fpc_evergreen_trees)),
    c("cell"),
    sum,
    na.rm = TRUE
  )
  fpc_grass_total <- apply(
    subset_array(fpc, list(band = fpc_grass)),
    c("cell"),
    sum,
    na.rm = TRUE
  )
  fpc_total <- apply(
    subset_array(fpc,
                           list(band = fpc_names[fpc_names != "natvegfrac"])),
    c("cell"),
    sum,
    na.rm = TRUE
  )
  max_share_trees <- apply(
    subset_array(fpc, list(band = fpc_trees)),
    c("cell"),
    max,
    na.rm = TRUE
  )
  # max_share <- apply(
  #   subset_array(fpc,
  #                          list(band = fpc_names[fpc_names != "natvegfrac"])),
  #   c("cell"),
  #   max,
  #   na.rm = TRUE
  # )
  fpc_tree_broadleaf <- fpc_tree_total - fpc_tree_needle

  if (savannaProxy == "natLAI") {
    #prepare natLAI array
    natLAI <- rowSums( pft_lai[,1:npft] * fpc[,2:(npft + 1)] * fpc[,1] )
  }

  # use vegc 7500 gC/m2 or natLAI 6 as proxy threshold for forest/savannah "boundary"
  if (savannaProxy == "vegc") {
    is_tropical_proxy <- vegc >= vegc_threshold
    is_savannah_proxy <- vegc < vegc_threshold
  } else if (savannaProxy == "natLAI") {
    is_tropical_proxy <- natLAI >= lai_threshold
    is_savannah_proxy <- natLAI < lai_threshold
  } else {
    stop(paste0("Unknown parameter savannaProxy = ",savannaProxy))
  }


  # Desert
  is_desert <- {
    fpc_total <= 0.05 &
      temp >= 0 #-2
  }

  # montane (for classification of montane grassland)
  is_montane <- {
    elevation > elevation_threshold
  }

  # high latitude (for classification of montane grassland)
  is_high_latitude <- {
    abs(latitudes) > latitude_threshold
  }

  # FORESTS ------------------------------------------------------------------ #
  is_boreal_forest <- {
    fpc_tree_total >= min_tree_cover[["boreal forest"]]
  }
  is_temperate_forest <- {
    fpc_tree_total >= min_tree_cover[["temperate forest"]]
  }
  is_tropical_forest <- {
    fpc_tree_total >= min_tree_cover[["tropical forest"]]
  }
  # Boreal Evergreen
  is_boreal_evergreen <- {
    is_boreal_forest &
      subset_array(fpc,
                             list(band = "Boreal Needleleaved Evergreen Tree")) ==
      max_share_trees &
      fpc_tree_broadleaf < (0.4 * fpc_tree_total)
  }

  if (npft == 9) {
    # Boreal Deciduous
    is_boreal_deciduous <- {
      is_boreal_forest &
        (subset_array(fpc,
                                list(band = "Boreal Broadleaved Summergreen Tree")) == # nolint
           max_share_trees) &
        fpc_tree_evergreen < (0.4 * fpc_tree_total)
    }
  }else if (npft == 11) {
    # Boreal Deciduous
    is_boreal_deciduous <- {
      is_boreal_forest &
        (subset_array(fpc,
                                list(band = "Boreal Broadleaved Summergreen Tree")) == # nolint
           max_share_trees |
           subset_array(fpc,
                                  list(band = "Boreal Needleleaved Summergreen Tree")) == # nolint
           max_share_trees) &
        fpc_tree_evergreen < (0.4 * fpc_tree_total)
    }
  }else{stop(paste("Unknown number of pfts:",npft))}

  # Temperate Coniferous Forest
  is_temperate_coniferous <- {
    is_temperate_forest &
      subset_array(fpc,
                             list(band = "Temperate Needleleaved Evergreen Tree")) == # nolint
      max_share_trees &
      fpc_tree_broadleaf < (0.4 * fpc_tree_total)
  }
  # Temperate Broadleaved Evergreen Forest
  is_temperate_broadleaved_evergreen <- { # nolint
    is_temperate_forest &
      subset_array(fpc,
                             list(band = "Temperate Broadleaved Evergreen Tree")) == # nolint
      max_share_trees &
      fpc_tree_tropical < (0.4 * fpc_tree_total) &
      fpc_tree_needle < (0.4 * fpc_tree_total)
  }
  # Temperate Broadleaved Deciduous Forest
  is_temperate_broadleaved_deciduous <- { # nolint
    is_temperate_forest &
      subset_array(fpc,
                             list(band = "Temperate Broadleaved Summergreen Tree")) == # nolint
      max_share_trees &
      fpc_tree_tropical < (0.4 * fpc_tree_total) &
      fpc_tree_needle < (0.4 * fpc_tree_total)
  }
  # Tropical Rainforest
  is_tropical_evergreen <- {
    is_tropical_forest &
      subset_array(fpc,
                             list(band = "Tropical Broadleaved Evergreen Tree")) == # nolint
      max_share_trees &
      (fpc_tree_boreal + fpc_tree_temperate) < (0.4 * fpc_tree_total) &
      is_tropical_proxy
  }
  # Tropical Seasonal & Deciduous Forest
  is_tropical_raingreen <- {
    is_tropical_forest &
      (subset_array(fpc,
                              list(band = "Tropical Broadleaved Raingreen Tree")) == # nolint
         max_share_trees) &
      (fpc_tree_boreal + fpc_tree_temperate) < (0.4 * fpc_tree_total) &
      is_tropical_proxy
  }
  # Warm Woody Savanna, Woodland & Shrubland
  is_tropical_forest_savannah <- {
    is_tropical_forest &
      (subset_array(fpc,
                              list(band = "Tropical Broadleaved Evergreen Tree")) == # nolint
         max_share_trees |
         subset_array(fpc,
                                list(band = "Tropical Broadleaved Raingreen Tree")) == # nolint
         max_share_trees) &
      (fpc_tree_boreal + fpc_tree_temperate) < (0.4 * fpc_tree_total) &
      is_savannah_proxy
  }
  is_mixed_forest <- {
    is_temperate_forest &
      !is_boreal_evergreen &
      !is_boreal_deciduous &
      !is_temperate_coniferous &
      !is_temperate_broadleaved_evergreen &
      !is_temperate_broadleaved_deciduous &
      !is_tropical_evergreen &
      !is_tropical_raingreen &
      !is_tropical_forest_savannah
  }

  # WOODY SAVANNAH ----------------------------------------------------------- #

  # Temperate Woody Savanna, Woodland & Shrubland
  is_temperate_woody_savannah <- {
    fpc_tree_total <= min_tree_cover[["temperate forest"]] &
    fpc_tree_total >= min_tree_cover[["temperate woodland"]] &
      subset_array(fpc, list(band = "Temperate C3 Grass")) >
      subset_array(fpc, list(band = "Tropical C4 Grass")) &
      temp >= 0 #-2 &
    #latitudes < 55
  }
  # Warm Woody Savanna, Woodland & Shrubland
  is_tropical_woody_savannah <- {
    fpc_tree_total <= min_tree_cover[["tropical forest"]] &
    fpc_tree_total >= min_tree_cover[["tropical woodland"]] &
      subset_array(fpc, list(band = "Temperate C3 Grass")) <
      subset_array(fpc, list(band = "Tropical C4 Grass"))
  }

  # OPEN SHRUBLAND / SAVANNAHS ----------------------------------------------- #

  # Temperate Savanna & Open Shrubland
  is_temperate_shrubland <- {
    fpc_tree_total <= min_tree_cover[["temperate woodland"]] &
    fpc_tree_total >= min_tree_cover[["temperate savannah"]] &
      subset_array(fpc, list(band = "Temperate C3 Grass")) >
      subset_array(fpc, list(band = "Tropical C4 Grass")) &
      temp >= 0 #-2 &
    #latitudes < 55
  }
  # Warm Savanna & Open Shrubland
  is_tropical_shrubland <- {
    fpc_tree_total <= min_tree_cover[["tropical woodland"]] &
    fpc_tree_total >= min_tree_cover[["tropical savannah"]] &
      subset_array(fpc, list(band = "Temperate C3 Grass")) <
      subset_array(fpc, list(band = "Tropical C4 Grass")) &
      temp >= 0 #-2
  }

  # GRASSLAND ---------------------------------------------------------------- #

  # Temperate Savanna & Open Shrubland
  is_temperate_grassland <- {
    fpc_total > 0.05 &
    fpc_tree_total <= min_tree_cover[["temperate savannah"]] &
      subset_array(fpc, list(band = "Temperate C3 Grass")) >
      subset_array(fpc, list(band = "Tropical C4 Grass")) &
      temp >= 0 #-2 &
    #latitudes < 55
  }
  # Warm Savanna & Open Shrubland
  is_tropical_grassland <- {
    fpc_total > 0.05 &
    fpc_tree_total <= min_tree_cover[["tropical savannah"]] &
      subset_array(fpc, list(band = "Temperate C3 Grass")) <
      subset_array(fpc, list(band = "Tropical C4 Grass")) &
      temp >= 0 #-2
  }

  # Arctic Tundra ------------------------------------------------------------ #
  is_arctic_tundra <- {
      (!is_boreal_evergreen &
      !is_boreal_deciduous &
      !is_temperate_forest &
      (temp < 0 |
      subset_array(fpc, list(band = "Temperate C3 Grass")) ==
      subset_array(fpc, list(band = "Tropical C4 Grass"))) &
      fpc_total > 0.05) |

      (temp < 0 & fpc_total < 0.05)
  }

  # Rocks and Ice
  is_rocks_and_ice <- {
    fpc_total == 0 &
      temp < 0 #-2
  }
  # Water body
  is_water <- {
    subset_array(fpc, list(band = "natvegfrac")) == 0
  }

  # CLASSIFY BIOMES ---------------------------------------------------------- #

  # initiate biome_class array
  biome_class <- array(NA, dim = c(ncell), dimnames = dimnames(fpc_total))

  biome_class[is_desert] <- biome_names["Desert"]

  # forests
  biome_class[is_boreal_evergreen] <- biome_names["Boreal Evergreen Forest"]
  biome_class[is_boreal_deciduous] <- biome_names["Boreal Deciduous Forest"]
  biome_class[is_temperate_coniferous] <- biome_names["Temperate Coniferous Forest"] # nolint
  biome_class[is_temperate_broadleaved_evergreen] <- biome_names["Temperate Broadleaved Evergreen Forest"] # nolint
  biome_class[is_temperate_broadleaved_deciduous] <- biome_names["Temperate Broadleaved Deciduous Forest"] # nolint
  biome_class[is_tropical_evergreen] <- biome_names["Tropical Rainforest"]
  biome_class[is_tropical_raingreen] <- biome_names["Tropical Seasonal & Deciduous Forest"] # nolint
  biome_class[is_tropical_forest_savannah] <- biome_names["Warm Woody Savanna, Woodland & Shrubland"] # nolint
  biome_class[is_mixed_forest] <- biome_names["Mixed Forest"]

  # woody savannah
  biome_class[is_temperate_woody_savannah] <- biome_names["Temperate Woody Savanna, Woodland & Shrubland"] # nolint
  biome_class[is_tropical_woody_savannah] <- biome_names["Warm Woody Savanna, Woodland & Shrubland"] # nolint

  # open shrubland / savannah
  biome_class[is_temperate_shrubland] <- biome_names["Temperate Savanna & Open Shrubland"] # nolint
  biome_class[is_tropical_shrubland] <- biome_names["Warm Savanna & Open Shrubland"] # nolint

  # grassland
  biome_class[is_temperate_grassland] <- biome_names["Temperate Grassland"]
  biome_class[is_tropical_grassland] <- biome_names["Warm Grassland"]

  biome_class[is_arctic_tundra] <- biome_names["Arctic Tundra"]
  if (montaneArcticProxy == "elevation") {
    biome_class[biome_class == biome_names["Arctic Tundra"] & is_montane] <- biome_names["Montane Grassland"]
  }else if (montaneArcticProxy == "latitude") {
    biome_class[biome_class == biome_names["Arctic Tundra"] & !is_high_latitude] <- biome_names["Montane Grassland"]
  }else {
    stop(paste0("Unknown value (",montaneArcticProxy,") for parameter montaneArcticProxy. Use 'elevation' or 'latitude'. Aborting."))
  }

  # other
  biome_class[is_rocks_and_ice] <- biome_names["Rocks and Ice"]
  biome_class[is_water] <- biome_names["Water"]

  return(list(biome_id = biome_class, biome_names = names(biome_names)))
}
