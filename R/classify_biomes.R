# written by Fabian Stenzel, based on work by Sebastian Ostberg, Johanna Braun, Jannes Breier
# 2022 - stenzel@pik-potsdam.de

#requirements:
require(lpjmliotools) # in at least version 0.2.17

#' Classify biomes
#'
#' Classify biomes based on foliage protected cover (FPC) and temperature
#' LPJmL output plus either vegetation carbon or pft_lai depending on
#' the savanna_proxy option and elevation if montane_arctic_proxy requires this
#'
#' @param folder to read default outputs from (fpc, grid, vegc, pft_lai, temp)
#' @param diff_output_files optional list for specification of output file names
#'        differing from default, which is list(grid = "grid.bin", fpc = "fpc.bin",
#'        vegc = "vegc.bin", pft_lai = "pft_lai.bin", temp = "temp.bin")
#' @param input_files list containing additional input (!) files, not in the
#'        folder, e.g. if temp was not written out:
#'        list(grid=..., temp = ..., elevation = ...)
#' @param file_ending replace default file ending. default: ".bin"
#' @param timespan as c(startyear,stopyear) to use for averaging outputs over

#' @param savanna_proxy "vegc", "natLAI" or NULL. Use vegetation carbon or LAI
#'        in natural vegetation as a proxy threshold to distinguish forests and
#'        savannahs. Set to NULL if no savanna proxy should be used
#'        - default: "natLAI"
#' @param montane_arctic_proxy "elevation" or "latitude". Use elevation or latitude
#'        as a proxy threshold to distinguish arctic tundra and montane grassland
#' @param elevation_threshold threshold in m above which ArcticTundra is
#'        classified as Montane Grassland, if montane_arctic_proxy is set to
#'        elevation - default: 1000
#' @param latitude_threshold threshold in degrees, south of which ArcticTundra is
#'        classified as Montane Grassland, if montane_arctic_proxy is set to
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
#'        "temperate savanna" = 0.1
#'        "tropical forest" = 0.6
#'        "tropical woodland" = 0.3
#'        "tropical savanna" = 0.1
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
#' classify_biomes(timespan = c(1982,2011)
#'      folder = "/p/projects/open/Fabian/runs/Gamma/output/historic_gamma/",
#'      files = list(grid = "grid.bin", fpc = "fpc.bin", vegc = "vegc.bin",
#'      pft_lai = "pft_lai.bin", temp = "temp.bin"))
#' }
#'
#' @export
classify_biomes <- function(folder = NULL, input_files = NULL, timespan = NULL,
                   diff_output_files = NULL, savanna_proxy = "natLAI", file_ending = ".bin",
                   montane_arctic_proxy = "elevation", lai_threshold = 6,
                   elevation_threshold = 1000, vegc_threshold = 7500,
                   latitude_threshold = 55, tree_cover_thresholds = list(),
                   lpjGridHeaderSize = 43, lpjCells = 67420) {

  require(lpjmliotools)
  output_files = list(grid = "grid.bin", fpc = "fpc.bin", vegc = "vegc.bin",
                      pft_lai = "pft_lai.bin",  temp = "temp.bin")

  # replace file_ending
  output_files <- gsub(".bin", file_ending, output_files)

  if (!is.null(diff_output_files)) {
    overwrite <- match(names(diff_output_files), names(output_files))
    output_files[c(overwrite)] <- diff_output_files
  }

  if ( is.null(folder) ) stop("Missing required parameter folder. Aborting.")
  # todo add consistence check for required inputs
  #if ( is.null(files) ) stop("Missing required parameter files Aborting.")
  if ( is.null(timespan) ) stop("Missing required parameter timespan. Aborting.")


  # define default minimum tree cover for forest / woodland / savanna
  min_tree_cover <- list("boreal forest" = 0.6, "temperate forest" = 0.6,
                         "temperate woodland" = 0.3, "temperate savanna" = 0.1,
                         "tropical forest" = 0.6, "tropical woodland" = 0.3,
                         "tropical savanna" = 0.1)

  # replace default values by values defined in tree_cover_thresholds
  # parameter
  overwrite <- match(names(tree_cover_thresholds), names(min_tree_cover))
  if (any(is.na(overwrite))) {
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
        min_tree_cover[["temperate savanna"]] |
      min_tree_cover[["tropical woodland"]] <=
        min_tree_cover[["tropical savanna"]] |
      min_tree_cover[["tropical forest"]] <=
        min_tree_cover[["tropical woodland"]]) {
    stop(paste0("Tree cover threshold for forest are not always higher than",
                "tree cover thresholds for woodland and savanna. Aborting."))
  }

  # test if savanna proxy is valid
  match.arg(savanna_proxy, c("vegc", "natLAI"))

  if (file.exists(output_files$grid)) {
    grid_ending <- tail(strsplit(output_files$grid,".", fixed = T)[[1]], n = 1)
    if (grid_ending %in% c("bin","clm","raw")) {
      grid <- lpjmliotools::autoReadMetaOutput(metaFile = paste0(folder,"/",files$grid,".json"))
      ncell <- length(grid)/2
      lon   <- grid[c(1:ncell)*2 - 1]
      lat   <- grid[c(1:ncell)*2]
    }else if (grid_ending %in% c("nc","cdf")) {
      print("Reading of netcdf output is still preliminary. Please specify LPJmL grid input.")
      grid <- readGridInputBin(inFile = input_files$grid, headersize = lpjGridHeaderSize, ncells = lpjCells)

    }else{
      stop(paste0("Unknown file ending (",grid_ending,"). Aborting."))
    }
  }else{
    stop(paste0("Output file ",output_files$grid, " does not exist. Make sure
                 the specified input folder is correct. If your file names
                 differ from the default, please use diff_output_files to
                 specify them. "))

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

    if (!is.null(savanna_proxy)) {
      if (savanna_proxy == "vegc") {
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
      } else if (savanna_proxy == "natLAI") {
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

    if (montane_arctic_proxy == "elevation") {
        elevation <- lpjmliotools::autoReadInput(inFile = input_files$elevation)[1,]
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



  # biome_names after biome classification in Ostberg et al. 2013
  # (https://doi.org/10.5194/esd-4-347-2013), Ostberg et al 2015
  # (https://doi.org/10.1088/1748-9326/10/4/044011) and Gerten et al. 2020
  # (https://doi.org/10.1038/s41893-019-0465-1)
  # biome names

  biome_mapping <- system.file("extdata",
                               "biomes.csv",
                               package = "pbfunctions") %>%
                   readr::read_delim(col_types = readr::cols(), delim = ";")
  biome_names <- biome_mapping$id
  names(biome_names) <- biome_mapping$name


  pft_categories <- system.file("extdata",
                                "pft_categories.csv",
                                package = "pbfunctions") %>%
    read_pft_categories() %>%
    dplyr::filter(npft_proxy == npft)

  fpc_names <- dplyr::filter(pft_categories, category == "natural")$pft

  # indices (when estimation only via npft possible) or names for pft subsets
  fpc_temperate_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & zone == "temperate" & category == "natural"
  ) %>% {
      if (any(is.na(.$npft_proxy))) .$pft else .$lpjml_index + 1
    }

  fpc_tropical_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & zone == "temperate" & category == "natural"
  ) %>% {
    if (any(is.na(.$npft_proxy))) .$pft else .$lpjml_index + 1
  }

  fpc_boreal_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & zone == "boreal" & category == "natural"
  ) %>% {
    if (any(is.na(.$npft_proxy))) .$pft else .$lpjml_index + 1
  }

  fpc_needle_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & category == "needle"
  ) %>% {
    if (any(is.na(.$npft_proxy))) .$pft else .$lpjml_index + 1
  }

  fpc_grass <- dplyr::filter(
    pft_categories,
    type == "grass" & category == "natural"
  ) %>% {
    if (any(is.na(.$npft_proxy))) .$pft else .$lpjml_index + 1
  }

  fpc_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & category == "natural"
  ) %>% {
    if (any(is.na(.$npft_proxy))) .$pft else .$lpjml_index + 1
  }

  #process grid
  dim(lpjml_grid) <- c(coordinate = 2, cell = ncell)
  dimnames(lpjml_grid) <- list(coordinate = c("lon", "lat"),
                               cell = seq_len(ncell))

  # process foliage projected cover (fpc)
  dim(fpc) <- c(cell = ncell,
                band = npft + 1)
  dimnames(fpc) <- list(cell = seq_len(ncell),
                        band = fpc_names)

  if (!is.null(savanna_proxy)) {
    if (savanna_proxy == "vegc") {
      # process vegetation carbon output
      dim(vegc) <- c(cell = ncell)
      dimnames(vegc) <- list(cell = seq_len(ncell))
    } else if (savanna_proxy == "natLAI") {
      # process pft_lai input
      di2 <- dim(pft_lai)
      dim(pft_lai) <- c(cell = di2[1],
                      band = di2[2])
      dimnames(pft_lai) <- list(cell = seq_len(di2[1]),
                              band = c(fpc_names[2:(npft + 1)],(npft + 1):di2[2]))
    }
  }

  # process temperature input
  dim(temp) <- c(cell = ncell)
  dimnames(temp) <- list(cell = seq_len(ncell))

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

  # use vegc 7500 gC/m2 or natLAI 6 as proxy threshold for forest/savannah "boundary"
  if (!is.null(savanna_proxy)) {
    if (savanna_proxy == "vegc") {
      is_tropical_proxy <- vegc >= vegc_threshold
      is_savannah_proxy <- vegc < vegc_threshold
    } else if (savanna_proxy == "natLAI") {
      #prepare natLAI array
      natLAI <- rowSums( pft_lai[,1:npft] * fpc[,2:(npft + 1)] * fpc[,1] )
      is_tropical_proxy <- natLAI >= lai_threshold
      is_savannah_proxy <- natLAI < lai_threshold
    }
  } else {
    is_tropical_proxy <- rep(TRUE, ncell)
    is_savannah_proxy <- rep(FALSE, ncell)
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
  is_tropical_forest_savanna <- {
    is_tropical_forest &
      (subset_array(fpc,
                              list(band = "Tropical Broadleaved Evergreen Tree")) == # nolint
         max_share_trees |
         subset_array(fpc,
                                list(band = "Tropical Broadleaved Raingreen Tree")) == # nolint
         max_share_trees) &
      (fpc_tree_boreal + fpc_tree_temperate) < (0.4 * fpc_tree_total) &
      is_savanna_proxy
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
      !is_tropical_forest_savanna
  }

  # WOODY savanna ----------------------------------------------------------- #

  # Temperate Woody Savanna, Woodland & Shrubland
  is_temperate_woody_savanna <- {
    fpc_tree_total <= min_tree_cover[["temperate forest"]] &
    fpc_tree_total >= min_tree_cover[["temperate woodland"]] &
      subset_array(fpc, list(band = "Temperate C3 Grass")) >
      subset_array(fpc, list(band = "Tropical C4 Grass")) &
      temp >= 0 #-2 &
    #latitudes < 55
  }
  # Warm Woody Savanna, Woodland & Shrubland
  is_tropical_woody_savanna <- {
    fpc_tree_total <= min_tree_cover[["tropical forest"]] &
    fpc_tree_total >= min_tree_cover[["tropical woodland"]] &
      subset_array(fpc, list(band = "Temperate C3 Grass")) <
      subset_array(fpc, list(band = "Tropical C4 Grass"))
  }

  # OPEN SHRUBLAND / SAVANNAS ----------------------------------------------- #

  # Temperate Savanna & Open Shrubland
  is_temperate_shrubland <- {
    fpc_tree_total <= min_tree_cover[["temperate woodland"]] &
    fpc_tree_total >= min_tree_cover[["temperate savanna"]] &
      subset_array(fpc, list(band = "Temperate C3 Grass")) >
      subset_array(fpc, list(band = "Tropical C4 Grass")) &
      temp >= 0 #-2 &
    #latitudes < 55
  }
  # Warm Savanna & Open Shrubland
  is_tropical_shrubland <- {
    fpc_tree_total <= min_tree_cover[["tropical woodland"]] &
    fpc_tree_total >= min_tree_cover[["tropical savanna"]] &
      subset_array(fpc, list(band = "Temperate C3 Grass")) <
      subset_array(fpc, list(band = "Tropical C4 Grass")) &
      temp >= 0 #-2
  }

  # GRASSLAND ---------------------------------------------------------------- #

  # Temperate Savanna & Open Shrubland
  is_temperate_grassland <- {
    fpc_total > 0.05 &
    fpc_tree_total <= min_tree_cover[["temperate savanna"]] &
      subset_array(fpc, list(band = "Temperate C3 Grass")) >
      subset_array(fpc, list(band = "Tropical C4 Grass")) &
      temp >= 0 #-2 &
    #latitudes < 55
  }
  # Warm Savanna & Open Shrubland
  is_tropical_grassland <- {
    fpc_total > 0.05 &
    fpc_tree_total <= min_tree_cover[["tropical savanna"]] &
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
  biome_class[is_tropical_forest_savanna] <- biome_names["Warm Woody Savanna, Woodland & Shrubland"] # nolint
  biome_class[is_mixed_forest] <- biome_names["Mixed Forest"]

  # woody savanna
  biome_class[is_temperate_woody_savanna] <- biome_names["Temperate Woody Savanna, Woodland & Shrubland"] # nolint
  biome_class[is_tropical_woody_savanna] <- biome_names["Warm Woody Savanna, Woodland & Shrubland"] # nolint

  # open shrubland / savanna
  biome_class[is_temperate_shrubland] <- biome_names["Temperate Savanna & Open Shrubland"] # nolint
  biome_class[is_tropical_shrubland] <- biome_names["Warm Savanna & Open Shrubland"] # nolint

  # grassland
  biome_class[is_temperate_grassland] <- biome_names["Temperate Grassland"]
  biome_class[is_tropical_grassland] <- biome_names["Warm Grassland"]

  biome_class[is_arctic_tundra] <- biome_names["Arctic Tundra"]
  if (montane_arctic_proxy == "elevation") {
    biome_class[biome_class == biome_names["Arctic Tundra"] & is_montane] <- biome_names["Montane Grassland"]
  }else if (montane_arctic_proxy == "latitude") {
    biome_class[biome_class == biome_names["Arctic Tundra"] & !is_high_latitude] <- biome_names["Montane Grassland"]
  }else {
    stop(paste0("Unknown value (",montane_arctic_proxy,") for parameter montane_arctic_proxy. Use 'elevation' or 'latitude'. Aborting."))
  }

  # other
  biome_class[is_rocks_and_ice] <- biome_names["Rocks and Ice"]
  biome_class[is_water] <- biome_names["Water"]

  return(list(biome_id = biome_class, biome_names = names(biome_names)))
}



read_pft_categories <- function(file_path) {
  # read_delim, col_types = readr::cols(), delim = ";")to suppress messages
  readr::read_delim(file_path, col_types = readr::cols(), delim = ";") %>%
    # change 1, 0.5, 0 values to TRUE and NAs (NA's can be dropped)
    dplyr::mutate_at(dplyr::vars(dplyr::starts_with(c("category_", "zone_"))),
                     function(x) ifelse(as.logical(x), TRUE, NA)) %>%
    # filter natural pfts
    dplyr::filter(category_natural) %>%
    # all binary zone columns (tropical, temperate, boreal) in one categorical
    #   zone column
    tidyr::pivot_longer(cols = starts_with("zone_"),
                 names_to = "zone",
                 names_prefix = "zone_",
                 values_to = "zone_value",
                 values_drop_na = TRUE) %>%
    # all binary category columns (natural, needle, evergreen) in one categorical # nolint
    #   category column
    tidyr::pivot_longer(cols = starts_with("category_"),
                 names_to = "category",
                 names_prefix = "category_",
                 values_to = "category_value",
                 values_drop_na = TRUE) %>%
    # delete side product - logical columns
    dplyr::select(-c("category_value", "zone_value")) %>%
    # values to lpjml_index, names to length of npft (convert to numeric)
    tidyr::pivot_longer(cols = starts_with("lpjml_index_npft_"),
                 values_to = "lpjml_index",
                 names_to = "npft_proxy",
                 names_transform = list(npft_proxy = function(x) suppressWarnings(as.numeric(x))), # nolint
                 names_prefix = "lpjml_index_npft_",
                 values_drop_na = TRUE) %>%
    return()
}
