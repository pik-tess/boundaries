#' Calculate the planetary boundary status for the land-system change boundary
#'
#' Calculate the PB status for the LSC (land-system change) boundary based
#' on a scenario LPJmL run and a reference LPJmL run.
#'
#' @param path_scenario output directory (character string) of the scenario
#'        LPJmL run where binary files (soon with metafiles) are written
#'
#' @param path_reference output directory (character string) of the reference
#'        LPJmL run where binary files (soon with metafiles) are written
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#'        as an integer vector, e.g. `1982:2011` (default)
#'
#' @param time_span_reference time span to be used for the scenario run, defined
#'        as an integer vector, e.g. `1901:1930`. Can differ in offset and
#'        length from `time_span_scenario`! If `NULL` value of
#'        `time_span_scenario` is used
#'
#' @param spatial_resolution character. Spatial resolution, available options
#'        are `"biome"` (default) and `"grid"`
#'
#' @param eurasia logical. If `spatial_resolution` = `"biome"` merge continents
#'        Europe and Asia to avoid arbitrary biome cut at europe/asia border.
#'        Defaults to `TRUE`
#'
#' @param biome_classification default to NULL, if biomes are to be calculated
#' within this function based on classify_biomes. Can be set to list with 
#' already defined biomes (xx$biome_id = numeric with ncells; xx$biome_names =
#' biome names) if biomes were already calculated with classify_biomes()
#'
#' @param avg_nyear_args list of arguments to be passed to
#'        \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#'        time series analysis
#' 
#' @param input_files optional `list` containing additional input (!) files,
#'        not in the `path_data`, e.g. if temp was not written out:
#'        `list(grid=..., temp = ..., elevation = ...)`
#' 
#' @param diff_output_files optional list for specification of output file names
#'        differing from default, which is list(grid = "grid.bin", fpc = "fpc.bin",# nolint
#'        vegc = "vegc.bin", pft_lai = "pft_lai.bin", temp = "temp.bin")
#' 
#' @param file_type replace default file type. default: ""raw
#' 
#' @param read_args list of arguments for reading input/output. only required
#'        for:
#'        nc output: specification of header_size and ncell to read in
#'                   lpjml grid input
#'        raw/clm output: specification of header_size, ncell, firstyear and
#'                   fpc_nbands (12 or 10)
#' @param ... arguments forwarded to \link[boundaries](classify_biomes)
#'
#' @examples
#' \dontrun{
#'  calc_lsc_status(path_scenario, path_reference)
#' }
#'
#' @md
#' @export
calc_lsc_status <- function(path_scenario,
                            path_reference,
                            time_span_scenario = c(1982, 2011),
                            time_span_reference = NULL,
                            spatial_resolution = "biome",
                            eurasia = TRUE,
                            biome_classification = NULL,
                            avg_nyear_args = list(),
                            input_files = list(),
                            diff_output_files = list(),
                            file_type = "raw",
                            read_args = list(
                              header_size = 0,
                              ncell = 67420,
                              firstyear = 1901,
                              fpc_nbands = 12,
                              size = 4
                            ),
                            ...
                            ) {

  # verify available temporal resolution
  spatial_resolution <- match.arg(spatial_resolution, c("biome",
                                                        "grid"))

  if (.Platform$OS.type == "windows") {
    future_plan <- future::plan("multisession")
  } else {
    future_plan <- future::plan("multicore")
  }
  on.exit(future::plan(future_plan))

  # check time_spans of scenario and reference runs
  if (is.null(time_span_reference)) {
    time_span_reference <- time_span_scenario
    nyear_ref <- NULL
  } else {
    if (diff(time_span_reference) > diff(time_span_scenario)) {
      stop(paste0("time_span_reference is longer than time_span_scenario.",
                  "Define a time_span_reference that is shorter than",
                  "time_span_scenario"))
    } else if (diff(time_span_reference) < diff(time_span_scenario)) {
      nyear_ref <- length(time_span_scenario[1]:time_span_scenario[2])
    } else {
      nyear_ref <- NULL
    }
  }

  # define file paths
  file_extension <- switch(file_type,
                           raw = ".bin",
                           meta = ".bin.json",
                           clm = ".clm",
                           nc = ".nc",
                           nc4 = ".nc4",
                           cdf = ".nc")
  # default output files with defined file_extension
  output_files <- list(grid = "grid",
                       fpc = "fpc") %>%
                  lapply(paste0, file_extension)
  #  concatenate path_data and output_files
  output_files_scenario <- lapply(output_files,
                                  function(x, path_data) {
                                  paste(path_data, x, sep = "/")
                                  },
                                  path_data = path_scenario)
  output_files_reference <- lapply(output_files,
                                  function(x, path_data) {
                                  paste(path_data, x, sep = "/")
                                  },
                                  path_data = path_reference)

  # classify biomes based on foliage projected cover (FPC) output
  if (is.null(biome_classification)) {
    biome_classes %<-% classify_biomes(
      path_data = path_reference,
      timespan = time_span_reference,
      avg_nyear_args = avg_nyear_args,
      input_files = input_files,
      diff_output_files = diff_output_files,
      file_type = file_type,
      read_args = read_args,
      ...
      )
    biome_classes <- biome_classes$biome_id
  } else {
    biome_classes <- biome_classification$biome_id
  }
  if (spatial_resolution == "biome") {
    # get continents mask - pass arg of whether to merge europe and asia
    continent_grid %<-% calc_continents_mask(path_reference, eurasia = eurasia)
  }

  # read in biome mapping
  biome_mapping <- read_biome_mapping(system.file("extdata",
                               "biomes.csv",
                               package = "boundaries"))

  # TO BE REPLACED BY lpjmlKit::read_io
  # read grid
  if (file.exists(output_files_scenario$grid) && file_type == "meta" &&
      is.null(input_files$grid)) {
      grid <- lpjmliotools::autoReadMetaOutput(
        metaFile = output_files_scenario$grid
      )
      ncell <- length(grid$lon)
      lon   <- grid$lon
      lat   <- grid$lat
  } else if (file_type %in% c("nc", "cdf", "nc4")) {
    # if nc output is defined, we need an lpjml grid to convert to 
    # the correct array size, this needs to be given in input_files$grid
    message("Reading of netcdf output is still preliminary. Please specify LPJmL grid input.") # nolint
    grid <- lpjmliotools::readGridInputBin(inFile = input_files$grid,
                             headersize = read_args$header_size,
                             ncells = read_args$ncell)
    lon <- grid$lon
    lat <- grid$lat
    ncell <- length(grid$lon)
  } else if (file.exists(output_files_scenario$grid) &&
             file_type %in% c("raw", "clm")) {
    grid <- lpjmliotools::readGridOutputBin(inFile = output_files_scenario$grid,
                                            headersize = read_args$header_size, #nolint
                                            ncells = read_args$ncell) %>%
            rename_step2month()
    lon <- grid$lon
    lat <- grid$lat
    ncell <- length(grid$lon)
  } else {
    stop(paste0("Output file ",
                output_files_scenario$grid,
                " does not exist. Make sure the specified input path_data is ",
                "correct. If your file names differ from the default, please ",
                "use diff_output_files to specify them. "))
  }
  # calculate cell area
  lpjml_grid <- rbind(lon, lat)
  cell_area <- lpjmliotools::cellarea
  # TODO replace with lpjmlkit function to calculate cellarea for the respective
  # grid
  # cell_area <-  lpjmlKit::calc_cellarea(
  # lpjmliotools::subset_array(lpjml_grid, list(coordinate = "lat"))
  # )


  # read fpc
  # TODO: convert to yearly if output is monthly
  if (file_type == "meta") {
    fpc_scenario <- lpjmliotools::autoReadMetaOutput(
      metaFile = output_files_scenario$fpc,
      getyearstart = time_span_scenario[1],
      getyearstop = time_span_reference[2]
    ) %>% rename_step2month()
    fpc_reference <- lpjmliotools::autoReadMetaOutput(
      metaFile = output_files_reference$fpc,
      getyearstart = time_span_reference[1],
      getyearstop = time_span_reference[2]
    ) %>% rename_step2month()
  } else if (file_type %in% c("nc", "nc4", "cdf")) {
    fpc_scenario <- lpjmliotools::netcdfCFT2lpjarray(
      ncInFile = output_files_scenario$fpc,
      var = "FPC",
      lon = lon,
      lat = lat
    ) %>% rename_step2month()
    fpc_reference <- lpjmliotools::netcdfCFT2lpjarray(
      ncInFile = output_files_reference$fpc,
      var = "FPC",
      lon = lon,
      lat = lat
    ) %>% rename_step2month()
  } else if (file_type %in% c("raw", "clm")) {
    fpc_scenario %<-% lpjmliotools::readCFToutput(inFile = output_files_scenario$fpc,
                                      startyear = read_args$firstyear,
                                      stopyear = time_span_scenario[2],
                                      size = read_args$size,
                                      headersize = read_args$header_size,
                                      getyearstart = time_span_scenario[1],
                                      getyearstop = time_span_scenario[2],
                                      ncells = read_args$ncell,
                                      bands = read_args$fpc_nbands) %>% rename_step2month() # nolint
    fpc_reference %<-% lpjmliotools::readCFToutput(inFile = output_files_reference$fpc,
                                      startyear = read_args$firstyear,
                                      stopyear = time_span_reference[2],
                                      size = read_args$size,
                                      headersize = read_args$header_size,
                                      getyearstart = time_span_reference[1],
                                      getyearstop = time_span_reference[2],
                                      ncells = read_args$ncell,
                                      bands = read_args$fpc_nbands) %>% rename_step2month() # nolint
  } else {
      stop(paste0("Unknown file ending (",
                  fpc_ending,
                  "). Aborting."))
  }

  fpc_nbands <- dim(fpc_scenario)[["band"]]
  npft <- fpc_nbands - 1

  pft_categories <- system.file("extdata",
                                "pft_categories.csv",
                                package = "boundaries") %>%
    read_pft_categories() %>%
    {
      if (file_type != "meta") {
        dplyr::filter(., npft_proxy == npft)
      } else {
        dplyr::filter(., is.na(npft_proxy))
      }
    }
  fpc_names <- dplyr::filter(pft_categories, category == "natural")$pft

  # TODO if lpjmlkit is used, fpc names are read in and do not have to be
  # defined here
  # if (file_type != "meta") {
  dimnames(fpc_scenario)$band <- c("natural stand fraction", fpc_names)
  dimnames(fpc_reference)$band <- c("natural stand fraction", fpc_names)

  # subset only tree pfts
  fpc_trees <- dplyr::filter(
    pft_categories,
    type == "tree" & category == "natural"
  )$pft

  tree_share_scenario <- lpjmliotools::subset_array(
    fpc_scenario,
    list(band = fpc_trees)
  )
  tree_share_reference <- lpjmliotools::subset_array(
    fpc_reference,
    list(band = fpc_trees)
  )
  # calculate actual tree cover area
  tree_cover_scenario <- (
    tree_share_scenario *
    array(lpjmliotools::subset_array(fpc_scenario,
                    list(band = c("natural stand fraction"))),
         dim = dim(tree_share_scenario)) * cell_area
  )
  tree_cover_reference <- (
    tree_share_reference *
    array(lpjmliotools::subset_array(fpc_reference,
                    list(band = c("natural stand fraction"))),
         dim = dim(tree_share_reference)) * cell_area
  )
  # sum tree pfts for forest cover
  all_tree_cover_scenario %<-% apply(tree_cover_scenario,
                                     c("cell", "year"),
                                     sum,
                                     na.rm = TRUE)
  all_tree_cover_reference %<-% apply(tree_cover_reference,
                                      c("cell", "year"),
                                      sum,
                                      na.rm = TRUE)
  # average forest over time
  avg_trees_scenario %<-% do.call(average_nyear_window,
                                  append(list(x = all_tree_cover_scenario),
                                         avg_nyear_args))
  if (!is.null(nyear_ref)) {
    avg_nyear_args["nyear_reference"] <- nyear_ref
  }

  # average forest over time
  avg_trees_reference %<-% do.call(average_nyear_window,
                                   append(list(x = all_tree_cover_reference),
                                          avg_nyear_args))

  # binary is forest biome - mask
  is_forest <- array(0,
                     dim = dim(avg_trees_reference),
                     dimnames = dimnames(avg_trees_reference))

  # binary is tropical forest biome - mask
  is_tropical_forest <- is_forest
  # binary is temperate forest biome - mask
  is_temperate_forest <- is_forest
  # binary is boreal forest biome - mask
  is_boreal_forest <- is_forest
  # init forest type mask (tropical:1, temperate:2, boreal:3)
  forest_type <- is_forest

  # 8 forest biomes
  is_forest[biome_classes %in% dplyr::filter(biome_mapping,
                                             category == "forest"
                                            )$id] <- 1

  # 2 tropical forest biomes
  is_tropical_forest[biome_classes %in% dplyr::filter(biome_mapping,
                                             category == "forest" &
                                             zone == "tropical"
                                            )$id] <- 1
  forest_type[which(is_tropical_forest == 1)] <- 1
  # 4 temperate forest biomes
  is_temperate_forest[biome_classes %in% dplyr::filter(biome_mapping,
                                             category == "forest" &
                                             zone == "temperate"
                                            )$id] <- 1
  forest_type[which(is_temperate_forest == 1)] <- 2

  # 2 boreal forest biomes
  is_boreal_forest[biome_classes %in% dplyr::filter(biome_mapping,
                                             category == "forest" &
                                             zone == "boreal"
                                            )$id] <- 1
  forest_type[which(is_boreal_forest == 1)] <- 3

  # calculate deforestation status by share of scenario and reference run
  deforestation <- avg_trees_scenario / (avg_trees_reference + 1e-9)
  deforestation[deforestation > 1] <- 1
  deforestation[deforestation > 0] <- 1 - deforestation[deforestation > 0]

  # get flexibly named time dimension
  third_dim <- names(dim(deforestation))[
    !names(dim(deforestation)) %in% c("cell")
  ]

  if (spatial_resolution == "biome") {
    # create space of combinations to loop over (even though not all make sense)
    comb <- expand.grid(
      continent = sort(unique(lpjmliotools::subset_array(continent_grid, list(coordinate = "continent")))), # nolint
      forest = sort(unique(factor(forest_type)))[-1]
    )
    for (idx in seq_len(nrow(comb))) {
      sub_deforest <- deforestation
      # replace for every combination a subset of cells with mean of each
      sub_cells <- {
        # match forest type for cells
        forest_type == comb$forest[idx] &
        # match continent for cells
        array(
          lpjmliotools::subset_array(continent_grid, list(coordinate = "continent"), drop = FALSE) == comb$continent[idx], # nolint
          dim = dim(deforestation),
          dimnames = dimnames(deforestation)
        )
      }
      if (!any(sub_cells)) {
        next
      }
      sub_deforest[!sub_cells] <- NA
      sub_deforest <- apply(
        sub_deforest,
        third_dim,
        function(x) {
          # mean over cell subset of forest type and continent
          # weighted by cell area
          return(rep(weighted.mean(x, cell_area, na.rm = TRUE), length(x)))
        }
      )
      deforestation[sub_cells] <- sub_deforest[sub_cells]
    }
  }
  # init pb_status array with 0 = no data and initial dimensions
  pb_status <- array(0,
                     dim = dim(deforestation),
                     dimnames = dimnames(deforestation))
  # boundaries for tropical forest and boreal forest after Steffen et al. 2015
  #   (https://doi.org/10.1126/science.1259855)
  #   share for safe space <= 0.15, uncertainty zone < 0.4 & >0.15 and
  #   for high risk >= 0.4
  # high risk
  pb_status[
    (is_tropical_forest |
    is_boreal_forest) &
    deforestation >= 0.4
  ] <- 3
  # uncertainty zone
  pb_status[
    (is_tropical_forest |
    is_boreal_forest) &
    deforestation < 0.4 &
    deforestation > 0.15
  ] <- 2
  # safe space
  pb_status[
    (is_tropical_forest |
    is_boreal_forest) &
    deforestation <= 0.15
  ] <- 1

  # boundaries for temperate forest after Steffen et al. 2015
  #   (https://doi.org/10.1126/science.1259855)
  #   share for safe space <= 0.5, uncertainty zone < 0.7 & >0.5 and
  #   for high risk >= 0.7
  # high risk
  pb_status[
    is_temperate_forest &
    deforestation >= 0.7
  ] <- 3
  # uncertainty zone
  pb_status[
    is_temperate_forest &
    deforestation < 0.7 &
    deforestation > 0.5
  ] <- 2
  # safe space
  pb_status[
    is_temperate_forest &
    deforestation <= 0.5
  ] <- 1

  return(pb_status)
}

read_biome_mapping <- function(file_path) {
  # read_delim, col_types = readr::cols(), delim = ";")to suppress messages
  readr::read_delim(file_path, col_types = readr::cols(), delim = ";") %>%
    # change 1, 0.5, 0 values to TRUE and NAs (NA's can be dropped)
    dplyr::mutate_at(dplyr::vars(dplyr::starts_with(c("category_", "zone_"))),
                     function(x) ifelse(as.logical(x), TRUE, NA)) %>%
    # all binary zone columns (tropical, temperate, boreal) in one categorical
    #   zone column
    tidyr::pivot_longer(cols = starts_with("zone_"),
                 names_to = "zone",
                 names_prefix = "zone_",
                 values_to = "zone_value",
                 values_drop_na = TRUE) %>%
    # all binary category columns (forest, xx) in one categorical
    #   category column
    tidyr::pivot_longer(cols = starts_with("category_"),
                 names_to = "category",
                 names_prefix = "category_",
                 values_to = "category_value",
                 values_drop_na = TRUE) %>%
    # delete side product - logical columns
    dplyr::select(-c("category_value", "zone_value")) %>%
    return()
}