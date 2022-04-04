#' Calculate the planetary boundary status for the land-system change boundary
#'
#' Calculate the PB status for the LSC (land-system change) boundary based
#' on a scenario LPJmL run and a reference LPJmL run.
#'
#' @param path_scenario output directory (character string) of the scenario
#' LPJmL run where binary files (soon with metafiles) are written
#'
#' @param path_reference output directory (character string) of the reference
#' LPJmL run where binary files (soon with metafiles) are written
#'
#' @param time_span_scenario time span to be used for the scenario run, defined
#' as an integer vector, e.g. `1982:2011` (default)
#'
#' @param time_span_reference time span to be used for the scenario run, defined
#' as an integer vector, e.g. `1901:1930`. Can differ in offset and length from
#' `time_span_scenario`! If `NULL` value of `time_span_scenario` is used
#'
#' @param method method (character string) to be used , currently available
#' method is `c("gerten2020")` based on
#' [Gerten et al. 2020](https://doi.org/10.1038/s41893-019-0465-1).
#'
#' @param temporal_resolution character. Temporal resolution, available options
#' are `"annual"` (default) and `"monthly"`
#'
#' @param cut_min double. Exclude boundary calculations for Q < cut_min
#'
#' @param prefix_monthly_output character. Provide a prefix if required for
#' monthly LPJmL output files, e.g. `"m"` for `"mdischarge.bin"` instead of
#' `"discharge.bin"` (default is `""`)
#'
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#' time series analysis
#'
#' @examples
#' \dontrun{
#'  calc_bluewater_status(path_scenario, path_reference)
#' }
#'
#' @md
#' @export
calc_lsc_status <- function(path_luinput,
                            path_scenario,
                            path_reference,
                            time_span_scenario = c(1982, 2011),
                            time_span_reference = NULL,
                            method = "gerten2020",
                            temporal_resolution = "annual",
                            # ???
                            cut_min = 0.0864,
                            vegc_proxy = TRUE,
                            # prefix_monthly_output = "",
                            avg_nyear_args = list(),
                            # to be replaced by lpjmlKit::read_output
                            start_year = 1901) {

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

  # TO BE REPLACED BY lpjmlKit::read_output -------------------------------- #
  # read grid
  ncell <- 67420
  size <- 2
  grid_file <- file(paste(path_reference, "grid.bin", sep = "/"), "rb")
  lpjml_grid <- readBin(grid_file, integer(), n = 2 * ncell, size = size) /
                100
  close(grid_file)
  dim(lpjml_grid) <- c(coordinate = 2, cell = ncell)
  dimnames(lpjml_grid) <- list(coordinate = c("lon", "lat"),
                               cell = seq_len(ncell))
  # ------------------------------------------------------------------------ #
  cell_area <-  lpjmlKit::calc_cellarea(
    lpjmlKit::subset_array(lpjml_grid, list(coordinate = "lat"))
  )
  # TO BE REPLACED BY lpjmlKit::read_input ----------------------------------- #
  #   hardcoded values to be internally replaced
  input_data_size <- 4
  header <- suppressWarnings(lpjmlKit::read_header(filename = path_luinput))
  headersize <- lpjmlKit::get_headersize(header)
  input_file <- file(path_luinput, "rb")
  seek(input_file,
       where = headersize +
               (time_span_scenario[1] -
               lpjmlKit::get_header_item(header, "firstyear")) *
               lpjmlKit::get_header_item(header, "nbands") *
               lpjmlKit::get_header_item(header, "ncell") *
               input_data_size,
       origin = "start")
  landuse <- readBin(input_file,
                     double(),
                     n = (time_span_scenario[2] - time_span_scenario[1] + 1) *
                         lpjmlKit::get_header_item(header, "ncell") *
                         lpjmlKit::get_header_item(header, "nbands"),
                     size = input_data_size) *
              lpjmlKit::get_header_item(header, "scalar")
  close(input_file)      #remove to save space
  dim(landuse) <- c(band = unname(lpjmlKit::get_header_item(header, "nbands")),
                    cell = unname(lpjmlKit::get_header_item(header, "ncell")),
                    year = length(
                        time_span_scenario[1] : time_span_scenario[2]
                    ))
  dimnames(landuse) <- list(
    band = unname(seq_len(lpjmlKit::get_header_item(header, "nbands"))),
    cell = unname(seq_len(lpjmlKit::get_header_item(header, "ncell"))),
    year = time_span_scenario[1] : time_span_scenario[2])

  # TO BE REPLACED BY lpjmlKit::read_output ---------------------------------- #
  #   hardcoded values to be internally replaced
  # read foliage protected cover (fpc) output
  ncell <- 67420
  fpc_names <- c("natvegfrac", # 1
                 "Tropical Broadleaved Evergreen Tree", # 2
                 "Tropical Broadleaved Raingreen Tree", # 3
                 "Temperate Needleleaved Evergreen Tree", # 4
                 "Temperate Broadleaved Evergreen Tree", # 5
                 "Temperate Broadleaved Summergreen Tree", # 6
                 "Boreal Needleleaved Evergreen Tree", # 7
                 "Boreal Broadleaved Summergreen Tree", # 8
                 "Boreal Needleleaved Summergreen Tree", # 9
                 "Tropical C4 grass", # 10
                 "Temperate C3 grass", # 11
                 "Polar C3 grass") # 12

  fpc_scenario %<-% tmp_read_pft_yearly(
    file_name = paste0(path_scenario,
                       "/",
                       "fpc.bin"),
    time_span = time_span_scenario,
    start_year = start_year,
    ncell = ncell,
    nbands = 12,
    band_names = fpc_names,
    size = 4
  )

  fpc_reference %<-% tmp_read_pft_yearly(
    file_name = paste0(path_reference,
                       "/",
                       "fpc.bin"),
    time_span = time_span_reference,
    start_year = start_year,
    ncell = ncell,
    nbands = 12,
    band_names = fpc_names,
    size = 4
  )
  # subset only tree pfts
  tree_share_scenario <- lpjmlKit::subset_array(
    fpc_scenario,
    list(band = c("Tropical Broadleaved Evergreen Tree", # 2
                  "Tropical Broadleaved Raingreen Tree", # 3
                  "Temperate Needleleaved Evergreen Tree", # 4
                  "Temperate Broadleaved Evergreen Tree", # 5
                  "Temperate Broadleaved Summergreen Tree", # 6
                  "Boreal Needleleaved Evergreen Tree", # 7
                  "Boreal Broadleaved Summergreen Tree", # 8
                  "Boreal Needleleaved Summergreen Tree"))
  )
  tree_cover_scenario <- (
    tree_share_scenario *
    array(lpjmlKit::subset_array(fpc_scenario,
                    list(band = c("natvegfrac"))),
         dim = dim(tree_cover_scenario)) * cell_area
  )
  tree_share_reference <- lpjmlKit::subset_array(
    fpc_reference,
    list(band = c("Tropical Broadleaved Evergreen Tree", # 2
                  "Tropical Broadleaved Raingreen Tree", # 3
                  "Temperate Needleleaved Evergreen Tree", # 4
                  "Temperate Broadleaved Evergreen Tree", # 5
                  "Temperate Broadleaved Summergreen Tree", # 6
                  "Boreal Needleleaved Evergreen Tree", # 7
                  "Boreal Broadleaved Summergreen Tree", # 8
                  "Boreal Needleleaved Summergreen Tree"))
  )
  tree_cover_reference <- (
    tree_share_reference *
    array(lpjmlKit::subset_array(fpc_reference,
                    list(band = c("natvegfrac"))),
         dim = dim(tree_cover_scenario)) * cell_area
  )
  # sum tree pfts for forest cover
  forest_cover_scenario %<-% apply(tree_cover_scenario,
                                   c("cell", "year"),
                                   sum,
                                   na.rm = TRUE) # * cell_area
  forest_cover_reference %<-% apply(tree_cover_reference,
                                    c("cell", "year"),
                                    sum,
                                    na.rm = TRUE)
  # average fpc
  avg_forest_scenario %<-% do.call(average_nyear_window,
                                append(list(x = forest_cover_scenario),
                                       avg_nyear_args))
  if (!is.null(nyear_ref)) {
    avg_nyear_args["nyear_reference"] <- nyear_ref
  }

  avg_forest_reference %<-% do.call(average_nyear_window,
                                append(list(x = forest_cover_reference),
                                       avg_nyear_args))

  # classify biomes based on foliage protected cover (FPC) output
  biome_classes <- classify_biomes(
      path_data = path_reference,
      time_span = time_span_reference,
      vegc_proxy = vegc_proxy,
      avg_nyear_args = avg_nyear_args,
      # to be replaced by lpjmlKit::read_output
      start_year = start_year)

  # forest biomes masks
  is_forest <- array(0,
                     dim = dim(biome_classes),
                     dimnames = dimnames(biome_classes))

  is_tropical_forest <- is_forest
  is_temperate_forest <- is_forest
  is_boreal_forest <- is_forest

  is_forest[biome_classes %in% seq_len(8)] <- 1

  # tropical forest mask
  is_tropical_forest[biome_classes %in% seq_len(2)] <- 1
  # temperate
  is_temperate_forest[biome_classes %in% seq(3, 6)] <- 1

  # boreal
  is_boreal_forest[biome_classes %in% seq(7, 8)] <- 1

  deforestation <- avg_forest_scenario / (avg_forest_reference + 1e-9)
  deforestation[deforestation > 1] <- 1
  deforestation[deforestation > 0] <- 1 - deforestation[deforestation > 0]

  third_dim <- names(dim(deforestation))[
    !names(dim(deforestation)) %in% c("cell")
  ]

  if (!spatial_resolution == "biome") {
    deforestation <- lpjmlKit::replace_array(
      deforestation,
      list(cell = which(is_tropical_forest == 1)),
      apply(
        lpjmlKit::subset_array(deforestation,
                               list(cell = which(is_tropical_forest == 1)),
                               drop = FALSE),
        third_dim,
        function(x) {
          return(rep(mean(x, na.rm = TRUE), length(x)))
        }
      )
    ) %>%
      lpjmlKit::replace_array(
        list(cell = which(is_boreal_forest == 1)),
        apply(
          lpjmlKit::subset_array(deforestation,
                                 list(cell = which(is_boreal_forest == 1)),
                                 drop = FALSE),
          third_dim,
          function(x) {
            return(rep(mean(x, na.rm = TRUE), length(x)))
          }
        )
      ) %>%
      lpjmlKit::replace_array(
        list(cell = which(is_temperate_forest == 1)),
        apply(
          lpjmlKit::subset_array(deforestation,
                                 list(cell = which(is_temperate_forest == 1)),
                                 drop = FALSE),
          third_dim,
          function(x) {
            return(rep(mean(x, na.rm = TRUE), length(x)))
          }
        )
      )
  }

  pb_status <- array(0,
                     dim = dim(deforestation),
                     dimnames = dimnames(deforestation))
  pb_status[
    is_tropical_forest |
    is_boreal_forest &
    deforestation >= 0.4
  ] <- 3
  pb_status[
    is_tropical_forest |
    is_boreal_forest &
    deforestation < 0.4 &
    deforestation > 0.15
  ] <- 2
  pb_status[
    is_tropical_forest |
    is_boreal_forest &
    deforestation <= 0.15
  ] <- 1

  pb_status[
    is_temperate_forest &
    deforestation >= 0.7
  ] <- 3
  pb_status[
    is_temperate_forest &
    deforestation < 0.7 &
    deforestation > 0.5
  ] <- 2
  pb_status[
    is_temperate_forest &
    deforestation <= 0.5
  ] <- 1

}