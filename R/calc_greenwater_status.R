#' Calculate the planetary boundary status for the greenwater boundary
#'
#' Calculate the PB status for the greenwater (former freshwater) boundary based
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
#' method is `c("wang-erlandsson2022")` based on
#' [Wang-Erlandsson et al. 2022](https://doi.org/10.1038/s43017-022-00287-8).
#'
#' @param temporal_resolution character. Temporal resolution, available options
#' are `"annual"` (default) and `"monthly"`. `"annual"` describes the mean only
#' over months that are transgressed.
#'
#' @param prefix_monthly_output character. Provide a prefix if required for
#' monthly LPJmL output files, e.g. `"m"` for `"mrootmoist.bin"` instead of
#' `"rootmoist.bin"` (default is `""`)
#'
#' @param avg_nyear_args list of arguments to be passed to
#' \link[pbfunctions]{average_nyear_window} (see for more info). To be used for
#' time series analysis
#'
#' @param startyear first year of simulation
#'
#'
#' @examples
#' \dontrun{
#'  calc_greenwater_status(path_scenario, path_reference)
#' }
#'
#' @md
#' @export
calc_greenwater_status <- function(path_scenario,
                                   path_reference,
                                   time_span_scenario = c(1982, 2011),
                                   time_span_reference = NULL,
                                   method = "wang-erlandsson2022",
                                   #temporal_resolution = "annual",
                                   # Q < 1m³/s
                                   prefix_monthly_output = "m",
                                   avg_nyear_args = list(),
                                   # to be replaced by lpjmlKit::read_output
                                   start_year = 1850) { #OBS: here start_year for scenario!
  # verify available methods
  method <- match.arg(method, c("wang-erlandsson2022"))
  # verify available temporal resolution
  #temporal_resolution <- match.arg(temporal_resolution, c("annual"))
  
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
  
  # reference rootmoisture
  # TO BE REPLACED BY lpjmlKit::read_output ---------------------------------- #
  #   hardcoded values to be internally replaced
  # file_reference %<-% tmp_read_monthly(
  #   file_name = paste0(path_reference,
  #                      "/",
  #                      prefix_monthly_output,
  #                      "rootmoist.bin"),
  #   time_span = time_span_reference,
  #   start_year = start_year,
  #   nstep = 12,
  #   ncell = 67420,
  #   nbands = 1,
  #   size = 4
  # )
  # 
  file_reference <- brick(paste0(path_reference,
                                      "/",
                                      prefix_monthly_output,
                                      "rootmoist.nc"))
  
  # -------------------------------------------------------------------------- #
  
  # scenario rootmoisture
  # TO BE REPLACED BY lpjmlKit::read_output ---------------------------------- #
  #   hardcoded values to be internally replaced
  # file_scenario %<-% tmp_read_monthly(
  #   file_name = paste0(path_scenario,
  #                      "/",
  #                      prefix_monthly_output,
  #                      "rootmoist.bin"),
  #   time_span = time_span_scenario,
  #   start_year = start_year,
  #   nstep = 12,
  #   ncell = 67420,
  #   nbands = 1,
  #   size = 4
  # )
  # 
  #OBS: I use picontrol (1600-1850) as the reference period
  file_scenario <- brick(paste0(path_scenario,
                                     "/",
                                     prefix_monthly_output,
                                     "rootmoist.nc"))
  
  file_scenario <- file_scenario[[((time_span_scenario[1]-start_year)*12+1):((time_span_scenario[2]+1-start_year)*12)]]
  # -------------------------------------------------------------------------- #
  
  #calculate the green water 5% and 95% quantiles of the baseline period
  quants <- calc_water_baseline(file_reference)
  
  # calculate GW dry & wet departures and return df with percentage 
  # of annual area of departure
  library(lpjmliotools) #to do: replace with lpjmlKit
  arearaster         <- raster(ncols=720,nrows=360)
  arearaster[cellFromXY(arearaster,cbind(lon,lat))] <- lpjmliotools::cellarea
  
  departure.df <- calc_water_depart(file_scenario,arearaster,area_sum,quants)
  
  pb_status <- departure.df
  return(pb_status)
}

# raster quantile functions
q5  <- function(x) {quantile(x, probs = 0.05, na.rm = T)}
q95 <- function(x) {quantile(x, probs = 0.95, na.rm = T)}  

# calculate the baseline quantiles
calc_water_baseline <- function(file_reference) {
  
  #nyr of reference baseline
  nyr_base <- nlayers(file_reference)/12
  
  #empty yr lists
  dry_base_yr <- list()
  wet_base_yr <- list()
  
  for (j in 1:nyr_base) {
    print(j)
    index_j  <- seq(12*j-11,12*j)
    file_j   <- file_reference[[index_j]]
    
    dry_j    <- min(file_j) #driest month per gridcell for year i
    dry_base_yr[j] <- dry_j
    
    wet_j    <- max(file_j) #wettest month per gridcell for year i
    wet_base_yr[j] <- wet_j
    if (j == nyr_base){
      #stack
      dry_base_yr <- stack(dry_base_yr)
      wet_base_yr <- stack(wet_base_yr)
    }
  }
  # calc 5% and 95% percentile for each cell 
  # for each year over all months over baseline period
  
  q5_base  <- calc(x = dry_base_yr, fun = q5)
  q5_base  <- brick(q5_base) #otherwise export is buggy
  
  q95_base <- calc(x = wet_base_yr, fun = q95)
  q95_base <- brick(q95_base)
  
  quants <- list(q5_base,q95_base)
  return(quants)
}

# calculate GW dry & wet departures and return df with annual area of departure
calc_water_depart <- function(file_scenario,arearaster,area_sum,quants) {
  
  q5_base <- brick(quants[1])
  q95_base <- brick(quants[2])
  
  #nyr variables
  nyr_scn <- nlayers(file_scenario)/12
  
  #create departure data.frame
  departure.df     <- matrix(nrow = nyr_scn,ncol = 3) 
  departure.df[,1] <- time_span_scenario[1] : (nyr_scn+time_span_scenario[1] -1) #create year col
  colnames(departure.df) <- c("year","dry","wet")
  
  # for each hist year calc area of rootmoist wet & dry departures
  for (i in 1:nyr_scn) {
    print(i)
    index_i  <- seq(12*i-11,12*i)
    file_i   <- file_scenario[[index_i]]
    
    
    #identify cells with dry departures
    dry_i       <- min(file_i) #driest month per gridcell for year i
    dry_i[dry_i >= q5_base] <- NA
    dry_i[dry_i >= 0] <- 1
    dry_i       <- dry_i*arearaster
    dry_i       <- cellStats(dry_i, stat='sum')
    dry_perc_i  <- dry_i/area_sum
    dry_perc_i  <- dry_perc_i*100
    departure.df[i,2] <- dry_perc_i #assign value to df
    
    
    #identify cells with wet departures
    wet_i       <- max(file_i) #wettest month per gridcell for year i -> ignores which month
    wet_i[wet_i <= q95_base] <- NA
    wet_i[wet_i >= 0] <- 1
    wet_i       <- wet_i*arearaster
    wet_i       <- cellStats(wet_i, stat='sum')
    wet_perc_i  <- wet_i/area_sum
    wet_perc_i  <- wet_perc_i*100
    departure.df[i,3] <- wet_perc_i #assign value to df
    
    if (i == nyr_scn) {
      departure.df <- data.frame(departure.df)
      return(departure.df)
    }
  }
}

