library(tibble)
library(lpjmlkit)
library(dplyr)

# paths
model_path <- "." # "/p/projects/open/Jannes/lpjml/default/LPJmL_internal/"
output_path <- "." # "/p/projects/open/Jannes/testing/boundaries/runs/"

#==================INPUT RUNS===================================================

params <- tibble(
  sim_name = c(
    "spinup",
    "lu_1500_2016",
    "pnv_1500_2016",
    "all_crops"
  ),
  sowing_date_option = "fixed_sdate",
  sdate_fixyear = 2000,
  others_to_crop = TRUE,
  tillage_type = "all",
  intercrop = FALSE,
  grazing = "livestock",
  prescribe_lsuha = TRUE,
  crop_resp_fix = TRUE,
  radiation = "radiation_lwdown",
  landuse = c("no", "yes","no","all_crops"),
  reservoir = c(FALSE, TRUE, FALSE, TRUE),
  wateruse = c("no", "all","no", "all"),
  irrigation = c(rep("lim",3), "pot"),
  with_nitrogen = c(rep("lim",3), "unlim"),
  input.prec.name  =  "/p/projects/lpjml/input/historical/GSWP3-W5E5/pr_gswp3-w5e5_1901-2016.clm",
  input.humid.name  =  "/p/projects/lpjml/input/historical/GSWP3-W5E5/huss_gswp3-w5e5_1901-2016.clm",
  input.temp.name  =  "/p/projects/lpjml/input/historical/GSWP3-W5E5/tas_gswp3-w5e5_1901-2016.clm",
  input.lwdown.name  =  "/p/projects/lpjml/input/historical/GSWP3-W5E5/rlds_gswp3-w5e5_1901-2016.clm",
  input.swdown.name  =  "/p/projects/lpjml/input/historical/GSWP3-W5E5/rsds_gswp3-w5e5_1901-2016.clm",
  input.wind.name  =  "/p/projects/lpjml/input/historical/GSWP3-W5E5/wind_gswp3-w5e5_1901-2016.clm",
  input.co2.name  =  "/p/projects/lpjml/input/historical/input_VERSION2/co2_1841-2018.dat",
  input.landuse.name  =  "/p/projects/lpjml/input/historical/input_toolbox_30arcmin/paper_version/cft_default_cft_aggregation_30min_1500-2017_fallow_to_others_64bands.clm",
  input.fertilizer_nr.name = "/p/projects/lpjml/input/historical/input_toolbox_30arcmin/paper_version/fert_N_default_cft_aggregation_30min_1860-2017.clm",
  input.manure_nr.name = "/p/projects/lpjml/input/historical/input_toolbox_30arcmin/paper_version/manure_N_default_cft_aggregation_30min_1860-2017.clm",
  store_climate = FALSE,
  output_metafile = TRUE,
  dependency = c(NA, rep("spinup",3)),
  wtime = c("10:00", rep("4:00",3)),
  outputyear = c(NA, rep(1500,3)),
  firstyear = 1901, 
  lastyear = c(1901, rep(2016,3)),
  restart_year = c(1900, rep(2016,3)),
  write_restart = TRUE,
  river_routing = FALSE,
  nspinup = c(10000, rep(401,3)),
  grid_scaled = TRUE,
  double_harvest = FALSE,
  startgrid = 27410,
  endgrid = 27411
)

# write config files
config_details <- write_config(
  x = params,
  model_path = model_path,
  js_filename = "lpjml_newoutputs.cjson",
  sim_path = output_path,
  debug = TRUE
)

# lpjcheck
check_config(config_details, model_path, output_path)

# multi lpjsubmit
job_details <- submit_lpjml(x = config_details,
                            model_path,
                            #constraint = "haswell",
                            sclass = "standby",
                            output_path,
                            group = "worldtrans",
                            ntasks = 1
                            )

run_lpjml(x = config_details,
          model_path = model_path,
          output_path = output_path)
