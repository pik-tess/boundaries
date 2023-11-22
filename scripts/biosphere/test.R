# test calc_biosphere_status
library(raster)
library(lpjmliotools)
devtools::load_all("/p/projects/open/Fabian/LPJbox/boundaries_development")
devtools::load_all("/p/projects/open/Fabian/LPJbox/biosphere_indicators_github")
inFol_lu <- "/p/projects/open/Fabian/runs/metrics_202308/output/lu_1500_2016/"
inFol_pnv <- "/p/projects/open/Fabian/runs/metrics_202308/output/pnv_1500_2016/"
varnames <- data.frame(row.names = c("grid",            "terr_area",         "npp",         "pft_npp",         "pft_harvest",             "pft_rharvest",             "firec",         "timber_harvest",          "cftfrac",         "fpc",        "pft_lai",          "vegc"),
                          outname = c("grid.bin.json",  "terr_area.bin.json","npp.bin.json","pft_npp.bin.json","pft_harvest.pft.bin.json","pft_rharvest.pft.bin.json","firec.bin.json","timber_harvestc.bin.json","cftfrac.bin.json","fpc.bin.json","pft_lai.bin.json","vegc.bin.json"),
                          timestep = c("Y",             "Y",                 "Y",           "Y",               "Y",                       "Y",                        "Y",             "Y",                       "Y",               "Y" ,          "Y",               "Y"))
files_scenario <- list(grid = paste0(inFol_lu,varnames["grid","outname"]),
                       terr_area = paste0(inFol_lu,varnames["terr_area","outname"]),
                       npp = paste0(inFol_lu,varnames["npp","outname"]),
                       pft_npp = paste0(inFol_lu,varnames["pft_npp","outname"]),
                       pft_harvestc = paste0(inFol_lu,varnames["pft_harvest","outname"]),
                       pft_rharvestc = paste0(inFol_lu,varnames["pft_rharvest","outname"]),
                       firec = paste0(inFol_lu,varnames["firec","outname"]),
                       timber_harvestc = paste0(inFol_lu,varnames["timber_harvest","outname"]),
                       cftfrac = paste0(inFol_lu,varnames["cftfrac","outname"]),
                       fpc = paste0(inFol_lu,varnames["fpc","outname"])
)
files_baseline <- list(grid = paste0(inFol_pnv,varnames["grid","outname"]),
                       terr_area = paste0(inFol_lu,varnames["terr_area","outname"]),
                       npp = paste0(inFol_pnv,varnames["npp","outname"]),
                       pft_npp = paste0(inFol_pnv,varnames["pft_npp","outname"]),
                       pft_harvestc = paste0(inFol_pnv,varnames["pft_harvest","outname"]),
                       pft_rharvestc = paste0(inFol_pnv,varnames["pft_rharvest","outname"]),
                       firec = paste0(inFol_pnv,varnames["firec","outname"]),
                       timber_harvestc = paste0(inFol_pnv,varnames["timber_harvest","outname"]),
                       cftfrac = paste0(inFol_pnv,varnames["cftfrac","outname"]),
                       fpc = paste0(inFol_pnv,varnames["fpc","outname"]),
                       temp = "/p/projects/lpjml/input/historical/GSWP3-W5E5/tas_gswp3-w5e5_1901-2016.clm",
                       pft_lai = "",
                       vegc = ""
)
files_reference <- list(npp = paste0(inFol_pnv,varnames["npp","outname"]))

pb_bi_grid <- calc_status(boundary = c("biosphere"),
                      path_scenario = inFol_lu,
                      path_reference = inFol_pnv,
                      path_baseline = inFol_pnv,
                      #files_baseline = list(npp = "/p/projects/open/Fabian/runs/metrics_202308/output/pnv_1500_2016/npp.bin.json"),
                      time_span_scenario = as.character(2000:2016),
                      time_span_baseline = as.character(2000:2016),
                      time_span_reference = as.character(1510:1520),
                      input_files = list(prec = "/p/projects/lpjml/input/historical/CRUDATA_TS3_23/gpcc_v7_cruts3_23_precip_1901_2013.clm",
                                         temp = "/p/projects/lpjml/input/historical/CRUDATA_TS3_23/cru_ts3.23.1901.2014.tmp.dat.clm"),
                      spatial_resolution = "grid",
                      savanna_proxy = list(vegc = 7500))

boundaries::plot_status(pb_bi_grid)

devtools::load_all("/p/projects/open/Fabian/LPJbox/boundaries_development")

pb_bi_subglobal <- calc_status(boundary = c("biosphere"),
                          path_scenario = inFol_lu,
                          path_reference = inFol_pnv,
                          path_baseline = inFol_pnv,
                          #files_baseline = list(npp = "/p/projects/open/Fabian/runs/metrics_202308/output/pnv_1500_2016/npp.bin.json"),
                          time_span_scenario = as.character(2009:2013),
                          time_span_baseline = as.character(2009:2013),
                          time_span_reference = as.character(1511:1515),
                          input_files = list(prec = "/p/projects/lpjml/input/historical/CRUDATA_TS3_23/gpcc_v7_cruts3_23_precip_1901_2013.clm",
                                             temp = "/p/projects/lpjml/input/historical/CRUDATA_TS3_23/cru_ts3.23.1901.2014.tmp.dat.clm",
                                             elevation = "/p/projects/lpjml/input/historical/input_VERSION2/elevation.bin"),
                          spatial_resolution = "subglobal",
                          savanna_proxy = list(vegc = 7500))

boundaries::plot_status(pb_bi_subglobal)

pb_bi_global <- calc_status(boundary = c("biosphere"),
                               path_scenario = inFol_lu,
                               path_reference = inFol_pnv,
                               path_baseline = inFol_pnv,
                               #files_baseline = list(npp = "/p/projects/open/Fabian/runs/metrics_202308/output/pnv_1500_2016/npp.bin.json"),
                               time_span_scenario = as.character(2000:2016),
                               time_span_baseline = as.character(2000:2016),
                               time_span_reference = as.character(1510:1520),
                               input_files = list(prec = "/p/projects/lpjml/input/historical/CRUDATA_TS3_23/gpcc_v7_cruts3_23_precip_1901_2013.clm",
                                                  temp = "/p/projects/lpjml/input/historical/CRUDATA_TS3_23/cru_ts3.23.1901.2014.tmp.dat.clm"),
                               spatial_resolution = "global")

boundaries::plot_status_global(pb_bi_global)

lpjmliotools::plotGlobalManToScreen(data = pb_bi_subglobal[,1], title = "pb biosphere subglobal",
                                    brks = c(0,0.1,0.2,1), palette = c("green","yellow","red"),
                                    legendtitle = "",legYes = T)
lpjmliotools::plotGlobalManToScreen(data = pb_bi_grid[,1], title = "pb biosphere grid",
                                    brks = c(0,0.1,0.2,1), palette = c("green","yellow","red"),
                                    legendtitle = "",legYes = T)
