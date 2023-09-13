# test calc_biosphere_status
library(raster)
devtools::load_all("/p/projects/open/Fabian/LPJbox/boundaries_development")
devtools::load_all("/p/projects/open/Fabian/LPJbox/biospheremetrics_review_paper/")
inFol_lu <- "/p/projects/open/Fabian/runs/metrics_202308/output/lu_1500_2016/"
inFol_pnv <- "/p/projects/open/Fabian/runs/metrics_202308/output/pnv_1500_2016/"
varnames <- data.frame(row.names = c("grid",            "npp",         "pft_npp",         "pft_harvest",             "pft_rharvest",             "firec",         "timber_harvest",          "cftfrac",         "fpc",        "pft_lai",          "vegc"),
                          outname = c("grid.bin.json",  "npp.bin.json","pft_npp.bin.json","pft_harvest.pft.bin.json","pft_rharvest.pft.bin.json","firec.bin.json","timber_harvestc.bin.json","cftfrac.bin.json","fpc.bin.json","pft_lai.bin.json","vegc.bin.json"),
                          timestep = c("Y",             "Y",           "Y",               "Y",                       "Y",                        "Y",             "Y",                       "Y",               "Y" ,          "Y",               "Y"))
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
files_baseline <- list(grid = paste0(inFol_pnv,varnames["grid","outname"]),
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
                       vegc =""
)
files_reference <- list(npp = paste0(inFol_pnv,varnames["npp","outname"]))

pb_bi <- calc_biosphere_status(files_scenario = files_scenario,
                               files_reference = files_reference,
                               files_baseline = files_baseline,
                               time_span_scenario = as.character(2000:2016),
                               time_span_baseline = as.character(2000:2016),
                               time_span_reference = as.character(1510:1520),
                               spatial_resolution = "subglobal",
                               thresholds = NULL,
                               gridbased = T,
                               npp_threshold = 20
)
boundaries::plot_status(status_data = list(biosphere=pb_bi[,17]))
str(pb_bi)
