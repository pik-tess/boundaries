# test calc_biosphere_status
devtools::load_all("/p/projects/open/Fabian/LPJbox/boundaries_development")
inFol_lu <- "/p/projects/open/Fabian/runs/metrics_202306/output/lu_1500_2014/"
inFol_pnv <- "/p/projects/open/Fabian/runs/metrics_202306/output/pnv_1500_2014/"
varnames <- data.frame(row.names = c("grid",         "npp",     "pft_npp",    "pft_harvest",        "pft_rharvest",        "firec",    "timber_harvest",     "cftfrac",    "fpc"),
                          outname = c("grid.bin.json",  "mnpp.bin.json","pft_npp.bin.json","pft_harvest.pft.bin.json","pft_rharvest.pft.bin.json","firec.bin.json","timber_harvestc.bin.json","cftfrac.bin.json","fpc.bin.json"),
                          timestep = c("Y",             "M",       "Y",          "Y",                  "Y",                   "Y",        "Y",                  "Y",          "Y"))
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
pb_bi <- calc_biosphere_status(files_scenario = files_scenario,
                               files_reference = files_reference,
                               time_span_scenario = as.character(2000:2014),
                               time_span_reference = as.character(1510:1520),
                               reference_npp_file = files_reference$npp,
                               spatial_resolution = "grid",
                               thresholds = NULL)
str(pb_bi)
