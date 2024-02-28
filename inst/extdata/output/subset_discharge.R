library(lpjmlkit)
library(tibble)

output_path <- "." # "/p/projects/open/Jannes/repos/boundaries/inst/extdata/output" # nolint

base_path_lu <- "." # "/p/projects/open/Johanna/boundaries/lpjml/output/lu_1500_2017_mg" # nolint
base_path_pnv <- "." # "/p/projects/open/Johanna/boundaries/lpjml/output/pnv_1500_2017_mg/" # nolint


input_files <- c(
  # discharge lu
  file.path(base_path_lu, "mdischarge.bin.json"),
  # discharge pnv
  file.path(base_path_pnv, "mdischarge.bin.json")
)

scenarios <- c("lu_1500_2016", "pnv_1500_2016")

for (input_file in input_files) {

  scen <- scenarios[1]
  scenarios <- scenarios[-1]
  # output name for new input
  output_meta_name <- strsplit(basename(input_file), "^m")[[1]][2]
  output_name <- strsplit(output_meta_name, ".json")[[1]][1]

  output_meta <- read_meta(input_file) |>
    as_list()

  output_meta$filename <- output_name
  output_meta$firstcell <- 27410
  output_meta$ncell <- 2

  jsonlite::write_json(
    output_meta,
    file.path(output_path, scen, basename(output_meta_name)),
    simplifyVector = TRUE,
    auto_unbox = TRUE,
    pretty = TRUE
  )


  input <- read_io(
    input_file,
    subset = list(cell = as.character(c(27410, 27411))),
  ) |>
    as_array() |>
    drop() |>
    as.vector()

  # write new input data to file with new header
  output <- file.path(output_path, scen, output_name)
  fo <- file(output, "wb")
  writeBin(input, fo, size = 4)
  close(fo)
}