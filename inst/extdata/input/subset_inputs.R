library(lpjmlkit)
library(tibble)

output_path <- "." # "/p/projects/open/Jannes/repos/boundaries/inst/extdata/input" # nolint
base_path <- "." # "/p/projects/lpjml/input/historical" # nolint

input_files <- c(
  # temperature
  file.path(base_path, "GSWP3-W5E5/tas_gswp3-w5e5_1901-2016.clm"),
  # precipitation
  file.path(base_path, "ISIMIP3av2/obsclim/GSWP3-W5E5/pr_gswp3-w5e5_obsclim_1901-2019.clm"), # nolint
  # drainage
  file.path(base_path, "input_VERSION2/drainagestn.bin"),
  # elevation
  file.path(base_path, "input_VERSION2/elevation.bin")
)

for (input_file in input_files) {

  # output name for new input
  output_name <- basename(input_file)

  # create header
  input_header <- read_header(filename = input_file)
  input_header$header["version"] <- 4
  input_header$header["ncell"] <- 2
  input_header$header["firstcell"] <- 27410
  input_header$header["datatype"] <- 4
  input_header$header["scalar"] <- 1


  # read, subset and format input data
  # daily climate input
  if (input_header$header[["nbands"]] == 365) {
    input <- read_io(
      input_file,
      subset = list(cell = as.character(c(27410, 27411))),
      nstep = 365,
      nbands = 1
    ) |>
      as_array() |>
      drop() |>
      t() |>
      as.vector()
    input_header$header["nstep"] <- 365
    input_header$header["nbands"] <- 1

    # drainage input
  } else if (input_header$name == "LPJDRAI") {
    input <- read_io(
      input_file,
      subset = list(cell = as.character(c(27410, 27411))),
      datatype = 2
    ) |>
      as_array() |>
      drop() |>
      t() |>
      as.vector()

    # elevation input
  } else if (input_header$name == "LPJELEV") {
    input <- read_io(
      input_file,
      subset = list(cell = as.character(c(27410, 27411))),
    ) |>
      as_array() |>
      drop()
  }

  # write new input data to file with new header
  output <-  file.path(output_path, output_name)
  write_header(filename = output, header = input_header, overwrite = TRUE)
  fo <- file(output, "ab")
  writeBin(input, fo)
  close(fo)
}