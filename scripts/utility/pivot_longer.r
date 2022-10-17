library(tidyverse)

data_file <-"/p/projects/open/Jannes/repos/boundaries/inst/extdata/pft_categories.csv" # nolint
data_file <- "/p/projects/open/Jannes/repos/boundaries/inst/extdata/biomes.csv"

# tidy solution for classify_biomes (can be done analogously for GAMMA)
pft_table <- readr::read_csv2(data_file) %>%
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
  # all binary category columns (natural, needle, evergreen) in one categorical
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
               names_prefix = "lpjml_index_npft_")
