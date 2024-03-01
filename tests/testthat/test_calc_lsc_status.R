test_that("test calc_lsc_status global", {

  timeframe <- as.character(1986:2016)

  test <- calc_status(
    boundary = "lsc",
    config_scenario = system.file(
      "extdata/output/lu_1500_2016/config_lu_1500_2016.json",
      package = "boundaries"
    ),
    config_reference = system.file(
      "extdata/output/pnv_1500_2016/config_pnv_1500_2016.json",
      package = "boundaries"
    ),
    time_span_scenario = timeframe,
    time_span_reference = timeframe,
    spatial_scale = "global",
    savanna_proxy = list(vegc = 7500),
    time_aggregation_args = c(1),
    in_parallel = FALSE,
  ) %>% suppressMessages()

  # test if lsc is the only attribute
  testthat::expect_true(
    all(attributes(test)$names == "lsc")
  )

  thresholds <- list_thresholds("lsc", "steffen2015", "global")

  # test for thresholds attributes in lsc
  testthat::expect_identical(
    attributes(test$lsc)$thresholds,
    thresholds
  )

  # test for expected control variable and class
  testthat::expect_true(
    attributes(test$lsc)$control_variable == "deforestation" &&
      attributes(test$lsc)$class == "control_variable"
  )

  # test for length of time series
  expect_true(
    length(test$lsc) == 31
  )

  # test for expected output
  testthat::expect_true(
    all(
      test$lsc > attributes(test$lsc)$thresholds$pb &
        test$lsc < attributes(test$lsc)$thresholds$highrisk
    )
  )

  # test for as_risk_level
  boundary_status <- as_risk_level(test)

  testthat::expect_true(
    class(boundary_status$lsc) == "boundary_status"
  )

  testthat::expect_true(
    all(
      boundary_status$lsc > attributes(boundary_status$lsc)$thresholds$pb &
        boundary_status$lsc < attributes(boundary_status$lsc)$thresholds$highrisk # nolint
    )
  )

})


test_that("test calc_lsc_status grid", {

  timeframe <- as.character(1986:2016)

  test <- calc_status(
    boundary = "lsc",
    config_scenario = system.file(
      "extdata/output/lu_1500_2016/config_lu_1500_2016.json",
      package = "boundaries"
    ),
    config_reference = system.file(
      "extdata/output/pnv_1500_2016/config_pnv_1500_2016.json",
      package = "boundaries"
    ),
    time_span_scenario = timeframe,
    time_span_reference = timeframe,
    spatial_scale = "grid",
    time_aggregation_args = c(1),
    in_parallel = FALSE,
  ) %>% suppressMessages()

  # test if lsc is the only attribute
  testthat::expect_true(
    all(attributes(test)$names == "lsc")
  )

  thresholds <- list_thresholds("lsc", "steffen2015", "grid")

  # test for expected control variable and class
  testthat::expect_true(
    attributes(test$lsc)$control_variable == "deforestation" &&
      attributes(test$lsc)$class == "control_variable"
  )

  # test for length of time series
  expect_true(
    all(dim(test$lsc) == c(2, 31))
  )

  # test for expected output
  testthat::expect_true(
    all(
      (test$lsc > lpjmlkit::asub(attributes(test$lsc)$thresholds, thresholds = "pb"))[2,] & # nolint
        test$lsc < lpjmlkit::asub(attributes(test$lsc)$thresholds, thresholds = "highrisk") # nolint
    )
  )

  # TODO: Error in thresholds[[, , "holocene"]] : subscript out of bounds
  # test for as_risk_level
  #   boundary_status <- as_risk_level(test)
  # 
  #   testthat::expect_true(
  #     all(
  #       boundary_status$lsc > attributes(boundary_status$lsc)$thresholds$pb &
  #         boundary_status$lsc < attributes(boundary_status$lsc)$thresholds$highrisk # nolint
  #     )
  #   )

})

test_that("test classify_biomes grid", {

  timeframe <- as.character(1986:2016)

  test <- classify_biomes(
    config_reference = system.file(
      "extdata/output/pnv_1500_2016/config_pnv_1500_2016.json",
      package = "boundaries"
    ),
    time_span_reference = timeframe,
    time_aggregation_args = c(1),
    savanna_proxy = list(vegc = NULL), # savanna_proxy = list(pft_lai = 6)
    input_files = list(
      temp = system.file(
        "extdata/input/tas_gswp3-w5e5_1901-2016.clm",
        package = "boundaries"
      ),
      elevation = system.file(
        "extdata/input/elevation.bin",
        package = "boundaries"
      )
    )
  ) %>% suppressMessages()

  # test if lsc is the only attribute
  testthat::expect_true(
    all(names(test) %in% c("biome_id", "biome_names"))
  )

  biome_mapping <- system.file("extdata", "biomes.csv",
                               package = "boundaries") %>%
    readr::read_delim(delim = ";", col_types = readr::cols())

  # test for expected biome names
  testthat::expect_true(
    all(test$biome_names %in% biome_mapping$name)
  )

  # test for dimension of array
  testthat::expect_true(
    all(dim(test$biome_id) == c(2, 31))
  )

  # test if biome id is 5 for these two cells
  #   (Temperate Needleleaved Evergreen Forest)
  testthat::expect_true(
    all(test$biome_id == 5)
  )
})
