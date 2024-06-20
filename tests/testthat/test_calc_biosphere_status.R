test_that("test biosphere_status global", {

  timeframe <- as.character(1986:2016)

  test <- calc_status(
    boundary = "biosphere",
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
    time_series_avg = c(1),
    in_parallel = FALSE,
    savanna_proxy = list(vegc = 7500),
    path_baseline = system.file(
      "extdata/output/pnv_1500_2016/",
      package = "boundaries"
    )
  ) %>% suppressMessages()

  # test if biosphere is the only attribute
  testthat::expect_true(
    all(attributes(test)$names == "biosphere")
  )

  thresholds <- list_thresholds("biosphere", "stenzel2023", "grid")

  # test for thresholds attributes in biosphere
  testthat::expect_identical(
    attributes(test$biosphere)$thresholds,
    thresholds
  )

  # test for expected control variable and class
  testthat::expect_true(
    attributes(test$biosphere)$control_variable == "BioCol" && # nolint
      attributes(test$biosphere)$class == "control_variable"
  )

  # test for length of time series
  expect_true(
    length(test$biosphere) == 31
  )

  # test for expected output
  testthat::expect_true(
    # almost all values are above the highrisk threshold
    mean(
      round(test$biosphere, digit = 1) >= attributes(test$biosphere)$thresholds$highrisk # nolint
    ) > 0.9
  )


  # test for as_risk_level
  boundary_status <- as_risk_level(test)

  testthat::expect_true(
    class(boundary_status$biosphere) == "boundary_status"
  )

  testthat::expect_true(
    # almost all values are above the highrisk threshold
    mean(
      boundary_status$biosphere > attributes(boundary_status$biosphere)$thresholds$highrisk
    ) > 0.9
  )
})


test_that("test biosphere_status grid", {

  timeframe <- as.character(1986:2016)

  test <- calc_status(
    boundary = "biosphere",
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
    time_series_avg = 1,
    in_parallel = FALSE,
    savanna_proxy = list(vegc = 7500),
    path_baseline = system.file(
      "extdata/output/pnv_1500_2016/",
      package = "boundaries"
    )
  ) %>% suppressMessages()

  # test if biosphere is the only attribute
  testthat::expect_true(
    all(attributes(test)$names == "biosphere")
  )

  thresholds <- list_thresholds("biosphere", "stenzel2023", "grid")

  # test for thresholds attributes in biosphere
  testthat::expect_identical(
    attributes(test$biosphere)$thresholds,
    thresholds
  )

  # test for expected control variable and class
  testthat::expect_true(
    attributes(test$biosphere)$control_variable == "BioCol" && # nolint
      attributes(test$biosphere)$class == "control_variable"
  )

  # test for length of time series
  expect_true(
    all(dim(test$biosphere) == c(2, 31))
  )

  # test for expected output
  testthat::expect_true(
    all(
      test$biosphere > attributes(test$biosphere)$thresholds$pb
    )
  )

  # test for as_risk_level
  boundary_status <- as_risk_level(test)

  testthat::expect_true(
    all(boundary_status$biosphere[1,] > attributes(boundary_status$biosphere)$thresholds$pb) && # nolint
    all(boundary_status$biosphere[2,] > attributes(boundary_status$biosphere)$thresholds$highrisk) # nolint
  )
})



test_that("test biosphere_status subglobal", {

  timeframe <- as.character(1986:2016)

  test <- calc_status(
    boundary = "biosphere",
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
    spatial_scale = "subglobal",
    time_series_avg = 1,
    in_parallel = FALSE,
    savanna_proxy = list(vegc = 7500),
    path_baseline = system.file(
      "extdata/output/pnv_1500_2016/",
      package = "boundaries"
    )
  ) %>% suppressMessages()


  # test for length of time series
  expect_true(
    all(dim(test$biosphere) == c(2, 31))
  )

  # test for expected output
  testthat::expect_true(
    mean(test$biosphere > attributes(test$biosphere)$thresholds$highrisk) > 0.9
  )

  # test for as_risk_level
  boundary_status <- as_risk_level(test)

  testthat::expect_true(
    mean(boundary_status$biosphere > attributes(boundary_status$biosphere)$thresholds$highrisk) > 0.9 # nolint
  )

})
