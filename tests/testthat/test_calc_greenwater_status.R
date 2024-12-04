test_that("test greenwater_status global", {

  timeframe <- as.character(1986:2016)

  test <- calc_status(
    boundary = "greenwater",
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
    approach = list(greenwater = "porkka2024"),
    time_series_avg = c(1),
    in_parallel = FALSE,
  ) %>% suppressMessages()

  # test if greenwater is the only attribute
  testthat::expect_true(
    all(attributes(test)$names == "greenwater")
  )

  # test for thresholds attributes (0 or manually set to 50)
  testthat::expect_true(
    round(attributes(test$greenwater)$thresholds$holocene) == 8 &&
      round(attributes(test$greenwater)$thresholds$pb) == 29
  )

  # test for expected control variable and class
  testthat::expect_true(
    attributes(test$greenwater)$control_variable == "area with wet/dry departures" && # nolint
      attributes(test$greenwater)$class == "control_variable"
  )

  # test for length of time series
  expect_true(
    length(test$greenwater) == 31
  )

  # test for expected output
  testthat::expect_true(
    sum(
      test$greenwater < attributes(test$greenwater)$thresholds$pb
    ) >= 9
  )
})

test_that("test greenwater_status regional", {

  timeframe <- as.character(1986:2016)

  test <- calc_status(
    boundary = "greenwater",
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
    spatial_scale = "regional",
    approach = list(greenwater = "porkka2024"),
    time_series_avg = c(1),
    in_parallel = FALSE,
  ) %>% suppressMessages()

  # test if greenwater is the only attribute
  testthat::expect_true(
    all(attributes(test)$names == "greenwater")
  )

  # test for thresholds attributes (0 or manually set to 50)
  testthat::expect_true(
    all(round(lpjmlkit::asub(attributes(test$greenwater)$thresholds, thresholds = "holocene"), 2) == 8.33) &&
      all(round(lpjmlkit::asub(attributes(test$greenwater)$thresholds, thresholds = "pb"), 2) == 33.33)
  )

  # test for expected control variable and class
  testthat::expect_true(
    attributes(test$greenwater)$control_variable == "area with wet/dry departures" && # nolint
      attributes(test$greenwater)$class == "control_variable"
  )

  # test for length of time series
  expect_true(
    all(dim(test$greenwater) == c(31, 2))
  )

  # test for expected output
  testthat::expect_true(
    all(
      test$greenwater[1, ] > lpjmlkit::asub(attributes(test$greenwater)$thresholds, thresholds = "holocene") # nolint 
    ) && # nolint
      # test for almost all (1 is FALSE) > pb
      all(test$greenwater[2,] >= lpjmlkit::asub(attributes(test$greenwater)$thresholds, thresholds = "pb")) # nolint 
  )
})
