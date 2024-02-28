test_that("test calc_greenwater_status global", {

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
    method = list(greenwater = "porkka2023"),
    time_aggregation_args = c(1),
    in_parallel = FALSE,
  ) %>% suppressMessages()

  # test if greenwater is the only attribute
  testthat::expect_true(
    all(attributes(test)$names == "greenwater")
  )

  # test for thresholds attributes (0 or manually set to 50)
  testthat::expect_true(
    attributes(test$greenwater)$thresholds$holocene == 0 &&
      attributes(test$greenwater)$thresholds$pb == 100
  )

  # test for expected control variable and class
  testthat::expect_true(
    attributes(test$greenwater)$control_variable == "area with wet/dry departures (%)" && # nolint
      attributes(test$greenwater)$class == "control_variable"
  )

  # test for length of time series
  expect_true(
    length(test$greenwater) == 31
  )

  # test for expected output
  testthat::expect_true(
    all(
      test$greenwater < attributes(test$greenwater)$thresholds$pb
    )
  )
})

test_that("test calc_greenwater_status subglobal", {

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
    spatial_scale = "subglobal",
    method = list(greenwater = "porkka2023"),
    time_aggregation_args = c(1),
    in_parallel = FALSE,
  ) %>% suppressMessages()

  # test if greenwater is the only attribute
  testthat::expect_true(
    all(attributes(test)$names == "greenwater")
  )

  # test for thresholds attributes (0 or manually set to 50)
  testthat::expect_true(
    all(attributes(test$greenwater)$thresholds$holocene == 100) &&
      all(attributes(test$greenwater)$thresholds$pb == 400)
  )

  # test for expected control variable and class
  testthat::expect_true(
    attributes(test$greenwater)$control_variable == "area with wet/dry departures (%)" && # nolint
      attributes(test$greenwater)$class == "control_variable"
  )

  # test for length of time series
  expect_true(
    all(dim(test$greenwater) == c(2, 31))
  )

  # test for expected output
  testthat::expect_true(
    all(
      test$greenwater[1,] > attributes(test$greenwater)$thresholds$holocene[1] && # nolint
        # test for almost all (1 is FALSE) > pb
        sum(test$greenwater[2,] >= attributes(test$greenwater)$thresholds$pb[2]) == 30 # nolint 
    )
  )
})
