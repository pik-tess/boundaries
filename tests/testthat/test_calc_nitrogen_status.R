test_that("test calc_nitrogen_status global", {

  timeframe <- as.character(1986:2016)

  test <- calc_status(
    boundary = "nitrogen",
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
    approach = list(nitrogen = "schulte_uebbing2022"),
    time_aggregation_args = c(1),
    in_parallel = FALSE,
  ) %>% suppressMessages()

  # test if nitrogen is the only attribute
  testthat::expect_true(
    all(attributes(test)$names == "nitrogen")
  )

  thresholds <- list_thresholds("nitrogen", "schulte_uebbing2022", "global")

  # test for thresholds attributes in lsc
  testthat::expect_identical(
    attributes(test$nitrogen)$thresholds,
    thresholds
  )

  # test for expected control variable and class
  testthat::expect_true(
    attributes(test$nitrogen)$control_variable == "nitrogen surplus on agricultural land" && # nolint
      attributes(test$nitrogen)$class == "control_variable"
  )

  # test for length of time series
  expect_true(
    length(test$nitrogen) == 31
  )

  # test for expected output
  # TODO: missing simulation data - new output files only inlcude dummy data
  testthat::expect_true(
    all(
      round(test$nitrogen, digit = 2) %in% c(0.01, 0.02)
    )
  )

  # test for as_risk_level
  boundary_status <- as_risk_level(test)

  testthat::expect_true(
    class(boundary_status$nitrogen) == "boundary_status"
  )

})


test_that("test calc_nitrogen_status grid", {

  timeframe <- as.character(1986:2016)

  test <- calc_status(
    boundary = "nitrogen",
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

  # test if nitrogen is the only attribute
  testthat::expect_true(
    all(attributes(test)$names == "nitrogen")
  )

  thresholds <- list_thresholds("nitrogen", "braun2022", "grid")

  # test for expected control variable and class
  testthat::expect_true(
    attributes(test$nitrogen)$control_variable == "nitrogen leaching in runoff to surface water" && # nolint
      attributes(test$nitrogen)$class == "control_variable"
  )

  # test for length of time series
  expect_true(
    all(dim(test$nitrogen) == c(2, 31))
  )

  # test for expected output
  testthat::expect_true(
    all(test$nitrogen[1, ] < attributes(test$nitrogen)$thresholds$pb) &&
      all(test$nitrogen[2, ] > attributes(test$nitrogen)$thresholds$pb)
  )

  # test for as_risk_level
  boundary_status <- as_risk_level(test)

  testthat::expect_true(
    class(boundary_status$nitrogen) == "boundary_status"
  )

  testthat::expect_true(
    all(
      boundary_status$nitrogen[1,] < attributes(boundary_status$nitrogen)$thresholds$pb & # nolint
        boundary_status$nitrogen[2,] > attributes(boundary_status$nitrogen)$thresholds$pb # nolint
    )
  )

})



test_that("test calc_nitrogen_status grid (minusref)", {

  timeframe <- as.character(1986:2016)

  test <- calc_status(
    boundary = "nitrogen",
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
    approach = list(nitrogen = "braun2022_minusref"),
    spatial_scale = "grid",
    time_aggregation_args = c(1),
    in_parallel = FALSE,
  ) %>% suppressMessages()

  # test for expected output
  testthat::expect_true(
    all(test$nitrogen[1, ] < attributes(test$nitrogen)$thresholds$pb) &&
      all(test$nitrogen[2, ] > attributes(test$nitrogen)$thresholds$pb)
  )
})
