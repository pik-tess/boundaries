test_that("test calc_bluewater_status global", {

  timeframe <- as.character(1986:2016)

  test <- calc_status(
    boundary = "bluewater",
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
    approach = list(bluewater = "porkka2024"),
    nyear_window = 1,
    in_parallel = FALSE,
  ) %>% suppressMessages()

  # test if bluewater is the only attribute
  testthat::expect_true(
    all(attributes(test)$names == "bluewater")
  )

  # test for thresholds attributes (0 or manually set to 50)
  testthat::expect_true(
    all(sapply(attributes(test$bluewater)$thresholds, function(x) { # nolint:undesirable_function_linter
      x <- round(x)
      x == 13 || x == 31 || x == 50
    }))
  )

  # test for expected control variable and class
  testthat::expect_true(
    attributes(test$bluewater)$control_variable == "area with wet/dry departures" && # nolint
      attributes(test$bluewater)$class == "control_variable"
  )

  # test for length of time series
  expect_true(
    length(test$bluewater) == 31
  )

  # test for expected output
  testthat::expect_true(
    all(
      test$bluewater < attributes(test$bluewater)$thresholds$highrisk
    )
  )

  # test for as_risk_level
  boundary_status <- as_risk_level(test)

  testthat::expect_true(
    class(boundary_status$bluewater) == "boundary_status"
  )

  testthat::expect_true(
    all(
      boundary_status$bluewater < attributes(boundary_status$bluewater)$thresholds$highrisk
    )
  )
})


test_that("test calc_bluewater_status grid", {

  timeframe <- as.character(1986:2016)

  test <- calc_status(
    boundary = "bluewater",
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
    nyear_window = c(1),
    in_parallel = FALSE,
  ) %>% suppressMessages()

  # test if bluewater is the only attribute
  testthat::expect_true(
    all(attributes(test)$names == "bluewater")
  )

  thresholds <- list_thresholds("bluewater", "gerten2020", "grid")

  # test for thresholds attributes in bluewater
  testthat::expect_identical(
    attributes(test$bluewater)$thresholds,
    thresholds
  )

  # test for expected control variable and class
  testthat::expect_true(
    attributes(test$bluewater)$control_variable == "EFR transgression to uncertainty ratio" && # nolint
      attributes(test$bluewater)$class == "control_variable"
  )

  # test for length of time series
  expect_true(
    all(dim(test$bluewater) == c(2, 31))
  )

  # test for expected output
  testthat::expect_true(
    any(
      test$bluewater > attributes(test$bluewater)$thresholds$highrisk
    )
  )

  # test for as_risk_level
  boundary_status <- as_risk_level(test)

  testthat::expect_true(
    any(
      boundary_status$bluewater > attributes(boundary_status$bluewater)$thresholds$highrisk # nolint
    )
  )
})



test_that("test calc_efrs", {

  timeframe <- as.character(1986:2016)

  discharge <- read_io_format(
    system.file(
      "extdata/output/lu_1500_2016/discharge.bin.json",
      package = "boundaries"
    ),
    timespan = timeframe
  ) %>%
    drop()

  efrs <- calc_efrs(
    discharge,
    approach = "steffen2015"
  )

  # test for length of time series
  expect_true(
    all(dim(efrs) == c(2, 12, 31))
  )

  # test for expected output
  testthat::expect_true(
    all(
      efrs > 0 & efrs < 8
    )
  )

})


test_that("test calc_efrs", {

  timeframe <- as.character(1986:2016)

  discharge <- read_io_format(
    system.file(
      "extdata/output/lu_1500_2016/discharge.bin.json",
      package = "boundaries"
    ),
    timespan = timeframe
  ) %>%
    drop()


  # VMF approach --------------------------------------------------------------- #
  efrs_vmf <- calc_efrs(
    discharge,
    approach = "vmf"
  )

  # test for length of time series
  expect_true(
    all(dim(efrs_vmf) == c(2, 12, 31))
  )

  # test for expected output
  testthat::expect_true(
    all(
      efrs_vmf > 0 & efrs_vmf < 5
    )
  )

  # steffen2015 approach ------------------------------------------------------- #
  efrs_steffen2015 <- calc_efrs(
    discharge,
    approach = "steffen2015"
  )

  # test for length of time series
  expect_true(
    all(dim(efrs_steffen2015) == c(2, 12, 31))
  )

  # test for expected output
  testthat::expect_true(
    all(
      efrs_vmf > 0 & efrs_vmf < 8
    )
  )

  avg_discharge <- aggregate_time(discharge)
  # q90q50 approach ------------------------------------------------------------ #
  efrs_q90q50 <- calc_efrs(
    avg_discharge,
    approach = "q90q50"
  )

  # test for length of time series
  expect_true(
    all(dim(efrs_q90q50) == c(2, 12, 1))
  )

  # test for expected output
  testthat::expect_true(
    all(
      efrs_q90q50 > 0 & efrs_q90q50 < 8
    )
  )
})
