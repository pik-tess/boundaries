
test_that("test aggregate_time", {

  timeframe <- as.character(1986:2016)

  test <- read_io_format(
    system.file(
      "extdata/output/lu_1500_2016/rootmoist.bin.json",
      package = "boundaries"
    ),
    timespan = timeframe
  ) %>%
    drop()

  # interpolate -------------------------------------------------------------- #
  aggr_test <- aggregate_time(
    test,
    nyear_window = 10,
    interpolate = TRUE
  )

  # check dimnames
  expect_true(
    all(names(dim(aggr_test)) == c("cell", "month", "year"))
  )

  # check dimensions
  expect_true(
    all(dim(aggr_test) == c(2, 12, 31))
  )

  # with time window of 5 years ---------------------------------------------- #
  aggr_test <- aggregate_time(
    test,
    nyear_window = 5
  )

  # check dimnames
  expect_true(
    all(names(dim(aggr_test)) == c("cell", "month", "window"))
  )

  # check dimensions
  expect_true(
    all(dim(aggr_test) == c(2, 12, 6))
  )

  # with time window of 5 years and moving average --------------------------- #
  aggr_test <- aggregate_time(
    test,
    nyear_window = 5,
    moving_average = TRUE
  )

  # check dimnames
  expect_true(
    all(names(dim(aggr_test)) == c("cell", "month", "year"))
  )

  # check dimensions
  expect_true(
    all(dim(aggr_test) == c(2, 12, 31))
  )


  # with time window of 5 years, moving average and nyear_reference ---------- #
  aggr_test <- aggregate_time(
    test,
    nyear_window = 5,
    moving_average = TRUE,
    nyear_reference = 30
  )
  # check dimnames
  expect_true(
    all(names(dim(aggr_test)) == c("cell", "month", "year"))
  )

  # check dimensions
  expect_true(
    all(dim(aggr_test) == c(2, 12, 31))
  )

  # with nyear_window and nyear_reference ------------------------------------ #
  aggr_test <- aggregate_time(
    test,
    nyear_window = 1,
    nyear_reference = 30
  )

  # check dimnames
  expect_true(
    all(names(dim(aggr_test)) == c("cell", "month", "window"))
  )

  # check dimensions
  expect_true(
    all(dim(aggr_test) == c(2, 12, 30))
  )

})