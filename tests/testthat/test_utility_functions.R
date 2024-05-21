
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
    all(names(dim(aggr_test)) == c("cell", "month", "year"))
  )

  # check dimensions
  expect_true(
    all(dim(aggr_test) == c(2, 12, 31))
  )
})