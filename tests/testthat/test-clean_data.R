test_that("clean_data returns the right elements", {
  captures <- system.file("extdata", "capture-example.csv", package = "mrmr")
  translocations <- system.file("extdata", "translocation-example.csv",
                                package = "mrmr")
  surveys <- system.file("extdata", "survey-example.csv", package = "mrmr")
  out <- clean_data(captures, translocations, surveys)
  expected_elements <- c("stan_d", "captures", "translocations", "surveys")
  expect_true(all(expected_elements %in% names(out)))
})
