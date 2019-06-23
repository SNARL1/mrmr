library(readr)

test_that("clean_data returns the right elements", {
  captures <- read_csv(system.file("extdata", "capture-example.csv",
                                   package = "mrmr"))
  translocations <- read_csv(system.file("extdata", "translocation-example.csv",
                                package = "mrmr"))
  surveys <- read_csv(system.file("extdata", "survey-example.csv",
                                  package = "mrmr"))
  out <- clean_data(captures, surveys, translocations)
  expected_elements <- c("stan_d", "captures", "translocations", "surveys")
  expect_true(all(expected_elements %in% names(out)))
})


test_that("formula specification results in the correct design matrix", {
  captures <- read_csv(system.file("extdata", "capture-example.csv",
                                   package = "mrmr"))
  surveys <- read_csv(system.file("extdata", "survey-example.csv",
                                  package = "mrmr"))
  out <- clean_data(captures, surveys,
                    capture_formula = ~ primary_period)
  expected_names <- c('(Intercept)', 'primary_period')
  expect_true(all(expected_names %in% colnames(out$stan_d$X_detect)))
})

test_that("invalid date formats raise errors", {
  captures <- read_csv(system.file("extdata", "capture-example.csv",
                                   package = "mrmr"))
  surveys <- read_csv(system.file("extdata", "survey-example.csv",
                                  package = "mrmr"))
  translocations <- read_csv(system.file("extdata", "translocation-example.csv",
                                         package = "mrmr"))

  # check for translocation and surveys
  translocations$release_date <- as.integer(translocations$release_date)
  expect_error(clean_data(captures, surveys, translocations), regexp = 'coerce')

  surveys$survey_date <- as.integer(surveys$survey_date)
  expect_error(clean_data(captures, surveys), regexp = 'coerce')
})