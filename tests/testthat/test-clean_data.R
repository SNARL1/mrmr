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

test_that("unnamed arguments to survival_fill_value raise errors", {
  expect_error(clean_data(captures = NA, surveys = NA,
                          survival_formula = ~ treatment,
                          survival_fill_value = c('treatment')),
               regexp = "must be a named")
})

test_that("providing survival_fill_value argument only works with formula", {
  expect_error(clean_data(captures = NA, surveys = NA,
                          survival_fill_value = c(treatment = "control")),
               regexp = "when a survival formula is specified")
})

test_that("providing a name in survival_fill_value w/out column match error", {
  captures <- read_csv(system.file("extdata", "capture-example.csv",
                                   package = "mrmr"))
  surveys <- read_csv(system.file("extdata", "survey-example.csv",
                                  package = "mrmr"))
  expect_error(clean_data(captures, surveys, survival_formula = ~treatment,
                          survival_fill_value = c(foobar = "control")),
               regexp = "must also be columns")
})

test_that("survival_fill_value fills values in the capture columns", {
  captures <- read_csv(system.file("extdata", "capture-example.csv",
                                   package = "mrmr"))
  surveys <- read_csv(system.file("extdata", "survey-example.csv",
                                  package = "mrmr"))
  d <- clean_data(captures, surveys, survival_formula = ~ treatment,
                  survival_fill_value = c(treatment = "filled_value"))
  pr_na <- mean(is.na(captures$treatment))
  pr_filled <- mean(d$captures$treatment == "filled_value")
  expect_identical(pr_na, pr_filled)
})


test_that("survival_fill_value fills values in the translocation columns", {
  captures <- read_csv(system.file("extdata", "capture-example.csv",
                                   package = "mrmr"))
  surveys <- read_csv(system.file("extdata", "survey-example.csv",
                                  package = "mrmr"))
  translocations <- read_csv(system.file("extdata", "translocation-example.csv",
                                         package = "mrmr"))

  d <- clean_data(captures, surveys, translocations,
                  survival_formula = ~ treatment,
                  survival_fill_value = c(treatment = "filled_value"))
  pr_na <- mean(is.na(translocations$treatment))
  pr_filled <- mean(d$translocations$treatment == "filled_value")
  expect_identical(pr_na, pr_filled)
})

test_that("surivival formulas create correct design matrices", {
  captures <- read_csv(system.file("extdata", "capture-example.csv",
                                   package = "mrmr"))
  surveys <- read_csv(system.file("extdata", "survey-example.csv",
                                  package = "mrmr"))
  translocations <- read_csv(system.file("extdata", "translocation-example.csv",
                                         package = "mrmr"))

  data <- clean_data(captures, surveys, translocations,
                     survival_formula = ~ treatment,
                     survival_fill_value = c(treatment = "wild-caught"))
  X_surv <- data$stan_d$X_surv
  expect_true("treatmentwild-caught" %in% colnames(X_surv))
})