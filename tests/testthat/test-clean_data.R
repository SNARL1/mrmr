library(readr)
library(dplyr)

captures <- read_csv(system.file("extdata", "capture-example.csv",
                                 package = "mrmr"))
translocations <- read_csv(system.file("extdata", "translocation-example.csv",
                                       package = "mrmr"))
surveys <- read_csv(system.file("extdata", "survey-example.csv",
                                package = "mrmr"))


test_that("clean_data returns the right elements", {
  out <- clean_data(captures, surveys, translocations)
  expected_elements <- c("stan_d", "captures", "translocations", "surveys")
  expect_true(all(expected_elements %in% names(out)))
})


test_that("formula specification results in the correct design matrix", {
  out <- clean_data(captures, surveys,
                    capture_formula = ~ primary_period)
  expected_names <- c('(Intercept)', 'primary_period')
  expect_true(all(expected_names %in% colnames(out$stan_d$X_detect)))
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
  expect_error(clean_data(captures, surveys, survival_formula = ~treatment,
                          survival_fill_value = c(foobar = "control")),
               regexp = "must also be columns")
})

test_that("survival_fill_value fills values in the capture columns", {
  d <- clean_data(captures, surveys,
                  survival_formula = ~ treatment,
                  survival_fill_value = c(treatment = "filled_value"))
  pr_na_in_data <- captures %>%
    distinct(pit_tag_id, treatment) %>%
    summarize(mean(is.na(treatment))) %>%
    unlist %>%
    as.numeric

  pr_filled_in_cleaned_data <- d$survival_covariate_df %>%
    filter(pit_tag_id %in% captures$pit_tag_id) %>%
    summarize(mean(treatment == "filled_value")) %>%
    unlist %>%
    as.numeric

  expect_identical(pr_na_in_data, pr_filled_in_cleaned_data)
})


test_that("survival_fill_value fills values in the translocation columns", {
  d <- clean_data(captures, surveys, translocations,
                  survival_formula = ~ treatment,
                  survival_fill_value = c(treatment = "filled_value"))

  individuals <- translocations %>%
    distinct(pit_tag_id, treatment) %>%
    full_join(distinct(captures, pit_tag_id, treatment))

  # be sure that each individual is only assigned one treatment
  expect_identical(length(unique(individuals$pit_tag_id)),
                   nrow(individuals))
  # Then check that the fraction of filled values is the same
  pr_na_in_data <- mean(is.na(individuals$treatment))

  pr_filled_in_cleaned_data <- d$survival_covariate_df %>%
    filter(!grepl("^aug", pit_tag_id)) %>%
    summarize(mean(treatment == "filled_value")) %>%
    unlist %>%
    as.numeric

  expect_identical(pr_na_in_data, pr_filled_in_cleaned_data)
})

test_that("surivival formulas create correct design matrices", {
  data <- clean_data(captures, surveys, translocations,
                     survival_formula = ~ treatment,
                     survival_fill_value = c(treatment = "wild-caught"))
  X_surv <- data$stan_d$X_surv
  expect_true("treatmentwild-caught" %in% colnames(X_surv))
})

test_that("dead captures raise errors", {
  captures <- system.file("extdata",
                          "equid-captures.csv",
                          package = "mrmr") %>%
    read_csv
  surveys <- system.file("extdata",
                         "equid-surveys.csv",
                         package = "mrmr") %>%
    read_csv
  removals <- system.file("extdata",
                          "equid-removals.csv",
                          package = "mrmr") %>%
    read_csv

  expect_error(clean_data(captures, surveys, removals = removals),
               regexp = "including dead animals encountered on surveys")

  # check case sensitivity
  capitalized_captures <- captures %>%
    mutate(capture_animal_state = tools::toTitleCase(capture_animal_state))
  expect_error(clean_data(capitalized_captures, surveys, removals = removals),
               regexp = "including dead animals encountered on surveys")
})

test_that("duplicate captures raise errors", {
  captures <- rbind(captures[1, ], captures)

  expect_error(clean_data(captures, surveys),
               regexp = "Duplicate entries")
})

test_that("recruitment warning is printed when no natural recruits", {
  captures <- captures %>%
    filter(pit_tag_id %in% translocations$pit_tag_id)

  expect_warning(d <- clean_data(captures, surveys, translocations),
               regexp = "natural recruitment")
  expect_equal(d$stan_d$any_recruitment, 0)
})

test_that("Misnamed removal date columns raise errors", {
  captures <- system.file("extdata",
                          "equid-captures.csv",
                          package = "mrmr") %>%
    read_csv %>%
    filter(capture_animal_state != "dead")
  surveys <- system.file("extdata",
                         "equid-surveys.csv",
                         package = "mrmr") %>%
    read_csv
  removals <- system.file("extdata",
                          "equid-removals.csv",
                          package = "mrmr") %>%
    read_csv %>%
    rename(remove_date = removal_date)

  expect_error(clean_data(captures, surveys, removals = removals),
               regexp = "was not found in the removal data")
})

test_that("Misnamed tag id columns raise errors", {
  captures <- system.file("extdata",
                          "equid-captures.csv",
                          package = "mrmr") %>%
    read_csv %>%
    filter(capture_animal_state != "dead")
  surveys <- system.file("extdata",
                         "equid-surveys.csv",
                         package = "mrmr") %>%
    read_csv
  removals <- system.file("extdata",
                          "equid-removals.csv",
                          package = "mrmr") %>%
    read_csv %>%
    rename(pit_tag_ref = pit_tag_id)

  expect_error(clean_data(captures, surveys, removals = removals),
               regexp = "was not found in the removal data")
})

test_that("Duplicate survey dates raise errors", {
  surveys$survey_date[2] <- surveys$survey_date[1]
  expect_error(clean_data(captures, surveys),
               regexp = "dates are duplicated")
})