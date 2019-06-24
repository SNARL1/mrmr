#' Load a mark-recapture data set
#'
#' This function loads and parses a mark recapture data set. It assumes that
#' three files are available, specifying capture, survey, and (optional)
#' translocation data.
#'
#' @param captures Data frame containing capture-recapture data. Necessary
#' columns include `pit_tag_id` and `survey_date`.
#' @param surveys Data frame containing survey data. Necessary columns include
#' `survey_date`, `primary_period`, and `secondary_period`.
#' @param translocations Optional data frame with translocation data. Necessary
#' columns include `pit_tag_id` and `release_date`. If nothing is provided
#' to this argument, the `clean_data` function assumes that there are no
#' translocations of individuals into the population.
#' @param capture_formula An optional formula specifying the structure of
#' survey-level capture probability covariates. Any variables in this formula
#' must be columns in the `surveys` data frame. The formula must start with
#' `~` and can be provided unquoted, e.g., `capture_formula = ~ temperature`.
#' It is advisable to ensure that any continuous covariates provided in this
#' formula are appropriately scaled (ideally, with mean = 0, and standard
#' deviation = 1).
#' @param survival_formula An optional formula specifying the structure of
#' individual-level survival covariates. Any variables in this formula
#' must be columns in the `translocations` and/or `captures` data frames.
#' The formula must start with `~` and can be provided unquoted.
#' It is advisable to ensure that any continuous covariates provided in this
#' formula are appropriately scaled (ideally, with mean = 0, and standard
#' deviation = 1). Variables specified in this formula cannot be time-varying.
#' They must be fixed for each individual over the entire study.
#' @param survival_fill_value A fill value to use for individual-level
#' covariates. This argument is only required when using the
#' `survival_formula` argument`.
#' @return A list containing the data frames resulting from the capture,
#' translocation, and survey data, along with a list of data formatted for
#' use in a mark recapture model (with name 'stan_d').
#' @examples
#' library(mrmr)
#' library(readr)
#' library(dplyr)
#'
#' captures <- system.file('extdata', 'capture-example.csv',
#'     package = 'mrmr') %>%
#'   read_csv
#' translocations <- system.file('extdata', 'translocation-example.csv',
#'     package = 'mrmr') %>%
#'   read_csv
#' surveys <- system.file('extdata', 'survey-example.csv', package = 'mrmr') %>%
#'   read_csv
#'
#' # read and clean the data using defaults
#' data <- clean_data(captures, surveys, translocations)
#'
#'\dontrun{
#' # (optional) specify a formula for detection probabilities, assuming
#' there is a column called "person_hours"
#' data <- clean_data(captures, surveys, translocations,
#'                    capture_formula = ~ person_hours)
#'}
#' @export
#' @importFrom readr read_csv parse_number
#' @importFrom dplyr %>% mutate group_by summarize ungroup full_join
#' arrange left_join select filter distinct anti_join lead n
#' @importFrom tibble tibble
#' @importFrom tidyr complete separate unite
#' @importFrom reshape2 acast
#' @importFrom rlang .data
#' @importFrom lubridate year
#' @importFrom stats model.matrix

clean_data <- function(captures, surveys,
                       translocations = NA,
                       capture_formula = ~ 1,
                       survival_formula = ~ 1,
                       survival_fill_value = NA) {

  # ---------------------------------
  # captures <- system.file('extdata', 'capture-example.csv',
  #     package = 'mrmr') %>%
  #   read_csv
  # translocations <- system.file('extdata', 'translocation-example.csv',
  #     package = 'mrmr') %>%
  #   read_csv
  # surveys <- system.file('extdata', 'survey-example.csv', package = 'mrmr') %>%
  #   read_csv
  # capture_formula = ~ 1
  # survival_formula = ~treatment
  # survival_fill_value = c(treatment = 'control')
  # ---------------------------------


  any_translocations <- 'data.frame' %in% class(translocations)

  survival_formula_specified <- !identical(survival_formula, ~1)
  survival_fill_named <- !is.null(names(survival_fill_value))
  if (!identical(survival_fill_value, NA) & !survival_formula_specified) {
    stop(paste("A survival fill value was provided, but no survival formula",
               "was specified. The survival_fill_value argument will only",
               "work when a survival formula is specified."))
  }
  if (survival_formula_specified) {
    # use fill values
    if (is.null(names(survival_fill_value))) {
      stop("The survival_fill_value argument must be a named vector, ",
           "e.g., survival_fill_value = c(treatment = 'control').")
    }
    for (n in names(survival_fill_value)) {

      name_in_captures <- n %in% names(captures)
      if (name_in_captures) {
        captures[[n]] <- ifelse(is.na(captures[[n]]),
                                survival_fill_value[n], captures[[n]])
      }

      name_in_translocations <- any_translocations & n %in% names(translocations)
      if (name_in_translocations) {
        translocations[[n]] <- ifelse(is.na(translocations[[n]]),
                                      survival_fill_value[n], translocations[[n]])
      }

      if (!name_in_captures & !name_in_translocations) {
        stop(paste("Named elements in the survival_fill_value argument",
                   "must also be columns in the capture and/or translocation",
                   "data, but", n, "was not found in either."))
      }
    }
  }

  # check that date columns can be parsed as Date objects
  surveys$survey_date <- parse_as_date(surveys$survey_date)
  captures$survey_date <- parse_as_date(captures$survey_date)
  surveys <- surveys %>%
    mutate(secondary_period = ifelse(.data$people == 0,
                                     0, .data$secondary_period),
           primary_period = .data$primary_period + 1) %>%
    group_by(.data$primary_period, .data$secondary_period) %>%
    summarize(survey_date = min(.data$survey_date)) %>%
    ungroup %>%
    mutate(year = year(.data$survey_date),
           is_overwinter = lead(.data$year) - .data$year == 1,
           is_overwinter = ifelse(is.na(.data$is_overwinter),
                                  FALSE, .data$is_overwinter))


  dummy_primary_period <- tibble(primary_period = 1,
                                 secondary_period = 0,
                                 survey_date = min(surveys$survey_date) - 14,
                                 year = min(surveys$year),
                                 is_overwinter = FALSE)
  surveys <- full_join(surveys, dummy_primary_period) %>%
    arrange(.data$primary_period, .data$secondary_period)



  if (any_translocations) {
    translocations$release_date <- parse_as_date(translocations$release_date)
    translocations <- translocations %>%
      mutate(pit_tag_id = as.character(.data$pit_tag_id),
             survey_date = as.Date(.data$release_date)) %>%
      left_join(surveys)
    transloc_tags <- translocations$pit_tag_id
  } else {
    transloc_tags <- c()
  }

  # find M, the superpopulation size
  n_obs <- length(unique(c(transloc_tags, captures$pit_tag_id)))
  n_aug <- n_obs * 2
  M <- n_obs + n_aug

  # find J, the number of secondary periods per primary period
  J <- surveys %>%
    group_by(.data$primary_period) %>%
    summarize(n_sec_periods = max(.data$secondary_period)) %>%
    select(.data$n_sec_periods) %>%
    unlist

  ever_detected <- transloc_tags %in% captures$pit_tag_id
  introduced_but_never_detected <- transloc_tags[!ever_detected]

  y_aug <- tibble(pit_tag_id = c(introduced_but_never_detected,
                                 paste0('aug', 1:n_aug)),
                  y = 1)

  captures$pit_tag_id <- as.character(captures$pit_tag_id)
  y_df <- captures %>%
    mutate(y = 2) %>%
    left_join(surveys) %>%
    select(.data$pit_tag_id, .data$primary_period, .data$secondary_period,
           .data$y) %>%
    unite("survey_id", .data$primary_period, .data$secondary_period,
          sep = '_') %>%
    # fill in implicit 'not detected' values y = 1
    complete(.data$pit_tag_id, .data$survey_id, fill = list(y = 1)) %>%
    separate(.data$survey_id,
             into = c('primary_period', 'secondary_period')) %>%
    mutate(primary_period = parse_number(.data$primary_period),
           secondary_period = parse_number(.data$secondary_period),
           pit_tag_id = as.character(.data$pit_tag_id)) %>%
    arrange(.data$pit_tag_id, .data$primary_period, .data$secondary_period)

  # augment y_df with new individuals
  y_df <- y_df %>%
    unite("survey_identifier", .data$primary_period,
          .data$secondary_period) %>%
    full_join(y_aug) %>%
    complete(.data$pit_tag_id, .data$survey_identifier, fill = list(y = 1)) %>%
    filter(!is.na(.data$survey_identifier)) %>%
    separate(.data$survey_identifier,
             into = c('primary_period', 'secondary_period')) %>%
    mutate(primary_period = parse_number(.data$primary_period),
           secondary_period = parse_number(.data$secondary_period)) %>%
    arrange(.data$pit_tag_id, .data$primary_period, .data$secondary_period)

  stopifnot(!any(is.na(y_df$y)))


  # augment y_df with y = 0 for all primary periods without any surveys
  # with y = 0 representing essentially an NA value
  no_survey_df <- tibble(
    primary_period = filter(surveys,
                            .data$secondary_period == 0)$primary_period,
    y = 0,
    secondary_period = 1)

  # also augment for surveys with no captures
  surveys_with_no_captures <- surveys %>%
    distinct(.data$primary_period, .data$secondary_period) %>%
    anti_join(captures %>% left_join(surveys)) %>%
    filter(.data$secondary_period != 0) %>%
    mutate(y = 1)

  y_df <- y_df %>%
    full_join(no_survey_df) %>%
    unite("survey_id", .data$primary_period, .data$secondary_period,
          sep = "_") %>%
    complete(.data$pit_tag_id, .data$survey_id, fill = list(y = 0)) %>%
    separate(.data$survey_id,
             into = c('primary_period', 'secondary_period')) %>%
    mutate(primary_period = parse_number(.data$primary_period),
           secondary_period = parse_number(.data$secondary_period)) %>%
    full_join(surveys_with_no_captures) %>%
    unite("survey_id", .data$primary_period, .data$secondary_period,
          sep = "_") %>%
    complete(.data$pit_tag_id, .data$survey_id, fill = list(y = 1)) %>%
    separate(.data$survey_id,
             into = c('primary_period', 'secondary_period')) %>%
    mutate(primary_period = parse_number(.data$primary_period),
           secondary_period = parse_number(.data$secondary_period)) %>%
    arrange(.data$pit_tag_id, .data$primary_period, .data$secondary_period) %>%
    filter(!is.na(.data$pit_tag_id))

  if (survival_formula_specified) {
    join_cols <- c('pit_tag_id', names(survival_fill_value))
    y_df <- y_df %>%
      left_join(distinct(captures[, join_cols]))

    if (any_translocations) {
      y_df <- y_df %>%
        left_join(distinct(translocations[, join_cols]))
    }

    # for each survival covariate, fill in missing values
    for (i in seq_along(survival_fill_value)) {
      nam <- names(survival_fill_value)[i]
      y_df[[nam]] <- ifelse(is.na(y_df[[nam]]),
                                  survival_fill_value[i],
                                  y_df[[nam]])
      # make sure there are no NA values after filling
      stopifnot(!any(is.na(y_df[[nam]])))
      # make sure each individual gets just one unique value
      i_counts <- y_df %>%
        group_by(.data$pit_tag_id) %>%
        summarize(nt = length(unique(!!nam)))
      stopifnot(all(i_counts$nt == 1))
    }
  }

  # assert that the only surveys where y = 0 correspond to primary periods with
  # no surveys
  stopifnot(
    all(
      y_df %>%
        filter(.data$y == 0) %>%
        distinct(.data$primary_period) %>%
        unlist == sort(filter(surveys,
                              .data$secondary_period == 0)$primary_period)
    )
  )


  Y <- acast(y_df[, !names(y_df) %in% names(survival_fill_value)],
          pit_tag_id ~ primary_period ~ secondary_period,
          fill = 0, value.var = "y")

  stopifnot(dim(Y)[1] == M)
  stopifnot(all(dimnames(Y)[[1]] == sort(dimnames(Y)[[1]])))
  stopifnot(dim(Y)[2] == max(surveys$primary_period))
  stopifnot(dim(Y)[3] == max(y_df$secondary_period))



  # Process data for introductions ------------------------------------------
  n_intro <- length(unique(transloc_tags))
  stopifnot(all(transloc_tags %in% dimnames(Y)[[1]]))
  introduced <- dimnames(Y)[[1]] %in% transloc_tags
  t_intro <- rep(0, M) # primary period index when the animal was introduced
  # 0 acts as an NA value
  for (i in 1:M) {
    if (introduced[i]) {
      translocate_df_row <- which(transloc_tags == dimnames(Y)[[1]][i])
      t_intro[i] <- translocations$primary_period[translocate_df_row]
    }
  }
  stopifnot(mean(t_intro > 0) == mean(introduced))

  # Deal with detection data ------------------------------------------------
  Jtot <- sum(J)

  # create a vector of indices for each secondary period
  jvec <- c()
  for (i in seq_along(J)) {
    if (J[i] > 0) {
      jvec <- c(jvec, 1:J[i])
    }
  }
  stopifnot(length(jvec) == Jtot)

  survey_number_df <- tibble(primary_period = 1:max(surveys$primary_period)) %>%
    mutate(secondary_idx = 1)

  p_df <- tibble(primary_period = rep(1:max(surveys$primary_period), J),
                 secondary_idx = jvec) %>%
    mutate(j_idx = 1:n())

  j_idx <- survey_number_df %>%
    dplyr::select(.data$primary_period, .data$secondary_idx) %>%
    full_join(p_df) %>%
    acast(primary_period ~ secondary_idx,
          fill = 0,
          value.var = 'j_idx')


  # generate detection design matrix
  X_detect <- model.matrix(capture_formula,
                           data = filter(surveys, .data$secondary_period > 0))
  stopifnot(nrow(X_detect) == Jtot)

  X_surv_df <- y_df[, c('pit_tag_id', names(survival_fill_value))] %>%
    distinct %>%
    arrange(.data$pit_tag_id)
  X_surv <- model.matrix(survival_formula, data = X_surv_df)
  stopifnot(nrow(X_surv) == M)

  # verify that the primary periods with no surveys have all zeros in j_idx
  stopifnot(identical(names(which(rowSums(j_idx) == 0)),
                        which(J == 0) %>% as.character))

  stan_d <- list(M = M,
                 T = max(surveys$primary_period),
                 maxJ = max(J),
                 J = J,
                 Jtot = sum(J),
                 Y = Y,
                 introduced = introduced,
                 t_intro = t_intro,
                 X_detect = X_detect,
                 m_detect = ncol(X_detect),
                 j_idx = j_idx,
                 any_surveys = ifelse(J > 0, 1, 0),
                 prim_idx = rep(1:max(surveys$primary_period), J),
                 m_surv = ncol(X_surv),
                 X_surv = X_surv)

  list(stan_d = stan_d,
       captures = captures,
       translocations = translocations,
       surveys = surveys)
}



parse_as_date <- function(date) {
  if (class(date) != "Date") {
    tryCatch(date <- as.Date(date),
             error = function(c) {
               stop(paste("Couldn't coerce date(s) to a Date object.",
                          "Try formatting date(s) as: %Y-%m-%d,",
                          "or convert all dates to the Date class",
                          "(see ?Date)."))
             }
    )
  }
  todays_date <- format(Sys.time(), "%Y-%m-%d")
  if (any(date > todays_date)) {
    stop("All provided dates must be <= the current date.")
  }
  date
}
