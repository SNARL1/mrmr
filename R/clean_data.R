#' Load a mark-recapture data set
#'
#' This function loads and parses a mark recapture data set. It assumes that
#' three files are available, specifying capture, translocation, and survey
#' data.
#'
#' @param captures Character string path to a capture-recapture detection
#' data csv file.
#' @param translocations Character string path to translocation data csv file.
#' @param surveys Character string path to survey data csv file.
#' @return A list containing the data frames resulting from the capture,
#' translocation, and survey data, along with a list of data formatted for
#' use in a mark recapture model (with name "stan_d").
#' @examples
#' captures <- system.file("extdata", "capture-example.csv", package = "mrmr")
#' translocations <- system.file("extdata", "translocation-example.csv",
#'                               package = "mrmr")
#' surveys <- system.file("extdata", "survey-example.csv", package = "mrmr")
#' clean_data(captures, translocations, surveys)
#' @export
#' @importFrom readr read_csv parse_number
#' @importFrom dplyr %>% mutate group_by summarize ungroup if_else full_join
#' arrange left_join select filter distinct anti_join lead n
#' @importFrom tibble tibble
#' @importFrom tidyr complete separate unite
#' @importFrom reshape2 acast
#' @importFrom rlang .data
#' @importFrom lubridate year
#' @importFrom stats model.matrix

clean_data <- function(captures, translocations, surveys) {
  translocations <- read_csv(translocations) %>%
    mutate(pit_tag_id = as.character(.data$pit_tag_id))
  captures <- read_csv(captures)
  surveys <- read_csv(surveys) %>%
    mutate(secondary_period = if_else(.data$people == 0,
                                      0, .data$secondary_period),
           primary_period = .data$primary_period + 1) %>%
    group_by(.data$primary_period, .data$secondary_period) %>%
    summarize(survey_date = min(.data$survey_date),
              duration = sum(.data$survey_duration),
              person_hours = sum(.data$person_hours)) %>%
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


  # find M, the superpopulation size
  n_obs <- length(unique(c(translocations$pit_tag_id, captures$pit_tag_id)))
  n_aug <- n_obs * 2
  M <- n_obs + n_aug

  # find J, the number of secondary periods per primary period
  J <- surveys %>%
    group_by(.data$primary_period) %>%
    summarize(n_sec_periods = max(.data$secondary_period)) %>%
    select(.data$n_sec_periods) %>%
    unlist

  ever_detected <- translocations$pit_tag_id %in% captures$pit_tag_id
  introduced_but_never_detected <- translocations$pit_tag_id[!ever_detected]

  y_aug <- tibble(pit_tag_id = c(introduced_but_never_detected,
                                 paste0('aug', 1:n_aug)),
                  y = 1)

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


  Y <- acast(y_df, pit_tag_id ~ primary_period ~ secondary_period,
             fill = 0,
             value.var = "y")

  stopifnot(dim(Y)[1] == M)
  stopifnot(dim(Y)[2] == max(surveys$primary_period))
  stopifnot(dim(Y)[3] == max(y_df$secondary_period))



  # Process data for introductions ------------------------------------------
  n_intro <- length(unique(translocations$pit_tag_id))
  stopifnot(all(translocations$pit_tag_id %in% dimnames(Y)[[1]]))
  introduced <- dimnames(Y)[[1]] %in% translocations$pit_tag_id
  translocations <- translocations %>%
    mutate(survey_date = as.Date(.data$release_date)) %>%
    left_join(surveys)
  t_intro <- rep(0, M) # primary period index when the animal was introduced
  # 0 acts as an NA value
  for (i in 1:M) {
    if (introduced[i]) {
      translocate_df_row <- which(translocations$pit_tag_id == dimnames(Y)[[1]][i])
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
  X_detect <- model.matrix(~ 1, data = filter(surveys,
                                              .data$secondary_period > 0))
  stopifnot(nrow(X_detect) == Jtot)

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
                 prim_idx = rep(1:max(surveys$primary_period), J))

  list(stan_d = stan_d,
       captures = captures,
       translocations = translocations,
       surveys = surveys)
}
