#' Create an observation matrix
#'
#' This function makes an observation matrix Y from a mark recapture dataset.
#'
#' @param mr The mark recapture data.
#' @param M The number of all zero pseudoindividual detection histories to add.
#' @return A matrix with a row for each individual (+ pseudo-individuals) and a
#' column for each secondary period, where the entries are binary, indicating
#' whether an individual was observed during a secondary period.
#' @export
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate_
#' @importFrom dplyr group_by_
#' @importFrom dplyr arrange_
#' @importFrom dplyr summarize_
#' @importFrom dplyr ungroup
#' @importFrom graphics image
getY <- function(mr, M) {
  stopifnot(class(mr) == "mr")

  # generate columns for primary and secondary periods in survey data
  mr$surveys <- mr$surveys %>%
    mutate_(datetime = ~lubridate::ymd(capture_date, tz = "PT"),
           yday = ~lubridate::yday(datetime)) %>%
    group_by_(~period_id) %>%
    mutate_(primary_period = ~get_period(period_id, 1),
           secondary_period = ~get_period(period_id, 2),
           primary_id = ~paste(survey_year, primary_period, sep = "-"),
           survey_id = ~paste(primary_id, secondary_period, sep = "-")) %>%
    arrange_(~survey_year, ~primary_period, ~secondary_period) %>%
    ungroup()

  # find the number of secondary periods for each primary period
  survey_df <- mr$surveys %>%
    group_by_(~primary_id) %>%
    summarize_(k = ~as.numeric(max(secondary_period)))

  ## define global variables for survey parameters -----------------------
  tags <- c(mr$translocations$pit_tag_id, mr$surveys$pit_tag_id) %>%
    unique()
  I <- length(tags)

  primary_periods <- survey_df$primary_id
  n_primary_periods <- length(primary_periods)
  stopifnot(nrow(survey_df) == n_primary_periods)

  secondary_periods <- unique(mr$surveys$survey_id)
  n_secondary_periods <- length(secondary_periods) # total number of secondary periods

  # first component: introduced individuals
  introduced_tags <- unique(mr$translocations$pit_tag_id)
  if (length(introduced_tags) > 1) {
    n_intro <- length(introduced_tags)
    y_introduced <- array(0, dim = c(n_intro, n_secondary_periods))
    rownames(y_introduced) <- introduced_tags
    for (i in 1:n_intro) {
      for (j in 1:n_secondary_periods) {
        seen <- sum(mr$surveys$pit_tag_id == introduced_tags[i] &
                      mr$surveys$survey_id == secondary_periods[j])
        if (!is.na(seen)) {
          y_introduced[i, j] <- seen
        }
      }
    }
  } else {
    y_introduced <- array(dim = c(0, n_secondary_periods))
  }

  # second component: individuals who recruited naturally
  n_not_introduced <- I - n_intro
  y_natural <- array(0, dim = c(n_not_introduced, ncol(y_introduced)))
  natural_tags <- tags[!(tags %in% introduced_tags)]
  rownames(y_natural) <- natural_tags
  for (i in 1:n_not_introduced) {
    for (j in 1:n_secondary_periods) {
      seen <- sum(mr$surveys$pit_tag_id == natural_tags[i] &
                    mr$surveys$survey_id == secondary_periods[j])
      if (!is.na(seen)) {
        y_natural[i, j] <- seen
      }
    }
  }

  # merge introduced and natural components
  y <- rbind(y_introduced, y_natural)
  colnames(y) <- secondary_periods

  # augment the matrix
  all_zero_detection_histories <- array(0, dim = c(M, ncol(y)))
  y <- rbind(y, all_zero_detection_histories)

  # include metadata
  class(y) <- "obs_mat"
  attr(y, "M") <- M
  attr(y, "I") <- I
  attr(y, "survey_df") <- survey_df
  attr(y, "n_primary_periods") <- n_primary_periods
  attr(y, "n_secondary_periods") <- n_secondary_periods
  attr(y, "mr") <- mr
  y
}

#' @export
plot.obs_mat <- function(x, ...) {
  I <- attr(x, "I")
  image(t(x[I:1, ]), axes = FALSE,
        col = c('white', 'blue3', 'blue', 'green'),
        xlab = 'Timesteps',
        ylab = 'Unique individuals',
        main = 'Observation matrix: all frogs', ...)
}


get_period <- function(x, which_period) {
  stopifnot(which_period == 1 | which_period == 2)
  # extracts primary or secondary period numbers from period_ids
  # which_period is either 1 (primary) or 2 (secondary)
  x <- x %>%
    strsplit("\\.") %>%
    unlist()
  x[which_period + 1]
}