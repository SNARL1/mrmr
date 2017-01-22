#' Create a partly observed occurrence matrix
#'
#' This function makes an occurrence matrix Z from an observation matrix.
#'
#' @param Y The observation matrix, of class obs_mat.
#' @return A matrix with a row for each individual (+ pseudo-individuals) and a
#' column for each primary period, where the entries are binary, indicating
#' whether an individual was present during a primary period. Unknown states
#' are represented by NA entries.
#' @export
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate_
#' @importFrom dplyr group_by_
#' @importFrom dplyr arrange_
#' @importFrom dplyr summarize_
#' @importFrom dplyr filter_
#' @importFrom dplyr select_
#' @importFrom dplyr distinct_
#' @importFrom dplyr full_join
#' @importFrom graphics image
#'
getZ <- function(Y) {
  stopifnot(class(Y) == "obs_mat")
  mr <- attr(Y, "mr")
  survey_df <- attr(Y, "survey_df")

  times <- mr$surveys %>%
    filter_(~secondary_period == '1') %>%
    select_(~capture_date) %>%
    distinct_(~capture_date) %>%
    mutate_(source = ~'observations',
            date = ~capture_date) %>%
    select_(~date, ~source)

  # determine whether there are introduction timesteps
  introduced_tags <- unique(mr$translocations$pit_tag_id)

  if (length(introduced_tags) > 0) {
    trans_times <- mr$translocations %>%
      select_(~date) %>%
      distinct_(~date) %>%
      mutate_(source = ~'translocations')
    times <- full_join(times, trans_times) %>%
      arrange_(~as.Date(date, format = "%m/%d/%Y"))
  }
  times <- arrange_(times, ~as.Date(date, format = "%m/%d/%Y"))

  times$index <- 1:nrow(times)
  stopifnot(nrow(times) != length(survey_df$primary_id))
  introduction_timesteps <- which(times$source == 'translocations')

  # construct the occurrence matrix
  n_primary_periods <- attr(Y, "n_primary_periods")
  n_z_cols <- n_primary_periods + length(introduction_timesteps)
  z <- array(NA, dim = c(nrow(Y), n_z_cols))
  rownames(z) <- rownames(Y)
  colnames(z) <- rep(NA, ncol(z))
  intro_dates <- unique(as.character(mr$translocations$date))
  stopifnot(length(introduction_timesteps) == length(intro_dates))

  new_z_column <- 1
  old_z_column <- 1
  intro_index <- 1
  while (new_z_column <= ncol(z)) {
    if (new_z_column %in% introduction_timesteps) {
      colnames(z)[new_z_column] <- intro_dates[intro_index]
      intro_index <- intro_index + 1
    } else {
      colnames(z)[new_z_column] <- survey_df$primary_id[old_z_column]
      # move to next primary period only if one used as column name
      old_z_column <- old_z_column + 1
    }
    new_z_column <- new_z_column + 1
  }
  if (length(intro_dates) > 0) {
    stopifnot(all(intro_dates %in% colnames(z)))
  }

  # fill in matrix entries
  # make a vector that indexes the surveys by primary period
  primary_match <- rep(NA, ncol(Y))
  for (i in 1:ncol(Y)) {
    colname <- colnames(Y)[i]
    primary_period <- colname %>%
      strsplit(split = "-") %>%
      unlist() %>%
      `[`(-3) %>%
      paste(collapse = "-")
    primary_match[i] <- which(colnames(z) == primary_period)
  }

  # verify that the introduction timesteps are not matched to detection histories
  if (length(introduction_timesteps) > 0) {
    stopifnot(!(introduction_timesteps %in% primary_match))
  }

  # find first and last observations and complete corresponding parts of z
  for (i in 1:nrow(Y)) {
    n_observations <- sum(Y[i, ] > 0, na.rm = TRUE)
    if (n_observations > 0) {
      first_observed <- which(Y[i, ] > 0)[1]
      last_observed <- which(Y[i, ] > 0)[n_observations]
      # find the timesteps (primary periods) associated with the first and last surveys
      first_timestep <- primary_match[first_observed]
      last_timestep <- primary_match[last_observed]

      # fill in known occurrence states
      z[i, first_timestep:last_timestep] <- 1
    }
  }

  # fill in introduction states
  if (length(introduction_timesteps) > 0) {
    t_intro <- rep(NA, nrow(mr$translocations))
    n_introductions <- length(unique(mr$translocations$date))
    n_added <- rep(0, n_introductions)
    for (i in seq_along(mr$translocations$pit_tag_id)) {
      # find row in z
      row_index <- which(rownames(z) == mr$translocations$pit_tag_id[i])

      # find column in z
      which_introduction <- which(colnames(z)[introduction_timesteps] ==
                                    as.character(mr$translocations$date[i]))

      col_index <- introduction_timesteps[which_introduction]
      t_intro[i] <- col_index
      n_added[which_introduction] <- n_added[which_introduction] + 1
      if (col_index == 1) {
        z[i, 1] <- 1
      } else {
        z[i, 1:(t_intro[i] - 1)] <- 0
        z[i, t_intro[i]] <- 1
      }
      # they were present at least through the last observation
      if (sum(Y[i, ], na.rm = T) > 0) {
        last_observed <- rev(which(Y[i, ] > 0))[1]
        last_timestep <- primary_match[last_observed]
        z[i, (which.max(z[i, ])[1]):last_timestep] <- 1
      }
    }
  }

  class(z) <- "occ_mat"
  attr(z, "Y") <- Y
  attr(z, "times") <- times
  z
}

#' @export
plot.occ_mat <- function(x, ...) {
  image(t(x[nrow(x):1, ]),
        xlab = 'Timesteps',
        ylab = 'Unique individuals',
        main = 'Partly observed state matrix',
        axes = FALSE,
        col = c('grey20', 'dodgerblue'), ...)
}
