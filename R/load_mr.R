#' Load a mark-recapture data set
#'
#' This function loads a mark recapture data set from an Excel file.
#'
#' @param path The MS Excel file containing data.
#' @return A data.frame
#' @export
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate_
load_mr <- function(path) {
  sheets <- readxl::excel_sheets(path)
  stopifnot("CaptureHistory" %in% sheets)
  names(sheets) <- tolower(names(sheets))

  # load and validate survey data
  surveys <- readxl::read_excel(path, sheet = "CaptureHistory") %>%
    mutate_(capture_date = ~as.Date(capture_date, origin = "1899-12-30"))
  survey_columns <- c("basin",
                      "site_id",
                      "capture_date",
                      "survey_year",
                      "period_id",
                      "pit_tag_id",
                      "pit_tag_source",
                      "pit_tag_insert",
                      "air_temp",
                      "sun_conditions",
                      "wind_conditions",
                      "survey_duration",
                      "snow_wc",
                      "frog_sex",
                      "frog_weight",
                      "frog_svl")
  validate_cols(survey_columns, names(surveys))

  # load and validate translocation data if it exists
  if ("Translocations" %in% sheets) {
    translocations <- readxl::read_excel(path, sheet = "Translocations") %>%
      mutate_(date = ~as.Date(date, origin = "1899-12-30"))
    transloc_columns <- c("site_id", "date", "pit_tag_id")
    validate_cols(transloc_columns, names(translocations))
  } else {
    translocations <- data.frame(site_id = NULL, date = NULL, pit_tag_id = NULL)
  }
  res <- list(surveys = surveys, translocations = translocations)
  class(res) <- "mr"
  res
}

validate_cols <- function(expected, observed) {
  all_columns_present <- all(expected %in% observed)
  if (!all_columns_present) {
    missing_cols <- expected[!(expected %in% observed)]
    stop(paste("The following expected columns were missing from sheet 1:",
               paste(missing_cols, collapse = ", ")))
  }
}