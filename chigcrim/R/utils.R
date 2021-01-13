# Script contains useful general purpose utility functions
#' @importFrom data.table as.data.table
#' @importFrom mltools one_hot
#' @import RSocrata
#' @import dplyr
#' @import magrittr
#' @import readr
#' @import lubridate
NULL

#' Parse matrix into useful format for ML algorithms
#'
#' If a matrix is passed, function checks that it is not a character matrix.
#' If a data frame is passed, the function similarly checks the types, and
#' one hot encodes factors.
#'
#' @param X matrix or data.frame
#' @return matrix
#' @export
parse_X <- function(X){
  if (is.matrix(X) & (typeof(X) == "character")) {
    stop("The type of the matrix should not be character")
  }
  if (inherits(X, "data.frame")){
    dtypes <- sapply(X, typeof)
    allowed_types <- c("factor", "integer", "logical", "double")
    stopifnot(all(dtypes %in% allowed_types))
    is_factor <- sapply(X, is.factor)
    if (any(is_factor)){
      X <- one_hot(as.data.table(X))
    }
  }
  X <- as.matrix(X)
  X
}


#' Download Chicago Crime data using Socrata API
#'
#' @description
#' Returns a data frame of Chicago Crime data from specified year, or entire dataset
#' from 2001 to present if no year is specified.
#' @param year integer in 2001:2020 specifying the desired year.
#' @param strings_as_factors Whether to convert `iucr`, `primary_type`,
#' `description`, `location_description` and `fbi_code` to factors.
#' @param drop_location Whether to drop location as it duplicates data in latitude
#' and longitude.
#' @param na_omit Whether to omit observations with NA values
#' @return tibble Data frame of Chicago Crime data.
load_data <- function(year = NULL, strings_as_factors = TRUE,
                      drop_location = TRUE, na_omit = FALSE) {
  # Only accept valid years
  if(!(year %in% 2001:2020)) return("Please choose a year between 2001 and 2020.")
  base_url <- "https://data.cityofchicago.org/resource/ijzp-q8t2.csv"
  # Download entire dataset if year not specified
  full_url <- ifelse(!is.null(year), paste0(base_url, "?Year=", as.character(year)), base_url)
  df <- read.socrata(full_url, app_token = "avsRhaZaZTqeJOGhBuFJieRzJ") %>% tibble()
  # Omit observations with NA values
  if (na_omit) df %<>% na.omit()
  # Convert dots in column names to underscore
  colnames(df) %<>% strsplit(split = "\\.") %>% lapply(paste, collapse = "_")
  # Convert character columns (excluding case_number and block) to factor or logical
  df %<>% type_convert(col_types = list(arrest = col_logical(),
                                        domestic = col_logical()))
  if (strings_as_factors){
    df %<>% type_convert(col_types = list(iucr = col_factor(),
                                          primary_type = col_factor(),
                                          description = col_factor(),
                                          location_description = col_factor(),
                                          fbi_code = col_factor()))
  }
  # Drop location
  if (drop_location) df %<>% select(-location)
  return(df)
}



#' Convert less common strings to other.
#'
#' @param string_vec Vector of strings.
#' @param n_threshold Threshold count below which will be converted to other.
#' @param print_summary Prints a summary of operation performed.
#'
#' @return vector of strings
#' @export
otherise <- function(string_vec, n_threshold, print_summary = TRUE){
  counts <- table(string_vec)
  other_names <- names(counts[counts < n_threshold])
  string_vec[string_vec %in% other_names] <- "OTHER"
  if (print_summary){
    print(paste(length(other_names), "out of", length(counts),
                "categories were converted to OTHER corresponding to",
                100*length(string_vec[string_vec == "OTHER"])/length(string_vec),
                "% of observations"))
  }
  string_vec
}



#' yday float
#'
#' Convert timestamp to float in 1-366
#' @param timestamp Date in POSIXct
#'
#' @return Floating point numbers
#' @export
#'
yday_float = function(timestamp){
  yday(timestamp) + hour(timestamp)/24 + minute(timestamp)/(24*60)
}

