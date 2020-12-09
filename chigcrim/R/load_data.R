# Function for downloading Chicago Crime data by year from Socrata API
#' @import RSocrata
#' @import dplyr
#' @import magrittr
NULL

#' Download Chicago Crime data
#' 
#' @description
#' Returns a data frame of Chicago Crime data from specified year, or entire dataset
#' from 2001 to present if no year is specified.
#' @param year integer in 2001:2020 specifying the desired year.
#' @return tibble Data frame of Chicago Crime data.

load_data <- function(year = NULL) {
  # Only accept valid years
  if(!(year %in% 2001:2020)) return("Please choose a year between 2001 and 2020.")
  base_url <- "https://data.cityofchicago.org/resource/ijzp-q8t2.csv"
  # Download entire dataset if year not specified
  full_url <- ifelse(!is.null(year), paste0(base_url, "?Year=", as.character(year)), base_url)
  df <- read.socrata(full_url, app_token = "avsRhaZaZTqeJOGhBuFJieRzJ") %>% tibble()
  # Convert dots in column names to underscore
  colnames(df) %<>% strsplit(split = "\\.") %>% lapply(paste, collapse = "_")
  # Convert character columns (excluding case_number and block) to factor or logical
  df %<>% type_convert(col_types = list(iucr = col_factor(),
                                        primary_type = col_factor(),
                                        description = col_factor(),
                                        location_description = col_factor(),
                                        arrest = col_logical(),
                                        domestic = col_logical()))
  return(df)
}
