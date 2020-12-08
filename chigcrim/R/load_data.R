load_data <- function(year = NULL) {
  # Only accept valid years
  stopifnot(year %in% 2001:2020)
  base_url <- "https://data.cityofchicago.org/resource/ijzp-q8t2.csv"
  # Download entire dataset if year not specified
  full_url <- ifelse(!is.null(year), paste0(base_url, "?Year=", as.character(year)), base_url)
  df <- read.socrata(full_url, app_token = "avsRhaZaZTqeJOGhBuFJieRzJ") %>% tibble()
  # Edit column names
  colnames(df) %<>% strsplit(split = "\\.") %>% lapply(paste, collapse = "_")
  # Convert character columns (excluding Case_Number and Block) to factor or logical
  df %>% type_convert(col_types = list(IUCR = col_factor(),
                                       Primary_Type = col_factor(),
                                       Description = col_factor(),
                                       Location_Description = col_factor(),
                                       Arrest = col_logical(),
                                       Domestic = col_logical()))
  return(df)
}