# Script contains useful general purpose utility functions
#' @importFrom stats na.omit
#' @importFrom data.table as.data.table
#' @importFrom mltools one_hot
#' @importFrom RSocrata read.socrata
#' @import dplyr
#' @import magrittr
#' @import readr
#' @import lubridate
#' @import doParallel
#' @import parallel
#' @import foreach
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
#' Returns a data frame of Chicago Crime data from specified year or range of years,
#' or returns entire dataset from 2001 to present if no years are specified.
#' @param year integer in 2001:2021 specifying the desired year, or two integers
#' specifying endpoints to filter between.
#' @param strings_as_factors Whether to convert `iucr`, `primary_type`,
#' `description`, `location_description` and `fbi_code` to factors.
#' @param drop_location Whether to drop location as it duplicates data in latitude
#' and longitude.
#' @param na_omit Whether to omit observations with NA values.
#' @param limit Limit set on how many rows of data to obtain.
#' @return tibble data frame of Chicago Crime data.
#' @export
load_data <- function(year = NULL, strings_as_factors = TRUE,
                      drop_location = TRUE, na_omit = FALSE, limit = NULL) {

  # Only accept NULL or valid years
  if (!is.null(year)) {
    if (!all(year %in% 2001:2021)) {
      return("Please choose a year(s) between 2001 and 2021.")
    }
    if (length(year) > 2) {
      return("Please choose either one year or upper and lower values.")
    }
  }
  base_url <- "https://data.cityofchicago.org/resource/ijzp-q8t2.csv"
  # Construct the URL for querying
  if (is.null(year)) { # No year specified, download entire dataset
    full_url <- base_url
  } else if (length(year) == 2) { # Endpoints chosen, filter between them
    year <- as.character(year)
    full_url <- paste0(base_url, "?$where=year between '", year[1], "' and '", year[2], "'")
  } else { # Single year chosen
    year <- as.character(year)
    full_url <- paste0(base_url, "?year=", year)
  }
  if (!is.null(limit)) {
    if(is.na(as.numeric(limit))) return("Limit must be numeric.")
    else { # Convert limit to a character integer
      limit %<>% as.integer() %>% as.character()
      # Add limit query to URL
      full_url <- paste0(full_url, ifelse(is.null(year), "?$", "&$"),
                            "limit=", limit)
    }
  }
  # Pull the data
  df <- read.socrata(full_url, app_token = "avsRhaZaZTqeJOGhBuFJieRzJ") %>% tibble()
  # Omit observations with NA values
  if (na_omit) df %<>% na.omit()
  # Convert dots in column names to underscore
  colnames(df) %<>% strsplit(split = "\\.") %>% lapply(paste, collapse = "_")
  # Convert character columns (excluding case_number and block) to factor or logical
  df %<>% type_convert(col_types = list(arrest = col_logical(),
                                        domestic = col_logical()))
  # Convert strings to factors
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

#' Indexed cross-validation for R6 class
#'
#' @description
#' Computes cross-validation error using specified indexes (rows) of the data set
#' as the test set.
#' @param object R6 class with fit and predict methods.
#' @param X Data matrix
#' @param y Vector of observations.
#' @param error Function to assess error.
#' @param index Vector of indexes corresponding to test data.
#' @param ... Additional arguments passed to fit
#' @return Test error.
#' @export
cv_eval <- function(object, X, y, error, index, ...){
  # Separate X and y into training and test data
  X_train <- X[-index, , drop = FALSE]
  X_test <- X[index, , drop = FALSE]
  y_train <- y[-index]
  y_test <- y[index]
  # Fit training data
  object$fit(X_train, y_train, ...)
  # Predict on test data
  y_hat <- object$predict(X_test)
  # Compute error
  error(y_hat, y_test)
}

#' K-fold cross-validation for R6 class with parallel computation
#'
#' @param object R6 object with fit and predict methods.
#' @param X Data matrix.
#' @param y Response vector.
#' @param error_funcs List of functions for calculating error.
#' @param k Number of folds.
#' @param n_reps Number of repeats.
#' @param parallel Whether to compute in parallel.
#' @param n_threads The number of parallel threads to use. If NULL, this is
#' chosen to be the number of cores minus one.
#' @param ... Additional arguments passed to cv_eval.
#' @return List of length equal to that of `error_funcs` with each element
#' containing a vector of length `n_reps` corresponding to the mean error
#' averaged over `k` folds.
#' @export
kfold_cv <- function(object, X, y, error_funcs, k, n_reps = 1000,
                     parallel = FALSE, n_threads = NULL, ...) {
  n <- nrow(X)
  m <- length(error_funcs)
  if (parallel) { # Parallel computations
    if (is.null(n_threads)) cluster <- makeCluster(detectCores(logical = TRUE) - 1)
    else cluster <- makeCluster(n_threads)
    registerDoParallel(cluster)
    # Load the required packages to the parallel sessions
    clusterEvalQ(cluster, {
      devtools::load_all()
      library(caret)
    })
    # Export the R objects in the current environment to the parallel sessions
    clusterExport(cluster, c("object", "X", "y", "error_funcs", "n", "k"),
                  envir = environment())
    error_list <- foreach (i = 1:n_reps) %dopar% {
      errors <- list()
      folds <- split(sample(1:n), 1:k)
      for (j in 1:k) {
        errors[[j]] <- mapply(cv_eval, error = error_funcs,
                            MoreArgs = list(object = object, X = X, y = y,
                                            index = folds[[j]], ...))
      }
      return(errors)
    }
    stopCluster(cluster)
  } else {
    error_list <- list()
    for (i in 1:n_reps) {
      errors <- list()
      folds <- split(sample(1:n), 1:k)
      for (j in 1:k) {
        errors[[j]] <- mapply(cv_eval, error = error_funcs,
                              MoreArgs = list(object = object, X = X, y = y,
                                              index = folds[[j]], ...))
      }
      error_list[[i]] <- errors
    }
  }
  # Average the list of errors
  mean_errors <- list()
  for (i in 1:m) {
    mean_errors[[i]] <- numeric(n_reps)
    for (j in 1:n_reps) {
      mean_errors[[i]][j] <- lapply(1:k, function(x) error_list[[j]][[x]][i]) %>%
        unlist() %>% mean()
    }
  }
  names(mean_errors) <- names(error_funcs)
  mean_errors
}

#' Convert dates to other formats in dataframe
#'
#' @param df Data frame containing a column with date-time entries to extract
#' instants from.
#' @param date_col The name of the column with date-time entries, defaults to "date".
#' @param as_factors Whether to convert `month`, `week`, `day` and `hour`
#' into factors.
#' @param exclude Specifies which instants to exclude from data frame returned.
#' @param filter_week Whether to filter out weeks > 53, as they  usually have <7 days
#' and can lead to spurious results.
#' @return Data frame `df` with entries of `date_col` converted into instants.
#' @export
convert_dates <- function(df, date_col = "date", as_factors = TRUE, 
                          exclude = NULL, filter_week = FALSE) {
  if (as_factors) { # Convert columns to factors
    df %<>%
      mutate(month = factor(month(get(date_col))),
             week = factor(week(get(date_col))),
             day = factor(day(get(date_col))),
             hour = factor(hour(get(date_col))),
             yday = as.integer(yday(get(date_col))),
             wday = as.integer(wday(get(date_col))),
             date = date(get(date_col)))
  } else { # No factors
    df %<>%
      mutate(month = month(get(date_col)),
             week = week(get(date_col)),
             day = day(get(date_col)),
             hour = hour(get(date_col)),
             yday = as.integer(yday(get(date_col))),
             wday = as.integer(wday(get(date_col))),
             date = date(get(date_col)))
  }
  # Remove specified rows to exclude
  if (!is.null(exclude)) {
    if (!all(exclude %in% c("month", "week", "day", "hour", "yday", "date"))) {
      return('Ensure that exclude is a valid selection from
             "month", "week", "day", "hour", "yday", "date".')
    } else df <- df[, !names(df) %in% exclude]
  }
  # Filter out weeks > 53 if required
  if (!("week" %in% exclude) && filter_week) {
    df %<>% filter(as.integer(week) < 53) %>% mutate(week = factor(week))
  }
  return(df)
}

#' Aggregate given dataset into count data.
#' 
#' @param df Dataset with columns corresponding to `time_period`, `region` and
#' `crime_type` (if desired), as well as an additional column `year`.
#' @param time_period The time period over which to aggregate.
#' @param region The region over which to aggregate.
#' @param crime_type The crime type over which to aggregate, if desired.
#' @return Data frame `count_data` containing count data aggregated over the desired groups.
#' @export
get_count_data = function(df, time_period = "week", region = "community_area", 
                          crime_type = NULL) {
  # Crime type data included
  if (!is.null(crime_type)) {
    count_data <- df %>%
      mutate(!!eval(region) := factor(get(region))) %>%
      count(get(region), get(time_period), 
            year, get(crime_type)) %>%
      rename(!!eval(time_period) := `get(time_period)`,
             !!eval(region) := `get(region)`,
             !!eval(crime_type) := `get(crime_type)`) %>%
      mutate(!!eval(crime_type) := fct_relevel(get(crime_type), sort)) %>%
      arrange(get(region), get(time_period), get(crime_type), year)
  } else { # Crime type data not included
    count_data <- df %>%
      mutate(!!eval(region) := factor(get(region))) %>%
      count(get(region), get(time_period), year) %>%
      rename(!!eval(time_period) := `get(time_period)`,
             !!eval(region) := `get(region)`) %>%
      arrange(get(time_period), get(region), year)
  }
  return(count_data)
}