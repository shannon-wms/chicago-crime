test_that("parse_X gives expected results", {
  # Factors one hot encoded
  df <- data.frame(num = c(1, 2), X = as.factor(c("a", "b")))
  my_mat <- parse_X(df)
  expected <- as.matrix(data.frame(num = c(1, 2), X_a = c(1, 0), X_b = c(0, 1)))
  expect_equal(my_mat, expected)
  
  expect_error(parse_X(matrix("a")))
  expect_error(parse_X(data.frame("a")))
})

test_that("load_data", {
  expect_error(load_data(year = 2030))
  expect_error(load_data(year = c(2012, 2022)))
  expect_error(load_data(year = c(2001, 2003, 2005)))
  # is.na(is.numeric("limit")) should give a warning as well as error
  expect_warning(expect_error(load_data(year = 2019, limit = "limit")))
  # Expect dataframe returned to be the same as querying the API directly and mutating the result
  df <- load_data(year = c(2019, 2020), limit = 50, strings_as_factors = TRUE, drop_location = FALSE)
  expected <- RSocrata::read.socrata("https://data.cityofchicago.org/resource/ijzp-q8t2.csv?$where=year between '2019' and '2020'&$limit=50",
                                     app_token = "avsRhaZaZTqeJOGhBuFJieRzJ") %>% 
    tibble() %>% type_convert(col_types = list(iucr = col_factor(),
                                               primary_type = col_factor(),
                                               description = col_factor(),
                                               location_description = col_factor(),
                                               fbi_code = col_factor()))
  expect_equal(df, expected)
})

test_that("cv_eval works for class Logistic Regression", {
  X <- test_data %>% select(fbi_code, yday)
  y <- test_data$arrest
  # Initialise LogisticRegression object  
  lr <- LogisticRegression$new(solver = "BFGS",
                               control = list(maxit = 1000, reltol = 1e-4),
                               round_y_hat = TRUE)
  error_funcs <- list(accuracy = classification_accuracy, 
                      sensitivity = classification_sensitivity,
                      specificity = classification_specificity)
  index <- sample(1:nrow(X), size = 8)
  # Perform cross-validation
  cv_errors <- cv_eval(lr, X, y, error_funcs, index) %>% unlist()
  # Do the same analysis
  X_train <- X[-index, , drop = FALSE]
  y_train <- y[-index]
  X_test <- X[index, , drop = FALSE]
  y_test <- y[index]
  lr$fit(X_train, y_train)
  y_hat <- lr$predict(X_test)
  expected <- c()
  for (i in 1:length(error_funcs)) expected[i] <- error_funcs[[i]](y_hat, y_test)
  # Expect mean cross-validation errors to be equal
  expect_equal(cv_errors, expected)
})

test_that("convert_dates", {
  expect_error(convert_dates(test_data, exclude = "weeek"))
  expect_error(convert_dates(test_data, date_col = "daate"))
})

test_that("get_count_data", {
  expect_error(get_count_data(test_data, time_period = "weeek"))
  expect_error(get_count_data(test_data, region = "region"))
  # Check the returned dataframe is correct
  df <- test_data %>% get_count_data(time_period = "month", region = "ward", 
                                     crime_type = "fbi_code")
  expected <- test_data %>% select(ward, month, year, fbi_code) %>% 
    count(ward, month, year, fbi_code) %>% 
    mutate(ward = factor(ward), 
           fbi_code = forcats::fct_relevel(fbi_code, sort)) %>%
    arrange(ward, month, year, fbi_code)
  expect_equal(df, expected)
})