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

