test_that("KernelRidge", {
  # Check kernel ridge with linear kernel == ridge regression

  # Make toy dataset
  n <- 200
  n_test <- 100
  X_train <- matrix(runif(n, -2,4), nrow = n, ncol = 1)
  y_train <- as.vector(sin(X_train)) + rnorm(n, sd = 0.3) + 10
  X_test <- matrix(seq(from = -2, to = 4, length.out = n_test), nrow = n_test, ncol = 1)

  kr <- KernelRidge$new("linear", lambda = 0.1)
  kr$fit(X_train = X_train, y_train)
  y_hat_linear <- kr$predict(X_test)

  # Ridge regression
  X_train <- cbind(1, X_train)
  X_test <- cbind(1, X_test)
  XXt_lamda_I <- crossprod(X_train)  + diag(0.1, nrow = ncol(X_train))
  Xyt <- crossprod(X_train, y_train)
  w <- solve(XXt_lamda_I, Xyt)
  result <- c(X_test %*% w)

  expect_equal(result, y_hat_linear)

  # Rough tests that other kernel functions work as expected
  my_mat <- matrix(c(1, 2, 3, 4), nrow=2)

  kr <- KernelRidge$new("polynomial", lambda = 1, 2, 2)
  similarity <- kr$kernel_function(my_mat[1,], my_mat)
  expect_equal(similarity, c(144, 256))

  kr <- KernelRidge$new("rbf", lambda = 1, 1)
  similarity <- kr$kernel_function(my_mat[1,], my_mat)
  expect_equal(similarity[1], 1)
})

test_that("PoissonGAM", {
  # Make a limited dataset from Chicago crime data
  data("test_data")
  # Initialise PoissonGAM object
  pg <- PoissonGAM$new(time_period = "week", region = "community_area", 
                       crime_type = "fbi_code", include_nb = TRUE,
                       filter_week = TRUE)
  pg$fit(df_train = test_data, convert = TRUE, n_threads = 7)
  pg$predict()
  # Do the same analysis
  data("community_bounds")
  nbd_list <- spdep::poly2nb(community_bounds, row.names = community_bounds$community)
  names(nbd_list) <- attr(nbd_list, "region.id")
  test_data %<>% convert_dates(exclude = "hour") %>% 
    filter(as.integer(week) < 53) %>%
    mutate(week = factor(week))
  count_data <- test_data %>%
    mutate(community_area = factor(community_area)) %>%
    count(week, community_area, year, fbi_code) %>%
    mutate(fbi_code = forcats::fct_relevel(fbi_code, sort)) %>%
    relocate(community_area, week) %>%
    arrange(community_area, week, fbi_code)
  # Expect the count data passed to the models to be equal
  expect_equal(count_data, pg$count_train)
  gam_fit <- mgcv::gam(n ~ s(as.numeric(week), bs = "cc") + 
                     s(community_area, bs = "mrf", xt = list(nb = nbd_list)) +
                     fbi_code, 
                   data = count_data, family = "poisson",
                   control = gam.control(nthreads = 7))
  gam_predict <- predict(gam_fit, type = "response")
  # Expect the fitted model summaries (except formula which will never be equal) to be equal
  expect_equal(summary(gam_fit)[-12], summary(pg$gam_fitted)[-12])
  # Expect the predicted values to be equal
  expect_equal(pg$predictions, gam_predict)
})
