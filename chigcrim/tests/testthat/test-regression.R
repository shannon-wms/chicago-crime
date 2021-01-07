test_that("KernelRidge", {
  # Check kernel ridge with linear kernel == ridge regression

  # Make toy dataset
  n = 200
  n_test = 100
  X_train = matrix(runif(n, -2,4), nrow=n, ncol=1)
  y_train = as.vector(sin(X_train)) + rnorm(n, sd = 0.3) + 10
  X_test = matrix(seq(from =-2, to=4, length.out = n_test), nrow=n_test, ncol=1)

  kr = KernelRidge$new("linear", lambda = 0.1)
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





