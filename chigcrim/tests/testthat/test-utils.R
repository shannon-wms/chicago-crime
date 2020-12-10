test_that("parse_X gives expected results", {
  # Factors one hot encoded
  df = data.frame(num = c(1,2), X = as.factor(c("a", "b")))
  my_mat = parse_X(df)
  expected_result = as.matrix(data.frame(num=c(1,2), X_a = c(1,0),
                                         X_b = c(0,1)))
  expect_equal(my_mat, expected_result)

  expect_error(parse_X(matrix("a")))
  expect_error(parse_X(data.frame("a")))
})

