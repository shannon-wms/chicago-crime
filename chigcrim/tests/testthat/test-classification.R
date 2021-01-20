# Tests can be run with ctrl + shift + t

test_that("Lostistic regression == glm logit", {
  # Use mtcars data set and classify as automatic/manual
  X = subset(mtcars, select = c("mpg", "wt"))
  y = mtcars$am

  lr <- LogisticRegression$new()
  lr$fit(X, y)

  glm_lr = glm(y ~ ., data=X, family="binomial")

  expect_equal(unname(glm_lr$coefficients), unname(lr$theta), tolerance = 0.05)
})
