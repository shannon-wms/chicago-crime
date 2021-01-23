# Script contains functions and classes for classification algorithms
#' @importFrom stats cor
#' @importFrom R6 R6Class
#' @importFrom e1071 svm
#' @import ggplot2
#' @importFrom randomForest randomForest
NULL

#' Sigmoid function
#'
#' Computes value of sigmoid function at `z`.
#'
#' @param z numeric Log ratio estimate.
#' @return float
#' @export
#' @examples sigmoid(0)
sigmoid <- function(z) 1 / (1 + exp(-z))

#' Logistic Regression.
#'
#' @description R6 class for Logistic Regression
#'
#' @field lambda Regularisation parameter (default = 0).
#' @field solver Solver to use in `optim`: one of "BFGS", "CG" or "L-BFGS-B".
#' @field X_train Training data matrix.
#' @field y_train Training response vector.
#' @field X_test Test data matrix.
#' @field theta The fitted parameters of the logistic regression model.
#' @field control List of control parameters to be passed to `optim`` control.
#' @field rounding Logical, whether to round predicted values to 0 and 1.
#' @field predictions The predicted values for `X_test`.
#' @export
#' @examples
#' X <- subset(mtcars, select = c("mpg", "wt"))
#' y <- mtcars$am
#' lr <- LogisticRegression$new()
#' lr$fit(X, y)
#' lr$predict(X)
LogisticRegression <- R6Class("LogisticRegression", public = list(
  lambda = "numeric",
  solver = "character",
  X_train = "matrix",
  y_train = "vector",
  X_test = "matrix",
  theta = "vector",
  control = "list",
  rounding = "logical",
  predictions = "vector",

  #' @description
  #' Create new `LogisticRegression` object.
  #' @param solver `"L-BFGS-B"`, `"BFGS"` or `"CG"`. Default is `"L-BFGS-B"`.
  #' @param lambda Regularisation parameter, defualt 0.
  #' @param control List of control parameters to be passed to `optim` control.
  #' @param rounding Whether to round predicted values.
  initialize = function(solver = "L-BFGS-B", lambda = 0,
                        control = NULL, rounding = FALSE){
    stopifnot(solver %in% c("BFGS", "CG", "L-BFGS-B"))
    stopifnot(lambda >= 0)
    stopifnot(is.logical(rounding))

    self$lambda <- lambda
    self$solver <- solver
    self$control <- control
    self$rounding <- rounding
  },
  
  #' @description
  #' Fit the object to training data X and y.
  #' @param X Training data (dataframe or matrix).
  #' @param y Training data (vector).
  #' @param ... Additional arguments passed to optim.
  fit = function(X, y, ...){
    self$X_train <- parse_X(X)
    self$y_train <- y
    # Add bias term to data matrix for fitting
    X <- cbind(bias = 1, self$X_train)
    # Initialise solution
    init_theta <- rep(0, ncol(X))
    optim_res <- optim(init_theta, fn = private$log_loss, gr = private$grad_log_loss,
                       method = self$solver, X = X, y = self$y_train, 
                       control = self$control, ...)
    if (optim_res$convergence == 1) warning("The algorithm did not converge.")
    # Obtain logistic regression parameters
    theta <- optim_res$par
    names(theta) <- colnames(X)
    self$theta <- theta
    invisible(self)
  },

  #' @description
  #' Predict on X.
  #' @param X X training or testing X.
  #' @param quiet If TRUE, the function will not return the predicted values and
  #' will only update the object.
  predict = function(X, quiet = FALSE) {
    self$X_test <- parse_X(X)
    X <- cbind(bias = 1, self$X_test)
    predictions <- as.vector(sigmoid(X %*% self$theta))
    names(predictions) <- rownames(X)
    if (self$rounding) predictions <- round(predictions)
    self$predictions <- predictions
    if (!quiet) return(self$predictions)
  }
), private = list(
  #' Computes log loss (with regularisation). Assumes one-hot-encoded
  #' categorical variables and bias term in X.
  log_loss = function(theta, X, y) {
    n <- length(y)
    y_hat <- sigmoid(X %*% theta)
    y_hat[which(y_hat == 1)] <- y_hat[which(y_hat == 1)] - 1e-15  # Fix numerical precision errors
    cost <- - (1 / n) * sum((y * log(y_hat) + (1 - y) * log(1 - y_hat)))
    reg <- (self$lambda / 2) * sum(theta[-1]^2)
    return(cost + reg)
  },
  
  #' Computes the gradient of the log loss function (with regularisation).
  #' Assumes one-hot-encoded categorical variables and bias term in X.
  grad_log_loss = function(theta, X, y) {
    n <- length(y)
    y_hat <- sigmoid(X %*% theta)
    grad <- (1 / n) * (t(X) %*% (y_hat - y)) + c(0, self$lambda * theta[-1])
    return(grad)
  }
))


#' Support vector machine wrapper
#'
#' @description
#' R6 class wrapper for support vector machine implementation from **e1071** package.
#'
#' @field kernel Kernel function to use: one of `"radial"`, `"linear"`, `"polynomial"`
#' or `"sigmoid"`.
#' @field X_train Training data matrix.
#' @field y_train Training response vector.
#' @field X_test Test data matrix. 
#' @field fitted_model `svm` object representing the fitted model.
#' @field predictions Vector of predicted response values for the test dataset.
#' @export
SupportVectorMachine <- R6Class("SupportVectorMachine", public = list(
  kernel = "character",
  X_train = "matrix",
  y_train = "vector",
  X_test = "matrix",
  fitted_model = "list",
  predictions = "vector",

  #' @description
  #' Create new SupportVectorMachine object.
  #' @param kernel Kernel function to use.
  initialize = function(kernel = "radial"){
    stopifnot(kernel %in% c("radial", "linear", "polynomial", "sigmoid"))
    self$kernel <- kernel
  },

  #' @description
  #' Fit the object to training data `X` and `y` via `svm` from **optim**.
  #' @param X Training data matrix.
  #' @param y Training response vector.
  #' @param ... Additional arguments to be passed to `svm`.
  fit = function(X, y, ...){
    self$X_train <- X
    self$y_train <- y
    df <- cbind(X, y)
    self$fitted_model <- svm(y ~ ., data = df, kernel = self$kernel,
                             type = "C-classification", cachesize = 120, ...)
    invisible(self)
  },

  #' @description
  #' Predict on X.
  #' @param X Test data matrix.
  #' @param logical Whether to store the predictions as logical values.
  #' @param quiet If TRUE, the function will not return the predicted values and
  #' will only update the object.
  predict = function(X, logical = TRUE, quiet = FALSE) {
    self$X_test <- X
    predictions <- predict(self$fitted_model, newdata = X)
    if (logical) self$predictions <- as.logical(predictions)
    else self$predictions <- predictions
    if(!quiet) return(self$predictions)
  }
))


#' Random forest
#' R6 class wrapper for random forest. Uses randomForest package.
#' 
#' @field X_train Training data matrix.
#' @field y_train Training response vector.
#' @field X_test Test data matrix. 
#' @field fitted_model Fitted randomForest object.
#' @field predictions Predicted response values for test dataset.
#' @export
RandomForest <- R6Class("RandomForest", list(
  X_train = "matrix",
  y_train = "vector",
  X_test = "matrix",
  fitted_model = "list",
  predictions = "vector",
  #' @description
  #' Create new RandomForest object.
  initialize = function() {
  },

  #' @description
  #' Fit the object to training data X and y.
  #' @param X Training data matrix.
  #' @param y Training response vector.
  #' @param ... Additional arguments to be passed to `randomForest`.
  fit = function(X, y, ...) {
    self$X_train <- X
    self$y_train <- y
    self$fitted_model <- randomForest(self$X_train, y = as.factor(self$y_train))
    invisible(self)
  },

  #' @description
  #' Predict on X.
  #' @param X Test data matrix.
  #' @param logical Whether to return the predictions as logical values.
  #' @param quiet If TRUE, the function will not return the predicted values and
  #' will only update the object.
  #' @return vector Vector of predicted values.
  predict = function(X, logical = TRUE, quiet = FALSE) {
    self$X_test <- X
    predictions <- predict(self$fitted_model, newdata = X)
    if (logical) self$predictions <- as.logical(predictions)
    else self$predictions <- predictions
    if (!quiet) return(self$predictions)
  }
))

#' Accuracy
#' 
#' @description
#' Computes the proportion of correctly classified predictions.
#' @param y_hat Binary vector of predictions.
#' @param y_test Binary vector of observations.
#'
#' @return float Proportion of correct predictions
#' @export
#'
#' @examples
#' classification_accuracy(c(0, 1, 0), c(0, 1, 1))
classification_accuracy <- function(y_hat, y_test){
  sum(y_hat == y_test) / length(y_hat)
}

#' Sensitivity
#' 
#' @description
#' Calculates the sensitivity (proportion of positives correctly identified).
#'
#' @param y_hat Binary vector of predictions.
#' @param y_test Binary vector of observations.
#'
#' @return float Sensitivity.
#' @export
#'
#' @examples
#' classification_sensitivity(c(0, 1, 0), c(0, 1, 1))
classification_sensitivity <- function(y_hat, y_test){
  sum(y_hat[y_test == 1] == 1) / length(y_hat[y_test == 1])
}

#' Specificity
#' 
#' @description
#' Calculates the specificity (proportion of negatives correctly identified).
#'
#' @param y_hat Binary vector of predictions.
#' @param y_test Binary vector of observations.
#'
#' @return float Specificity
#' @export
#'
#' @examples
#' classification_specificity(c(0, 1, 0), c(0, 1, 1))
classification_specificity <- function(y_hat, y_test){
  sum(y_hat[y_test == 0] == 0) / length(y_hat[y_test == 0])
}



