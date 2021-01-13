# Script contains functions and classes for regression algorithms
#' @importFrom R6 R6Class
NULL

#' Kernel Ridge Regression
#'
#' @description
#' R6 class for Kernel Ridge Regression. Kernel is chosen using a string in
#' 'linear', 'polynomial' or 'rbf' (for a radial basis function kernel). With a
#' linear kernel, hyperparameters A and B should not be defined. With a
#' polynomial kernel, A and B should be defined  as the constant and degree of
#' the polynomial respectively. With an rbf kernel, A should be defined as the
#' bandwidth, parameter, and B should not be defined. Note that the X matrix
#' should be scaled before use (e.g. to mean 0 and sd=1).
#'
#' A bias column is assumed not to be present in the X matrices, and is
#' added internally.
#'
#' @field kernel String in 'linear', 'polynomial' or 'rbf' for kernel choice.
#' @field lambda Regularisation parameter.
#' @field A First hyperparameter (see description).
#' @field B Second hyperparameter (see description).
#' @field kernel_function The kernel function.
#' @field n Integer length of training data.
#' @field big_k The n by n "K" matrix formed by kernel evaluations on X_train.
#' @field X_train Training X matrix.
#' @field y_train Training y vector.
#' @field X_test Testing X matrix.
#' @field y_test Testing y vector.
#' @field prediction The predicted values.
#' @export
#' @examples
#' n <- 200
#' n_test <- 100
#' X_train <- matrix(runif(n, -2,4), nrow=n, ncol=1)
#' X_test <- matrix(seq(from =-2, to=4, length.out = n_test), nrow=n_test, ncol=1)
#' y_train <- as.vector(sin(X_train)) + rnorm(n, sd = 0.3)
#' kr <- KernelRidge$new("rbf", lambda = 1, 3)
#' kr$fit(X_train = X_train, y_train)
#' y_hat_rbf <- kr$predict(X_test)
KernelRidge = R6Class("KernelRidge", public = list(
  kernel = "character",
  lambda = "numeric",
  A = "numeric",
  B = "numeric",
  kernel_function = "function",
  n = "integer",
  big_k = "matrix",
  X_train = "matrix",
  X_test = "matrix",
  y_train = "vector",
  y_test = "vector",
  prediction = "vector",

  #' @description
  #' Create new KernelRidge object.
  #' @param kernel Kernel choice, 'linear', 'polynomial' or 'rbf'.
  #' @param lambda Regularisation parameter.
  #' @param A Hyperparameters for polynomial and rbf kernels
  #' @param B Hyperparameters for polynomial and rbf kernels
  initialize = function(kernel, lambda, A = numeric(0), B = numeric(0)){
    self$kernel = kernel
    self$lambda = lambda
    self$A = A
    self$B = B

    private$check_params()
    self$kernel_function = self$get_kernel_function()
  },



  #' @description
  #' Gives the kernel function with hyperparameters set.
  get_kernel_function = function(){
    if (self$kernel == "linear"){
      k = private$linear_kernel
      return(k)
    }

    else if (self$kernel == "polynomial"){
      k = function(x1, X2){private$polynomial_kernel(x1, X2, self$A, self$B)}
      return(k)
    }

    else if (self$kernel == "rbf"){
      return(function(x1, X2){private$rbf_kernel(x1, X2, self$A)})
    }
  },

  #' @description
  #' Fits the big K matrix. n by n matrix of kernel evaluations on X_train.
  #' @param X_train Training X matrix.
  #' @param y_train Training y vector
  fit = function(X_train, y_train){
    stopifnot(length(y_train) == nrow(X_train))

    if (self$kernel %in% c("linear", "polynomial")){
      X_train <- cbind(1, X_train)  # RBF doesn't need bias term.
    }

    self$y_train = y_train  # Not really used here, just for consistency
    self$X_train = X_train
    self$n = nrow(X_train)
    K = matrix(nrow = self$n, ncol = self$n)

    for (j in 1:self$n){
      K[, j] = self$kernel_function(X_train[j, ], X_train)
    }
    self$big_k = K
  },

  #' @description
  #' Predicts on an X matrix (test or train)
  #' @param X_test Test X matrix.
  predict = function(X_test){

    if (self$kernel %in% c("linear", "polynomial")){
      X_test <- cbind(1, X_test)  # RBF doesn't need bias term.
    }

    self$X_test = X_test
    I = diag(nrow=self$n, ncol=self$n)
    predictions_ = c()

    for (i in 1:nrow(X_test)){
      little_k = self$kernel_function(X_test[i, ], self$X_train)
      pred = t(little_k) %*% solve(self$big_k + self$lambda*I) %*% self$y_train
      predictions_ = c(predictions_, pred)
    }
    self$prediction = predictions_
    return(predictions_)
  }
), private = list(
  # Checks arguments are defined correctly
  check_params = function(){
    stopifnot(self$kernel %in% c("linear", "polynomial", "rbf"))
    if (self$kernel == "linear" & (length(self$A) != 0 | length(self$B) != 0)){
      stop("With a linear kernel, A and B should not be defined")
    }
    else if (self$kernel == "polynomial" & (length(self$A) == 0 | length(self$B) == 0)){
      stop(paste("With a polynomial kernel, A and B should be defined as the",
                 "constant and degree of the polynomial respectively."))
    }
    else if (self$kernel == "rbf" & (length(self$A) == 0 | length(self$B) != 0)){
      stop(paste("With an rbf kernel, A should be defined as the bandwidth",
                 "parameter, and B should not be defined"))
    }
  },

  #' Radial basis function kernel function (compares matrix rows to vector)
  rbf_kernel = function(x1, X2, bandwidth){
    difs = t(X2)-x1
    dot_prods = colSums(difs^2)
    exp(-(sqrt(dot_prods))/(2*bandwidth^2))
  },

  #' Polynomial kernel function (compares matrix rows to vector)
  polynomial_kernel  = function(x1, X2, constant, degree){
    c((X2 %*% x1 + constant)^degree)
  },

  #' Linear kernel function (compares matrix rows to vector)
  linear_kernel  = function(x1, X2){
    c(X2 %*% x1)
  }
)
)


