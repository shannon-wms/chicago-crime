# Script contains functions and classes for regression algorithms
#' @importFrom R6 R6Class
#' @import mgcv
#' @import spdep
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
#' @field inv_big_k_lambda_y_train Matrix containing (K + lambda I)^-1 . y_train
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
#' kr$fit(X_train, y_train)
#' y_hat_rbf <- kr$predict(X_test)
KernelRidge <- R6Class("KernelRidge", public = list(
  kernel = "character",
  lambda = "numeric",
  A = "numeric",
  B = "numeric",
  kernel_function = "function",
  n = "integer",
  inv_big_k_lambda_y_train = "matrix",
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
    self$kernel <- kernel
    self$lambda <- lambda
    self$A <- A
    self$B <- B

    private$check_params()
    self$kernel_function <- self$get_kernel_function()
  },



  #' @description
  #' Gives the kernel function with hyperparameters set.
  get_kernel_function = function(){
    if (self$kernel == "linear"){
      k <- private$linear_kernel
      return(k)
    }

    else if (self$kernel == "polynomial"){
      k <- function(x1, X2){private$polynomial_kernel(x1, X2, self$A, self$B)}
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

    self$y_train <- y_train  # Not really used here, just for consistency
    self$X_train <- X_train
    self$n <- nrow(X_train)
    K <- matrix(nrow = self$n, ncol = self$n)

    for (j in 1:self$n){
      K[, j] <- self$kernel_function(X_train[j, ], X_train)
    }
    I <- diag(nrow=self$n, ncol=self$n)
    self$inv_big_k_lambda_y_train <- solve(K + self$lambda*I,
                                           self$y_train)
  },

  #' @description
  #' Predicts on an X matrix (test or train)
  #' @param X_test Test X matrix.
  predict = function(X_test){

    if (self$kernel %in% c("linear", "polynomial")){
      X_test <- cbind(1, X_test)  # RBF doesn't need bias term.
    }

    self$X_test <- X_test
    I <- diag(nrow<-self$n, ncol=self$n)
    predictions <- c()

    for (i in 1:nrow(X_test)){
      little_k <- self$kernel_function(X_test[i, ], self$X_train)
      pred <- t(little_k) %*% solve(self$big_k + self$lambda*I) %*% self$y_train
      predictions <- c(predictions, pred)
    }
    self$prediction <- predictions
    return(predictions)
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
    difs <- t(X2)-x1
    dot_prods <- colSums(difs^2)
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

#' Mean squared error loss
#'
#' @param y_hat Vector of predictions
#' @param y Vector of observed values
#'
#' @return Float of the mean squared error loss
#' @export
#'
#' @examples
#' squared_error_loss(c(1,2,3), c(2,3,2))
squared_error_loss <- function(y_hat, y){
  mean((y_hat - y)^2)
}

#' Poisson GAM Regression
#' 
#' @description 
#' R6 class for the Poisson family of GAMs.
#' 
#' @field df_train The training dataset.
#' @field df_test The test dataset.
#' @field time_period Specifies the time period: must be one of "days", "weeks" or
#' "months".
#' @field region Specifies the spatial region to use: must be one of "beat", "district",
#' "ward" or "community".
#' @field include_nb Whether to include a neighbours list for the discrete spatial 
#' regions and then use a Markov random field smoother for fitting.
#' @field include_crimetype Whether to include crime types as a predictor in the model.
#' @field n_threads Number of threads to use when fitting the GAM in parallel.
PoissonGAM <- R6Class("PoissonGAM", public = list(
  df_train = "data.frame",
  df_test = "data.frame",
  time_period = "character",
  region = "character",
  include_nb = "logical",
  include_crimetype = "logical",
  n_threads = "integer",
  nbd_list = "list",
  count_train = "data.frame",
  count_test = "data.frame",
  gam_fitted = "list",
  predictions = "numeric",
  
  #' @description 
  #' Create new PoissonGAM object.
  initialize = function(time_period = "week", region = "community_area", 
                        include_nb = FALSE, include_crimetype = FALSE) {
    self$time_period <- time_period
    self$region <- region
    self$include_nb <- include_nb
    self$include_crimetype <- include_crimetype
    
    private$check_params()
    if (include_nb) self$nbd_list <- private$get_nbd_list(region)
  },
  
  #' @description
  #' Function for fitting a GAM to the training dataset.
  #' @param df_train The training dataset.
  #' @param n_threads The number of threads to use for parallel smoothing
  #' parameter selection methods in `gam`
  #' @param ... Additional arguments to be passed to `gam`
  fit = function(df_train, n_threads = 1, ...) {
    self$count_train <- count_data(df_train)
    ctrl <- gam.control(nthreads = n_threads)
    # Construct formula for GAM
    f <- formula(n ~ s(as.numeric(get(self$time_period)), bs = "cc"))
    if (self$include_nb) {
      f %<>% update(~ . + s(get(self$region), bs = "mrf", xt = list(nb = self$nbd_list)))
    } else f %<>% update(~ . + get(self$region))
    # Add crime type to formula
    if (self$include_crimetype) f %<>% update(~ . + fbi_code)
    # Fit GAM
    self$gam_fitted <- gam(f, data = self$count_train, family = "poisson",
                           control = ctrl, ...)
  },
  #' @description 
  #' Function for predicting new dataset 
  predict = function(df_test = NULL, ...) {
    if (!is.null(df_test)) {
      self$count_test <- count_data(df_test)
      gam_predict <- predict(self$gam_fitted, newdata = self$count_test)
    } else gam_predict <- predict(self$gam_fitted)
    self$predictions <- gam_predict
  }
), private = list(
  # Check parameters are defined properly
  check_params = function() {
    if (!self$time_period %in% c("week", "day", "month")) {
      stop('Time period must be one of "yday", "week" or "month".')
    } else if (!self$region %in% c("beat", "community_area")) {
      stop('Region must be one of "beat" or "community_area".')
    }
  },
  # Aggregate given dataset into count data
  get_count_data = function(df) {
    # Crime type data included
    if (self$include_crimetype) {
      # Ensure FBI code is sorted alphabetically
      df$fbi_code %<>% forcats::fct_relevel(sort)
      count_data <- df %>% 
        mutate(!!eval(self$region) := factor(get(self$region))) %>%
        count(get(self$region), get(self$time_period), fbi_code) %>%
        rename(!!eval(self$time_period) := `get(self$time_period)`,
               !!eval(self$region) := `get(self$region)`) %>%
        arrange(eval(self$time_period), eval(self$region), fbi_code)
    } else { # Crime type data not included
      count_data <- df %>% 
        mutate(!!eval(self$region) := factor(get(self$region))) %>%
        count(get(self$region), get(self$time_period)) %>%
        rename(!!eval(self$time_period) := `get(self$time_period)`,
               !!eval(self$region) := `get(self$region)`) %>%
        arrange(eval(self$time_period), eval(self$region))
    }
  },
  # Construct a neighbourhood list for beats or community areas
  get_nbd_list = function() {
    if (self$region == "community") {
      bounds <- data(beat_bounds)
      nbd_list <- spdep::poly2nb(bounds, row.names = bounds$beat)
    }
    else {
      bounds <- data(community_bounds)
      nbd_list <- spdep::poly2nb(bounds, row.names = bounds$community)
    }
    # Name the rows in the list
    names(nbd_list) <- attr(nbd_list, "region.id")
    return(nbd_list)
  }
)
)