# Script contains functions and classes for regression algorithms
#' @importFrom R6 R6Class
#' @import mgcv
#' @importFrom spdep poly2nb
#' @importFrom forcats fct_relevel
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
#' @field n Length of training data.
#' @field pred_vector Vector containing (K + lambda I)^-1 * y_train, the 
#' pre-computed component of the prediction function.
#' @field X_train Training X matrix.
#' @field y_train Training y vector.
#' @field X_test Testing X matrix.
#' @field y_test Testing y vector.
#' @field predictions The predicted values.
#' @export
#' @examples
#' n <- 200
#' n_test <- 100
#' X_train <- matrix(runif(n, -2, 4), nrow = n, ncol = 1)
#' X_test <- matrix(seq(from = -2, to = 4, length.out = n_test), nrow = n_test, ncol = 1)
#' y_train <- as.vector(sin(X_train)) + rnorm(n, sd = 0.3)
#' kr <- KernelRidge$new("rbf", lambda = 1, A = 3)
#' kr$fit(X_train, y_train)
#' y_hat_rbf <- kr$predict(X_test)
KernelRidge <- R6Class("KernelRidge", public = list(
  kernel = "character",
  lambda = "numeric",
  A = "numeric",
  B = "numeric",
  kernel_function = "function",
  n = "integer",
  pred_vector = "vector",
  X_train = "matrix",
  X_test = "matrix",
  y_train = "vector",
  y_test = "vector",
  predictions = "vector",

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
    self$kernel_function <- private$get_kernel_function()
  },

  #' @description
  #' Fits the Gram matrix, an n*n matrix of kernel evaluations on X_train.
  #' @param X_train Training X matrix.
  #' @param y_train Training y vector.
  fit = function(X_train, y_train){
    stopifnot(length(y_train) == nrow(X_train))

    if (self$kernel %in% c("linear", "polynomial")){
      X_train <- cbind(1, X_train)  # RBF doesn't require bias term
    }
    self$y_train <- y_train  # Include only for consistency
    self$X_train <- X_train
    self$n <- nrow(X_train)
    K <- matrix(nrow = self$n, ncol = self$n)
    for (j in 1:self$n){
      K[, j] <- self$kernel_function(X_train[j, ], X_train)
    }
    I <- diag(nrow = self$n, ncol = self$n)
    self$pred_vector <- solve(K + self$lambda*I, self$y_train)
  },

  #' @description
  #' Predicts on an X matrix (test or train)
  #' @param X_test Test X matrix.
  #' @return predictions A vector of predicted values.
  predict = function(X_test){
    if (self$kernel %in% c("linear", "polynomial")){
      X_test <- cbind(1, X_test)  # RBF doesn't need bias term.
    }
    self$X_test <- X_test
    predictions <- c()
    for (i in 1:nrow(X_test)){
      little_k <- self$kernel_function(X_test[i, ], self$X_train)
      pred <- t(little_k) %*% self$pred_vector
      predictions <- c(predictions, pred)
    }
    self$predictions <- predictions
    return(self$predictions)
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
  
  #' Gives the kernel function with hyperparameters set.
  get_kernel_function = function(){
    if (self$kernel == "linear"){
      return(private$linear_kernel)
    }
    
    else if (self$kernel == "polynomial"){
      return(function(x1, X2) private$polynomial_kernel(x1, X2, self$A, self$B))
    }
    
    else if (self$kernel == "rbf"){
      return(function(x1, X2) private$rbf_kernel(x1, X2, self$A))
    }
  },
  
  #' Radial basis function kernel function (compares matrix rows to vector)
  rbf_kernel = function(x1, X2, bandwidth){
    difs <- t(X2)-x1
    dot_prods <- colSums(difs^2)
    exp(-(sqrt(dot_prods)) / (2 * bandwidth^2))
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
#' squared_error_loss(c(1, 2, 3), c(2, 3, 2))
squared_error_loss <- function(y_hat, y){
  mean((y_hat - y)^2)
}

#' Root mean squared error loss
#'
#' @param y_hat Vector of predicted values.
#' @param y Vector of observed values.
#'
#' @return Float root mean squared error loss
#' @export
#'
#' @examples
#' rmse_loss(c(1,2,3), c(2,3,2))
rmse_loss <- function(y_hat, y){
  sqrt(mean((y_hat - y)^2))
}

#' R-Squared
#' @description
#' Calculates the R2 value.
#'
#' @param y_hat Vector of predictions.
#' @param y Vector of observed values.
#'
#' @return Float R2 value between 0 and 1.
#' @export
#'
#' @examples
#' r_squared(c(1,2,3), c(2,3,2))
r_squared <- function(y_hat, y){
  cor(y_hat, y)^2
}

#' Poisson GAM Regression
#'
#' @description
#' R6 class wrapper for the Poisson family of GAMs, fitted via `gam` from the package
#' **mgcv** and tailored for use with a subset of the Chicago Crime dataset. To initialise
#' the object, the user supplies the variables over which to aggregate the data and, if 
#' "community_area" is chosen as `region`, whether to include a neighbourhood list as part
#' of the model. When `fit` is provided with data containing columns corresponding to 
#' `region`, `time_period` and `crime_type` as well as a column "date", it will aggregate
#' counts over these variables. It then fits a GAM with a cyclic cubic regression smoother for
#' `time_period` and treats `region` and `crime_type` as factors (unless a neighbourhood list
#' is included, in which case a Markov random field smoother is applied to `region`). For more
#' details on model fitting methods, see **mgcv**.
#'
#' @field X_train The training data matrix.
#' @field y_train The training response vector.
#' @field X_test The test data matrix.
#' @field time_period Specifies the time period: must be one of "days", "weeks" or
#' "months".
#' @field region Specifies the spatial region to use: must be one of "beat", "district",
#' "ward" or "community".
#' @field crime_type The crime type to include in the model as predictor. If `NULL`, 
#' this will not be included in the model.
#' @field include_nb If "community_area" is the specified region, whether to include a 
#' neighbourhood list and then use a Markov random field smoother for fitting.
#' @field nb_list Neighbourhood list for the community areas.
#' @field gam_fitted The GAM, fitted using package **mgcv**.
#' @field fit_summary Summary of the fitted GAM.
#' @field predictions Vector of predicted response values for the test dataset.
#' 
#' @export
#' 
#' @examples
#' df <- convert_dates(load_data(2019, na_omit = TRUE), filter_week = TRUE)
#' count_df <- get_count_data(df)
#' test_index <- sample(1:nrow(count_df), size = nrow(count_df) %/% 10)
#' X_train <- count_df[-test_index, 1:3]
#' y_train <- count_df[-test_index, 4]
#' X_test <- count_df[test_index, 1:3]
#' y_test <- data[test_index, 4]
#' pg <- PoissonGAM$new(time_period = "week", region = "community_area",
#'                      crime_type = NULL, include_nb = FALSE)
#' pg$fit(X_train, y_train, n_threads = 1)
#' pred <- pg$predict(X_test, y_test)
PoissonGAM <- R6Class("PoissonGAM", public = list(
  X_train = "data.frame",
  X_test = "data.frame",
  y_train = "vector",
  y_test = "vector",
  time_period = "character",
  region = "character",
  crime_type = "character",
  include_nb = "logical",
  nb_list = "list",
  gam_fitted = "list",
  fit_summary = "list",
  predictions = "numeric",
  
  #' @description
  #' Create new PoissonGAM object.
  #' @param time_period One of "week", "month" or "yday" specifying the time period
  #' over which counts are aggregated.
  #' @param region One of "beat" or "community_area", specifying the region over
  #' which counts are aggregated.
  #' @param crime_type crime_type The crime type to include in the model as predictor. 
  #' If `NULL`, this will not be included in the model.
  #' @param include_nb Whether to include a neighbourhood list for the regions,
  #' and hence whether to use a Markov random field smoother for the region.
  #' @param filter_week Whether to filter out the 53rd week.
  initialize = function(time_period = "week", region = "community_area",
                        crime_type = NULL, include_nb = FALSE) {
    self$time_period <- time_period
    self$region <- region
    self$crime_type <- crime_type
    self$include_nb <- include_nb
    
    private$check_params()
    if (include_nb) self$nb_list <- private$get_nb_list()
  },
  
  #' @description
  #' Function for fitting a GAM to the training dataset.
  #' @param X_train The training dataset.
  #' @param convert Whether the column `date` should be converted to instants
  #' or they are already present.
  #' @param n_threads The number of threads to use for parallel smoothing
  #' parameter selection methods in `gam`.
  #' @param ... Additional arguments to be passed to `gam`.
  fit = function(X_train, y_train, n_threads = 1, ...) {
    self$X_train <- X_train
    self$y_train <- y_train
    private$df_train <- cbind(self$X_train, n = self$y_train)
    ctrl <- gam.control(nthreads = n_threads)
    # Construct formula for GAM
    f <- formula(paste("n ~ s(as.numeric(", self$time_period, "), bs = 'cc')"))
    # Add crime type to formula
    if (!is.null(self$crime_type)) {
      f <- formula(paste(deparse(f), "+ ", self$crime_type))
    }
    # Add region and neighbourhood list if required
    if (self$include_nb) {
      f <- formula(paste(deparse(f), "+ s(", self$region,  
                         ", bs = 'mrf', xt = list(nb = self$nb_list))"))
    } else f <- formula(paste(deparse(f), "+", self$region))
    # Fit GAM using mgcv::gam
    self$gam_fitted <- gam(f, data = private$df_train, family = "poisson",
                           control = ctrl, ...)
    self$fit_summary <- summary(self$gam_fitted)
  },
  #' @description
  #' Prediction for a new dataset
  #' @param df_test The training dataset.
  #' @param convert Whether the column `date` should be converted to instants
  #' or they are already present.
  #' @return Vector of predicted values.
  predict = function(X_test = NULL, convert = FALSE, quiet = FALSE) {
    if (!is.null(X_test)) {
      self$X_test <- X_test
      # Prediction on test data
      self$predictions <- predict(self$gam_fitted, newdata = self$X_test, 
                                  type = "response")
    } else { # No test data, predict on training data
      self$predictions <- predict(self$gam_fitted, type = "response")
    }
    if (!quiet) return(self$predictions)
  }
), private = list(
  df_train = "data.frame",
  #' Check parameters are defined properly
  check_params = function() {
    if (!self$time_period %in% c("week", "day", "month")) {
      stop('Time period must be one of "yday", "week" or "month".')
    } else if (!self$region %in% c("beat", "ward", "district", "community_area")) {
      stop('Region must be one of "beat", "ward", "district", or "community_area".')
    } else if (!self$crime_type %in% c("primary_type", "description", "fbi_code")) {
      stop('Crime type must be one of "primary_type", "description", or "fbi_code".')
    } else if (self$include_nb && !(self$region == "community_area")) {
      stop('Neighbourhood list is compatible only with "community_area" type region.')
    }
  },
  #' Construct a neighbourhood list for community areas
  get_nb_list = function() {
    data("community_bounds")
    nb_list <- poly2nb(community_bounds, row.names = community_bounds$community)
    # Name the rows in the list
    names(nb_list) <- attr(nb_list, "region.id")
    return(nb_list)
  }
)
)