# Script contains functions and classes for classification algorithms
#' @importFrom R6 R6Class
#' @import ggplot2
NULL

#' Sigmoid function
#'
#' Computes sigmoid function using 1/(1 + exp(-z)).
#'
#' @param z numeric Log ratio estimate.
#' @return float
#' @export
#' @examples sigmoid(0)
sigmoid <- function(z){
  1/(1 + exp(-z))
}

#' Logistic Regression.
#'
#' @description R6 class for Logistic Regression
#'
#' @field lambda Regularisation parameter (defualt=0).
#' @field solver "BFGS", "CG" or "L-BFGS-B".
#' @field X Training X (dataframe or matrix).
#' @field y Training y vector.
#' @field theta The fitted parameters.
#' @export
#' @examples
#' X <- subset(mtcars, select = c("mpg", "wt"))
#' y <- mtcars$am
#' lr <- LogisticRegression$new()
#' lr$fit(X, y)
#' lr$predict(X)
LogisticRegression <- R6Class("LogisticRegression", list(
  lambda = NULL,
  solver = NULL,
  X = NULL,
  y = NULL,
  theta = NULL,

  #' @description
  #' Create new LogisticRegression object.
  #' @param solver "L-BFGS-B", "BFGS" or "CG". Default "L-BFGS-B".
  #' @param lambda regularisation parameter, defualt 0.
  initialize = function(solver = "L-BFGS-B", lambda=0){
    stopifnot(solver %in% c("BFGS", "CG", "L-BFGS-B"))
    stopifnot(lambda >= 0)

    self$lambda <- lambda
    self$solver <- solver
  },

  #' @description
  #' Computes log loss (with regularisation). Assumes one-hot-encoded
  #' categorical variables and bias term in X.
  #' @param theta Coefficients.
  #' @param X matrix with first column as bias.
  #' @param y 0,1 vector.
  log_loss = function(theta, X, y){
    n <- length(y)
    y_hat <- sigmoid(X %*% theta)
    y_hat[which(y_hat == 1)] <- y_hat[which(y_hat == 1)] - 1e-15  # Fix numerical precision errors
    cost <- -1/n * sum((y*log(y_hat) + (1-y)*log(1-y_hat)))
    reg <- self$lambda/2 * sum(theta[-1]^2)
    cost + reg
  },

  #' @description
  #' Computes the gradient of the log loss function (with regularisation).
  #' Assumes one-hot-encoded categorical variables and bias term in X.
  #' @param theta Coefficients.
  #' @param X matrix with first column as bias.
  #' @param y 0,1 Vector.
  grad_log_loss = function(theta, X, y){
    n <- length(y)
    y_hat <- sigmoid(X %*% theta)
    1/n * (t(X) %*% (y_hat - y)) + c(0, self$lambda * theta[-1])
  },

  #' @description
  #' Fit the object to training data X and y.
  #' @param X Training data (dataframe or matrix).
  #' @param y Training data (vector).
  #' @param ... Additional arguments passed to optim.
  fit = function(X, y, ...){
    X <- parse_X(X)
    self$X <- X
    self$y <- y
    bias <- 1
    X <- cbind(bias, X)
    init_theta <- rep(0, ncol(X))
    optim_res <- optim(init_theta, fn = self$log_loss, gr = self$grad_log_loss,
                   method = self$solver, X = X, y=y,
                   control = list(maxit=10000), ...)
    if (optim_res$convergence == 1){
      warning("The algorithm did not converge.")
    }
    theta <- optim_res$par
    names(theta) <- colnames(X)
    self$theta <- theta
    invisible(self)
  },

  #' @description
  #' Predict on X.
  #' @param X X training or testing X.
  predict = function(X){
    bias <- 1
    X <- cbind(bias, parse_X(X))
    predictions <- as.vector(sigmoid(X %*% self$theta))
    names(predictions) <- rownames(X)
    predictions
  },

  #' @description
  #' Print LogisticRegression object.
  #' @param ... Additional arguments passed to print.
  print = function(...) {
    cat("LogisticRegression: \n")
    cat("  head(X) = ", head(self$X, 1), "\n", sep = " ")
    cat("  head(y) = ", head(self$y, 1), "\n", sep = " ")
    cat("  solver = ", self$solver, "\n", sep = " ")
    cat("  lambda = ", self$lambda, "\n", sep = " ")
    cat("  theta = ", self$theta, "\n", sep = " ")
    invisible(self)
  }
))
