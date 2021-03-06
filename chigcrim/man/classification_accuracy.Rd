% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R
\name{classification_accuracy}
\alias{classification_accuracy}
\title{Accuracy}
\usage{
classification_accuracy(y_hat, y_test)
}
\arguments{
\item{y_hat}{Binary vector of predictions.}

\item{y_test}{Binary vector of observations.}
}
\value{
float Proportion of correct predictions
}
\description{
Computes the proportion of correctly classified predictions.
}
\examples{
classification_accuracy(c(0, 1, 0), c(0, 1, 1))
}
