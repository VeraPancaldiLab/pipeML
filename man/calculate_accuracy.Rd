% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/machine_learning.R
\name{calculate_accuracy}
\alias{calculate_accuracy}
\title{Calculates accuracy values from prediction}
\usage{
calculate_accuracy(metrics, target)
}
\arguments{
\item{metrics}{A data frame with metrics obtained using \code{get_sensitivity_specificity()}, containing
at least two columns: "Sensitivity" and "Specificity".}

\item{target}{A character vector containing the true values from the target variable.
It should have the same length as the predictions.}
}
\value{
A numeric vector representing the accuracy values.
The result is the fraction of correct predictions out of all predictions.
}
\description{
This function calculates the accuracy of a model based on the provided metrics and the true target values.
The accuracy is computed as the ratio of correct predictions (both true positives and true negatives)
to the total number of predictions.
}
