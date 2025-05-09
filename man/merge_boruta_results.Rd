% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/machine_learning.R
\name{merge_boruta_results}
\alias{merge_boruta_results}
\title{Merge Boruta Results}
\usage{
merge_boruta_results(
  importance_values,
  decisions,
  file_name,
  iterations,
  threshold,
  return = TRUE
)
}
\arguments{
\item{importance_values}{A list of data frames with feature importance values from each iteration.}

\item{decisions}{A list of character vectors with the decision labels from each iteration.}

\item{file_name}{A string used for naming the output plot file.}

\item{iterations}{Integer. The number of Boruta iterations performed.}

\item{threshold}{A numeric value between 0 and 1. Features labeled as 'Confirmed' or 'Tentative' in more than \code{threshold * iterations} will be retained.}

\item{return}{Logical. Whether to save plots in the "Results/" directory.}
}
\value{
A list containing:
\itemize{
\item A vector of confirmed features.
\item A vector of tentative features.
\item A data frame with median importance values and final decisions.
}
}
\description{
Merge results from multiple Boruta runs to identify robust feature selections.
}
