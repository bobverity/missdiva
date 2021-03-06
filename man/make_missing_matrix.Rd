% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{make_missing_matrix}
\alias{make_missing_matrix}
\title{Make a matrix of missingness using coverage matrix}
\usage{
make_missing_matrix(coverage = NULL, threshold = NULL)
}
\arguments{
\item{coverage}{coverage matrix, should be numeric values, can contain NAs.}

\item{threshold}{coverage values below this threshold are considered missing
data.}
}
\description{
Use coverage values extracted by the function
  \code{get_coverage_matrix()} and, using a variable threshold, make a binary
  matrix of missingness.
}
