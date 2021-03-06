% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{get_missing_distance}
\alias{get_missing_distance}
\title{Get pairwise missingness distances based on discordance}
\usage{
get_missing_distance(missing_matrix)
}
\arguments{
\item{missing_matrix}{matrix of missingness encoded as a binary variable.}
}
\description{
Discordance here is defined as the proportion of loci that do
  not match in their missingness status, i.e. one sample is missing at that
  locus and the other is not. This is equivalent to the Manhattan distance
  between rows of the missingness matrix.
}
