% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utills.R
\name{matrix2list}
\alias{matrix2list}
\title{Convert a matrix to a list}
\usage{
matrix2list(X, grid = NULL)
}
\arguments{
\item{X}{a n x p matrix containing functional trajectories.
n is the number of curves, and p is the number of timepoints.}

\item{grid}{a vector containing observed timepoints.
Default is setting to equal grids between 0 and 1.}
}
\value{
a list containing 2 list (Lt, Ly).
}
\description{
Convert a matrix to a list
}
