% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCA.R
\name{predict.funPCA}
\alias{predict.funPCA}
\title{Reconstruction via functional PCA
newdata should be a list containing Lt and Ly.}
\usage{
\method{predict}{funPCA}(object, newdata = NULL, K = NULL, ...)
}
\arguments{
\item{object}{a \code{funPCA} object from \code{funPCA()}}

\item{newdata}{a list containing \code{Lt} and \code{Ly}}

\item{K}{a number of PCs for reconstruction}

\item{...}{Not used}
}
\description{
Reconstruction via functional PCA
newdata should be a list containing Lt and Ly.
}
