% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_corr.R
\name{sim_corr}
\alias{sim_corr}
\title{Generate partially observed functional data with spatially correlated case}
\usage{
sim_corr(
  n = 403,
  type = c("partial", "snippet", "dense"),
  out.prop = 0.2,
  out.type = 1,
  dist = "normal",
  noise = 0,
  dist.mat = dist.mat,
  d = 1.4,
  f = 0.2,
  r.par = 200
)
}
\arguments{
\item{n}{a number of curves}

\item{type}{the type of generated data. "partial" means the option for partially observed data, "snippet" is short fragmented data, and "dense" means the fully observed curves.}

\item{out.prop}{a proportion of outlying curves of total n curves. Only used for dist = "normal".}

\item{out.type}{a outlier type, 1~3 are supported. Only used for dist = "normal".}

\item{dist}{a distribution which the data is generated. "normal"(Normal distribution) and "tdist"(t-distribution) are supported. If dist = "tdist", the option of \code{out.prop} and \code{out.type} are ignored.}

\item{noise}{a numeric value which is added random gaussian noises. Default is 0(No random noise).}

\item{dist.mat}{a 403 spatial locations and their distances matrix}

\item{d}{a parameter for missingness when \code{type} is "partial" (See Kraus(2015))}

\item{f}{a parameter for missingness when \code{type} is "partial" (See Kraus(2015))}

\item{r.par}{a parameter that control the spatial correlation}
}
\value{
a list contatining as follows:
\item{Ly}{a list of n vectors containing the observed values for each individual.}
\item{Lt}{a list of n vectors containing the observation time points for each individual corresponding to \code{Ly}}
\item{out.ind}{a vector containing outlier index. 0 is non-outlier and 1 is the outlier.}
\item{x.full}{a n x 51 dense matrix with n observations per 51 timepoints before making partially observed.}
}
\description{
Partially observed functional data is generated with 51 regular grids from the spatially correlated model
}
\examples{
set.seed(100)
n <- 100
x.list <- sim_corr(n = n,
                   type = "partial",
                   out.prop = 0.2,
                   dist = "normal",
                   dist.mat = dist.mat[1:n, 1:n])
x <- list2matrix(x.list)
matplot(t(x), type = "l")

}
\references{
\cite{Kraus, D. (2015). Components and completion of partially observed functional data. Journal of the Royal Statistical Society: Series B: Statistical Methodology, 777-801.}
}
