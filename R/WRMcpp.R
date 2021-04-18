#' Weighted repeated median smoothing
#'
#' It is the same algorithm of \code{wrm.smooth()} in \code{robfilter} pacakges.
#' For fast computation, transform the R code to C++.
#' The parameters of a function are the same with "wrm.smooth()" except on "kernel".
#'
#' @param x a list of vectors containing time points for each curve
#' @param y a list of vectors containing observations for each curve
#' @param h bandwidth
#' @param xgrid a vector containing time points to estimate
#' @param kernel "epanechnikov" and "gauss" are supported.
#'
#' @return a smoothed vector
#'
#' @references
#' \cite{Fried, R., Einbeck, J., & Gather, U. (2007). Weighted repeated median smoothing and filtering. Journal of the American Statistical Association, 102(480), 1300-1308.}
#'
#' @export
wrm_smooth <- function(x, y, h, xgrid = NULL, kernel = "epanechnikov") {

  if (is.null(xgrid)) {
    xgrid <- sort(x)
    if (length(xgrid) > 100) {
      xgrid <- seq(min(x), max(x), length.out = 100)
    }
  }

  res <- wrm_smooth_cpp(x, y, h, xgrid, kernel)

  return(res)
}


