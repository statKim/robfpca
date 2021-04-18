### Weighted repeated median smoothing
# It is the same algorithm of "wrm.smooth()" in robfilter pacakges.
# For fast computation, transform the R code to C++.
# The parameters of a function are the same with "wrm.smooth()" except on "kernel".
# "kernel" option is provided "gauss" and "epanechnikov" only.
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


