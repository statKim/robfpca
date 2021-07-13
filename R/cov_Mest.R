#################################################################
### Marginal M-estimator for partially observed functional data
#################################################################

#' Marignal M-estimator for Mean function (mean estimator)
#'
#' Obtain the marginal M-estimator for mean function.
#'
#' It is performed by C++ code for fast computation.
#'
#' @param x a n x p matrix or a list containing Lt and Ly
#'
#' @return a p-dimensional mean function
#'
#' @export
#' @useDynLib robfpca
mean_Mest <- function(x) {
    if (is.list(x)) {
        x <- list2matrix(x)
    }
    mu <- mean_Mest_cpp(x)

    return(mu)
}


#' Marginal M-estimator for Covaraince function (covariance estimator)
#'
#' Obtain the marginal M-estimator for covariance function.
#'
#' @param x a n x p matrix or a list containing Lt and Ly
#' @param smooth If it is \code{TRUE}, the bivariate Nadaraya-Watson smoothing is performed for the covariance estimator. (Default is \code{FALSE}.)
#' @param bw an optional parameter when \code{smooth = TRUE} for a bandwidth of the bivariate smoothing.
#' @param make.pos.semidef a logical value whether made a covarariance estimator positive semi-definite (Default is \code{TRUE}.)
#' @param noise.var a numeric value of the noise variance
#'
#' @return a p x p covariance function
#'
#' @export
#' @useDynLib robfpca
### Marginal M-estimator for covaraince function
## It performed by C++ code.
cov_Mest <- function(x,
                     smooth = F,
                     bw = 0.1,
                     make.pos.semidef = TRUE,
                     noise.var = 0) {
    if (is.list(x)) {
        x <- list2matrix(x)
    }

    n <- nrow(x)
    p <- ncol(x)

    # obtain the covariance based on marignal M-estimator via C++ code
    rob.var = cov_Mest_cpp(x)

    # subtract noise variance
    diag(rob.var) <- diag(rob.var) - noise.var

    # 2-dimensional smoothing
    if (smooth == T) {
        gr <- seq(0, 1, length.out = p)
        rob.var <- fields::smooth.2d(as.numeric(rob.var),
                                     x = expand.grid(gr, gr), surface = F,
                                     theta = bw, nrow = p, ncol = p)
    }

    # make positive-semi-definite
    if (isTRUE(make.pos.semidef)) {
        eig <- eigen(rob.var)
        k <- which(eig$values > 0)
        if (length(k) > 1) {
            rob.var <- eig$vectors[, k] %*% diag(eig$values[k]) %*% t(eig$vectors[, k])
        } else {
            rob.var <- eig$values[k] * (eig$vectors[, k] %*% t(eig$vectors[, k]))
        }
    }

    return(rob.var)
}

