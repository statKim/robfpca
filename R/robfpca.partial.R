#' Robust functional principal component analysis for partially observed functional data
#'
#' Robust functional principal component analysis (FPCA) for partially observed functional data is performed based on the pairwise robust covariance function estimation and eigenanalysis.
#'
#' The location and scale functions are computed via pointwise M-estimator, and the covariance function is obtained via robust pairwise computation based on Orthogonalized Gnanadesikan-Kettenring (OGK) estimation.
#' Additionally, bivariate Nadaraya-Watson smoothing is applied for smoothed covariance surfaces.
#' To deal with the missing segments, FPCA is performed via PACE (Principal Analysis via Conditional Expectation).
#'
#' @param X  a n x p matrix. It allows NA.
#' @param type  the option for robust dispersion estimator. "huber", "bisquare", and "tdist" are supported.
#' @param grid a vector containing the observed timepoints
#' @param K a number of FPCs. If K is NULL, K is selected by PVE.
#' @param PVE a proportion of variance explained
#' @param MM the option for M-scale estimator in GK identity. If it is FALSE, the same method using \code{type} is used, that is the iterative algorithm. The closed form solution using method of moments can be used when \code{MM == TRUE}. Defalut is TRUE.
#' @param smooth If it is TRUE, bivariate Nadaraya-Watson smoothing is performed using \code{fields::smooth2d()}. Default is TRUE.
#' @param bw a bandwidth when \code{smooth = TRUE}.
#' @param cv If it is TRUE, K-fold cross-validation is performed for the bandwidth selection when \code{smooth = TRUE}.
#' @param df the degrees of freedm when \code{type = "tdist"}.
#' @param cv_optns the options of K-fold cross-validation when \code{cv = TRUE}. See Details.
#'
#' @details The options of \code{cv_optns}:
#' \describe{
#' \item{bw_cand}{a vector contains the candidates of bandwidths for bivariate smoothing.}
#' \item{K}{the number of folds for K-fold cross validation.}
#' \item{ncores}{the number of cores on \code{foreach} for parallel computing.}
#' }
#'
#' @return a list contatining as follows:
#' \item{data}{a matrix which is the input X}
#' \item{lambda}{the first K eigenvalues}
#' \item{eig.fun}{the first K eigenvectors}
#' \item{pc.score}{the first K FPC scores}
#' \item{K}{a number of FPCs}
#' \item{PVE}{a proportion of variance explained}
#' \item{work.grid}{a work grid}
#' \item{eig.obj}{an object of the eigenanalysis}
#' \item{mu}{a mean function}
#' \item{cov}{a covariance function}
#' \item{sig2}{a noise variance}
#' \item{cov.obj}{the object from cov_ogk(). See ?cov_ogk().}
#'
#' @examples
#' set.seed(100)
#' x.list <- sim_delaigle(n = 100,
#'                        type = "partial",
#'                        out.prop = 0.2,
#'                        out.type = 1,
#'                        dist = "normal")
#' x <- list2matrix(x.list)
#'
#' # Given bandwidth
#' fpca.obj <- robfpca.partial(x,
#'                             type = "huber",
#'                             PVE = 0.95,
#'                             bw = 0.1)
#' fpc.score <- fpca.obj$pc.score
#'
#' # Using 5-fold cross-validation
#' bw_cand <- seq(0.02, 0.3, length.out = 10)
#' fpca.obj <- robfpca.partial(x,
#'                             type = "huber",
#'                             PVE = 0.95,
#'                             cv = TRUE,
#'                             cv_optns = list(bw_cand = bw_cand,
#'                                             K = 5,
#'                                             ncores = 1))
#' fpc.score <- fpca.obj$pc.score
#'
#' @references
#' \cite{Park, Y., Kim, H., & Lim, Y. (2022+). Functional principal component analysis for partially observed elliptical process, Under review.}
#'
#' \cite{Maronna, R. A., & Zamar, R. H. (2002). Robust estimates of location and dispersion for high-dimensional datasets. Technometrics, 44(4), 307-317.}
#'
#' \cite{Yao, F., Müller, H. G., & Wang, J. L. (2005). Functional data analysis for sparse longitudinal data. Journal of the American statistical association, 100(470), 577-590.}
#'
#' @import RobStatTM stats
#'
#' @export
robfpca.partial <- function(X,
                            type = c("huber","bisquare","tdist"),
                            grid = NULL,
                            K = NULL,
                            PVE = 0.99,
                            MM = TRUE,
                            smooth = TRUE,
                            bw = NULL,
                            cv = FALSE,
                            df = 3,
                            cv_optns = list(bw_cand = NULL,
                                            K = 5,
                                            ncores = 1)) {

    if (!is.matrix(X)) {
        stop("X should be matrix type!")
    }

    n <- nrow(X)   # number of curves
    p <- ncol(X)   # number of timepoints

    # Obtain work.grid
    if (is.null(grid)) {
        work.grid <- seq(0, 1, length.out = p)
    } else {
        if (length(grid) != p) {
            stop("length(grid) does not equal to ncol(X).")
        }
        work.grid <- grid
    }


    # Estimate robust covariance function using proposed method
    cov.obj <- cov_ogk(X,
                       type = type,
                       MM = MM,
                       smooth = smooth,
                       grid = work.grid,
                       bw = bw,
                       cv = cv,
                       df = df,
                       cv_optns = cv_optns)
    mu.ogk <- cov.obj$mean
    cov.ogk <- cov.obj$cov
    # noise.ogk <- cov.obj$noise.var

    # Transform matrix type to list
    X.list <- matrix2list(X, grid = work.grid)
    Lt <- X.list$Lt
    Ly <- X.list$Ly

    # Functional PCA based on PACE
    fpca.obj <- funPCA(Lt = Lt,
                       Ly = Ly,
                       mu = mu.ogk,
                       cov = cov.ogk,
                       sig2 = 0,
                       work.grid = work.grid,
                       K = K,
                       PVE = PVE)

    fpca.obj$data <- X   # Input data matrix
    fpca.obj$cov.obj <- cov.obj   # Object of covariance estimation

    class(fpca.obj) <- "robfpca.partial"

    return(fpca.obj)
}



#' Print a robfpca.partial object
#'
#' @param x a \code{robfpca.partial} object from \code{robfpca.partial()}
#' @param ... Not used
#'
#' @method print robfpca.partial
#'
#' @export
print.robfpca.partial <- function(x, ...) {
    obj <- x
    cat(paste0("Robust functional principal component analysis for partially observed functional data, ",
               class(obj), "object\n"))
    cat(
        paste0("The number of functional principal components selected is: ", obj$K,
               "and\n its proportion of variance explained is: ", round(obj$PVE, 3))
    )

}



#' Predict FPC scores, reconstruction and completion for a new data
#'
#' @param object a \code{robfpca.partial} object from \code{robfpca.partial()}
#' @param type "score" gives FPC scores, "reconstr" gives reconstruction of each curves, and "comp" gives completion of each curves.
#' @param newdata a n x p matrix containing n curves observed at p timepoints
#' @param K a number of FPCs
#' @param ... Not used
#'
#' @method predict robfpca.partial
#'
#' @export
predict.robfpca.partial <- function(object,
                                    type = c("score","reconstr","comp"),
                                    newdata = NULL,
                                    K = NULL, ...) {
    if (class(object) != "robfpca.partial") {
        stop("Check the class of input object! class name 'robfpca.partial' is only supported!")
    }

    if (is.null(K)) {
        K <- object$K
    }
    if (K > object$K) {
        stop(paste0("Selected number of PCs from robfpca object is less than K."))
    }

    # Prediction of FPC scores
    if (is.null(newdata)) {
        pc.score <- matrix(object$pc.score[, 1:K],
                           ncol = K)
        n <- nrow(pc.score)
        newdata <- object$data
    } else {
        Lt <- newdata$Lt
        Ly <- newdata$Ly
        n <- length(Lt)

        pc.score <- matrix(NA, n, K)
        for (i in 1:n) {
            pc.score[i, ] <- get_CE_score(Lt[[i]],
                                          Ly[[i]],
                                          object$mu,
                                          object$cov,
                                          object$sig2,
                                          object$eig.obj,
                                          K,
                                          object$work.grid)
        }
    }

    if (type == "score") {
        # Predicted FPC score
        pred <- pc.score
    } else {
        # Reconstruction
        mu <- matrix(rep(object$mu, n),
                     nrow = n, byrow = TRUE)
        eig.fun <- matrix(object$eig.fun[, 1:K],
                          ncol = K)
        pred <- mu + pc.score %*% t(eig.fun)   # reconstructed curves

        # Completion
        if (type == "comp") {
            pred[which(is.na(newdata), arr.ind = TRUE)] <- NA
        }
    }

    return(pred)
}


plot.robfpca.partial <- function(object) {

}
