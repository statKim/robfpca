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
#' @param smooth If it is \code{TRUE}, smoothing spline is performed in \code{stats::smooth.spline()}. (Default is \code{FALSE}.)
#'
#' @return a p-dimensional mean function
#'
#' @export
#' @useDynLib robfpca
mean_Mest <- function(x, smooth = FALSE) {
    if (is.list(x)) {
        x <- list2matrix(x)
    }
    mu <- mean_Mest_cpp(x, smooth)

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
                     smooth = FALSE,
                     bw = NULL,
                     make.pos.semidef = TRUE,
                     noise.var = 0) {
    if (is.list(x)) {
        x <- list2matrix(x)
    }

    n <- nrow(x)
    p <- ncol(x)

    # obtain the covariance based on marignal M-estimator via C++ code
    rob.var <- cov_Mest_cpp(x, smooth)

    # 2-dimensional smoothing - does not need to adjust noise variance
    if (smooth == T) {
        gr <- seq(0, 1, length.out = p)
        cov.sm.obj <- refund::fbps(rob.var, list(x = gr,
                                                 z = gr))
        rob.var <- cov.sm.obj$Yhat
        # if (is.null(bw)) {
        #     cv.obj <- cv.cov_Mest(x, ncores = 1)   # 5-fold CV is performed
        #     bw <- cv.obj$selected_bw
        #     cat(paste("Optimal bandwidth=", round(bw, 3), "is selected! \n"))
        # }
        #
        # # element-wise covariances
        # gr <- seq(0, 1, length.out = p)
        # st <- expand.grid(gr, gr)
        # cov_st <- as.numeric(rob.var)
        #
        # # # remove diagonal parts from raw covariance (See Yao et al.(2005))
        # # ind <- which(st[, 1] == st[, 2], arr.ind = T)
        # # st <- st[-ind, ]
        # # cov_st <- cov_st[-ind]
        #
        # rob.var <- fields::smooth.2d(cov_st,
        #                              x = st,
        #                              surface = F,
        #                              theta = bw,
        #                              nrow = p,
        #                              ncol = p)
    }
    # else {
    #     # subtract noise variance - Need for not smoothing
    #     diag(rob.var) <- diag(rob.var) - noise.var
    # }

    # make positive-semi-definite
    if (isTRUE(make.pos.semidef)) {
        eig <- eigen(rob.var)
        k <- which(eig$values > 0)
        lambda <- eig$values[k]
        phi <- matrix(eig$vectors[, k],
                      ncol = length(k))

        rob.var <- phi %*% diag(lambda, ncol = length(k)) %*% t(phi)

        # if (length(k) > 1) {
        #     rob.var <- eig$vectors[, k] %*% diag(eig$values[k]) %*% t(eig$vectors[, k])
        # } else {
        #     rob.var <- eig$values[k] * (eig$vectors[, k] %*% t(eig$vectors[, k]))
        # }
    }

    # subtract noise variance
    diag(rob.var) <- diag(rob.var) - noise.var

    return(rob.var)
}



#' K-fold Cross-Validation for Bivariate Nadaraya-Watson smoothing for Covariance
#'
#' Perform K-fold cross-validation for Bivariate Nadaraya-Watson smoothing for Covariance.
#'
#' @param x a n x p matrix
#' @param bw_cand a numeric vector of bandwidth candidate. (Default is \code{NULL}, and it is selected automatically.)
#' @param K a number of folds for K-fold cross-validation. (Default is 5.)
#' @param ncores a number of cores to select bandwidth from robust 5-fold cross-validation. (Defalut is 1.)
#' @param noise.var the numeric value containing noise variance. (Defalut is 0.)
#'
#' @return a list contatining as follows:
#' \item{selected_bw}{the optimal bandwidth selected from the robust K-fold cross-validation.}
#' \item{cv.error}{a matrix containing CV error per bandwidth candidates.}
#'
#' @export
### - Not exactly observation-wise cross-validation
### - It is conducted for element-wise covariance
cv.cov_Mest <- function(x,
                        bw_cand = NULL,
                        K = 5,
                        ncores = 1,
                        noise.var = 0) {

    if (is.list(x)) {
        gr <- sort(unique(unlist(x$Lt)))
        x <- list2matrix(x)
    } else {
        gr <- seq(0, 1, length.out = ncol(x))
    }

    # n <- nrow(x)
    p <- ncol(x)

    # bandwidth candidates
    if (is.null(bw_cand)) {
        a <- min(gr)
        b <- max(gr)
        bw_cand <- 10^seq(-2, 0, length.out = 10) * (b - a)/3
    }

    # obtain the raw covariance (Not smoothed)
    cov_hat <- cov_Mest(x,
                        smooth = FALSE,
                        make.pos.semidef = FALSE)

    # subtract noise variance
    diag(cov_hat) <- diag(cov_hat) - noise.var

    # element-wise covariances
    st <- expand.grid(gr, gr)
    cov_st <- as.numeric(cov_hat)

    # # remove diagonal parts from raw covariance (See Yao et al.(2005))
    # ind <- which(st[, 1] == st[, 2], arr.ind = T)
    # st <- st[-ind, ]
    # cov_st <- cov_st[-ind]

    # get index for each folds
    n <- length(st)
    folds <- list()
    fold_num <- n %/% K   # the number of curves for each folds
    fold_sort <- sample(1:n, n)
    for (k in 1:K) {
        ind <- (fold_num*(k-1)+1):(fold_num*k)
        if (k == K) {
            ind <- (fold_num*(k-1)+1):n
        }
        folds[[k]] <- fold_sort[ind]
    }

    # K-fold cross validation
    if (ncores > 1) {
        # Parallel computing setting
        if (ncores > parallel::detectCores()) {
            ncores <- parallel::detectCores() - 1
            warning(paste0("ncores is too large. We now use ", ncores, " cores."))
        }
        # ncores <- detectCores() - 3
        cl <- parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)

        # matrix of bw_cand and fold
        bw_fold_mat <- data.frame(bw_cand = rep(bw_cand, each = K),
                                  fold = rep(1:K, length(bw_cand)))

        cv_error <- foreach::foreach(i = 1:nrow(bw_fold_mat),
                                     .combine = "c",
                                     # .export = c("local_kern_smooth"),
                                     .packages = c("robfpca"),
                                     .errorhandling = "pass") %dopar% {
            bw <- bw_fold_mat$bw_cand[i]   # bandwidth candidate
            k <- bw_fold_mat$fold[i]   # fold for K-fold CV

            # data of kth fold
            st_train <- st[-folds[[k]], ]
            st_test <- st[folds[[k]], ]
            cov_train <- cov_st[-folds[[k]]]
            cov_test <- cov_st[folds[[k]]]

            # Bivariate Nadaraya-Watson smoothing
            cov_hat_sm <- fields::smooth.2d(cov_train,
                                            x = st_train,
                                            surface = F,
                                            theta = bw,
                                            nrow = p,
                                            ncol = p)
            cov_hat_sm <- as.numeric(cov_hat_sm)[folds[[k]]]
            err <- sum((cov_test - cov_hat_sm)^2)   # squared errors

            return(err)
        }
        parallel::stopCluster(cl)

        bw_fold_mat$cv_error <- cv_error
        cv_obj <- bw_fold_mat %>%
            dplyr::group_by(bw_cand) %>%
            dplyr::summarise(cv_error = sum(cv_error))

        bw <- list(selected_bw = cv_obj$bw_cand[ which.min(cv_obj$cv_error) ],
                   cv.error = as.data.frame(cv_obj))
    } else {
        cv_error <- rep(0, length(bw_cand))
        for (k in 1:K) {
            # data of kth fold
            st_train <- st[-folds[[k]], ]
            st_test <- st[folds[[k]], ]
            cov_train <- cov_st[-folds[[k]]]
            cov_test <- cov_st[folds[[k]]]

            for (i in 1:length(bw_cand)) {
                # Bivariate Nadaraya-Watson smoothing
                cov_hat_sm <- fields::smooth.2d(cov_train,
                                                x = st_train,
                                                surface = F,
                                                theta = bw_cand[i],
                                                nrow = p,
                                                ncol = p)
                cov_hat_sm <- as.numeric(cov_hat_sm)[folds[[k]]]
                cv_error[i] <- cv_error[i] + sum((cov_test - cov_hat_sm)^2)   # squared errors
            }
        }

        bw <- list(selected_bw = bw_cand[ which.min(cv_error) ],
                   cv.error = data.frame(bw = bw_cand,
                                         error = cv_error))
    }

    return(bw)
}
