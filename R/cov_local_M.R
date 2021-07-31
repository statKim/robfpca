#################################################################
### Local M-estimator for partially observed functional data
#################################################################


get_kernel <- function(t_vec, t, bw, method = "epan") {
    u <- (t_vec - t) / bw
    kern <- (1-u^2) * (3/4)

    return(kern)
}


#' Local M-estimator for Mean function (mean estimator)
#'
#' Obtain the local M-estimator for mean function.
#'
#' @param X a n x p matrix
#' @param bw a numeric value of bandwidth. (Defalut is (max(gr) - min(gr)) / 5.)
#' @param gr a vector containing observed timepoints. (Default is \code{NULL}, and it is set to equal grid between 0 and 1.)
#' @param cv If it is \code{TRUE}, robust 5-fold cross-validation(using Huber loss) is performed to select the optimal bandwidth. (Default is \code{FALSE}.)
#' @param ncores a number of cores to select bandwidth from robust 5-fold cross-validation. (Defalut is 1.)
#'
#' @return a vector containing the mean function from local M-estimation.
#'
#' @export
#' @useDynLib robfpca
mean_local_M <- function(X,
                         bw = NULL,
                         gr = NULL,
                         cv = FALSE,
                         ncores = 1) {
    n <- nrow(X)
    p <- ncol(X)

    if (is.null(gr)) {
        gr <- seq(0, 1, length.out = p)
    }

    if (is.null(bw)) {
        bw <- (max(gr) - min(gr)) / 5
        cat(paste("bw is not specified. ", round(bw, 3), "is used. \n"))
    }

    if (cv == TRUE) {
        cv.obj <- cv.mean_local_M(X, ncores = ncores)   # 5-fold CV is performed
        bw <- cv.obj$selected_bw
        cat(paste("Optimal bandwidth=", round(bw, 3), "is selected! \n"))
    }

    # transform to observed value per timepoint
    x <- as.numeric(X)
    t <- rep(gr, each = n)

    # remove NA
    ind_NA <- which(is.na(x))
    x <- x[-ind_NA]
    t <- t[-ind_NA]

    # call C++ function
    mu <- local_M_1D(x = x,
                     t = t,
                     new_t = gr,
                     h = bw)

    return(mu)
}


cv.mean_local_M <- function(X,
                            ncores = 1) {

}



#' Get Raw Covariances for Functional Data
#'
#' Obtain the raw covariances $C(s,t)$.
#'
#' It provides 2 options for implementing R or C++.
#'
#' @param X a n x p matrix
#' @param mu a vector containing the mean function. (Default is \code{NULL}, and it is calculated automatically.)
#' @param gr a vector containing observed timepoints. (Default is \code{NULL}, and it is set to equal grid between 0 and 1.)
#' @param diag a logical value whether containing diagonal parts of the covariance. (Default is \code{FALSE}.)
#' @param engine the option for implementation. (Defalut is \code{C++}, and it is much faster than \code{R}.)
#'
#' @return a dataframe contatining as follows:
#' \item{cov}{the raw covariance $C(s,t)$ between time \code{s} and \code{t} for \code{ind}th curve}
#' \item{s}{time 1}
#' \item{t}{time 2}
#' \item{ind}{the index of each curve}
#'
#' @export
#' @useDynLib robfpca
get_raw_cov <- function(X,
                        mu = NULL,
                        gr = NULL,
                        diag = FALSE,
                        engine = "C++") {
    if (is.null(gr)) {
        gr <- seq(0, 1, length.out = ncol(X))
    }
    n <- nrow(X)
    p <- ncol(X)

    if (is.null(mu)) {
        mu <- mean_Mest(X)
    }

    if (engine == "C++") {
        raw_cov <- get_raw_cov_cpp(X, mu, gr, diag)
        raw_cov <- as.data.frame(raw_cov)
    } else if (engine == "R") {
        for (i in 1:n) {
            ind <- which(!is.na(X[i, ]))
            gr_sub <- gr[ind]
            x_sub <- X[i, ind] - mu[ind]

            exp_grid <- expand.grid(gr_sub, gr_sub)
            exp_x <- expand.grid(x_sub, x_sub)

            if (diag == FALSE) {
                # eliminate diagonal parts
                idx <- which(exp_grid[, 1] == exp_grid[, 2])
                cov_i <- exp_x[-idx, 1] * exp_x[-idx, 2]

                if (i == 1) {
                    raw_cov <- cbind(cov_i,
                                     exp_grid[-idx, ],
                                     i)
                } else {
                    raw_cov <- rbind(raw_cov,
                                     cbind(cov_i,
                                           exp_grid[-idx, ],
                                           i))
                }
            } else {
                cov_i <- exp_x[, 1] * exp_x[, 2]

                if (i == 1) {
                    raw_cov <- cbind(cov_i,
                                     exp_grid,
                                     i)
                } else {
                    raw_cov <- rbind(raw_cov,
                                     cbind(cov_i,
                                           exp_grid,
                                           i))
                }
            }
        }
    }

    colnames(raw_cov) <- c("cov","s","t","ind")

    return(raw_cov)
}



#' Bivariate Local M-estimator for Covariance function (covariance estimator)
#'
#' Obtain the bivariate local M-estimator for covariance function.
#'
#' It provides 2 options for implementing R or C++.
#'
#' @param X a n x p matrix
#' @param bw a numeric value of bandwidth. (Defalut is (max(gr) - min(gr)) / 5.)
#' @param gr a vector containing observed timepoints. (Default is \code{NULL}, and it is set to equal grid between 0 and 1.)
#' @param diag logical value whether containing diagonal parts of the covariance. (Default is \code{FALSE}.)
#' @param cv If it is \code{TRUE}, robust 5-fold cross-validation(using Huber loss) is performed to select the optimal bandwidth. (Default is \code{FALSE}.)
#' @param ncores a number of cores to select bandwidth from robust 5-fold cross-validation. (Defalut is 1.)
#' @param raw_cov a dataframe which is the return value from \code{get_raw_cov()}. (Default is \code{NULL}, and it is calculated automatically.)
#' @param engine the option for implementation. (Defalut is \code{C++}, and it is much faster than \code{R}.)
#'
#' @return a p x p matrix containing the covariance function from local M-estimation.
#'
#' @export
#' @useDynLib robfpca
cov_local_M <- function(X,
                        bw = NULL,
                        gr = NULL,
                        diag = FALSE,
                        cv = FALSE,
                        ncores = 1,
                        raw_cov = NULL,
                        engine = "C++") {
    if (is.null(gr)) {
        gr <- seq(0, 1, length.out = ncol(X))
    }
    p <- length(gr)


    if (is.null(bw)) {
        bw <- (max(gr) - min(gr)) / 5
        cat(paste("bw is not specified. ", round(bw, 3), "is used. \n"))
    }

    if (cv == TRUE) {
        cv.obj <- cv.cov_local_M(X, ncores = ncores)   # 5-fold CV is performed
        bw <- cv.obj$selected_bw
        cat(paste("Optimal bandwidth=", round(bw, 3), "is selected! \n"))
    }

    # raw covariance object from "get_raw_cov" function
    if (is.null(raw_cov)) {
        raw_cov <- get_raw_cov(X,
                               mu = NULL,
                               gr = gr,
                               diag = diag,
                               engine = engine)
    }
    s <- raw_cov$s
    t <- raw_cov$t
    raw_cov <- raw_cov$cov

    if (engine == "C++") {
        cov_mat <- local_M_2D(raw_cov = raw_cov,
                              s = s,
                              t = t,
                              gr = gr,
                              h = bw)
    } else if (engine == "R") {
        cov_mat <- matrix(NA, p, p)
        for (i in 1:p) {
            i_neighbor <- which(s < gr[i] + bw & s > gr[i] - bw)
            # kern_s <- get_kernel(s, gr[i], h)
            # i_neighbor <- which(kern_s > 0)
            for (j in 1:p) {
                if (i <= j) {
                    j_neighbor <- which(t < gr[j] + bw & t > gr[j] - bw)

                    ind <- intersect(i_neighbor, j_neighbor)

                    cov_mat[i, j] <- MASS::huber(raw_cov[ind])$mu
                } else {
                    cov_mat[i, j] <- cov_mat[j, i]
                }

                # kern_t <- get_kernel(t, gr[j], h)
                # j_neighbor <- which(kern_t > 0)
                #
                # ind <- intersect(i_neighbor, j_neighbor)
                # # kern_weig <- kern_s[ind] * kern_t[ind]
                # # kern_weig <- kern_weig / max(kern_weig)
                # #
                # # cov_weig <- kern_weig * raw_cov$cov[ind]
                # cov_weig <- raw_cov$cov[ind]
                #
                # cov_mat[i, j] <- MASS::huber(cov_weig)$mu
            }
        }
    }

    # make positive definite
    eig <- get_eigen(cov_mat, gr)
    lambda <- eig$lambda
    idx_positive <- which(lambda > 0)
    lambda <- lambda[idx_positive]
    phi <- eig$phi[, idx_positive]
    cov_mat <- phi %*% diag(lambda, ncol = length(idx_positive)) %*% t(phi)

    return(cov_mat)
}



#' Robust K-fold Cross-Validation for the Optimal Bandwidth of Bivariate Local M-estimator for Covariance function
#'
#' Perform robust K-fold cross-validation to select the optimal bandwidth of bivariate local M-estimator for covariance function.
#'
#' It provides 2 options for implementing R or C++.
#'
#' @param X a n x p matrix
#' @param bw_cand a numeric vector of bandwidth candidate. (Default is \code{NULL}, and it is selected automatically.)
#' @param K a number of folds for K-fold cross-validation. (Default is 5.)
#' @param ncores a number of cores to select bandwidth from robust 5-fold cross-validation. (Defalut is 1.)
#' @param delta a cut-off value in Huber loss function. (Default is 1.345.)
#' @param engine the option for implementation. (Defalut is \code{C++}, and it is much faster than \code{R}.)
#'
#' @return a list contatining as follows:
#' \item{selected_bw}{the optimal bandwidth selected from the robust K-fold cross-validation.}
#' \item{cv.error}{a matrix containing CV error per bandwidth candidates.}
#'
#' @export
cv.cov_local_M <- function(X,
                           bw_cand = NULL,
                           K = 5,
                           ncores = 1,
                           delta = 1.345,
                           engine = "C++") {

    if (is.list(X)) {
        gr <- sort(unique(unlist(X$Lt)))
        X <- list2matrix(X)
    } else {
        gr <- seq(0, 1, length.out = ncol(X))
    }

    n <- nrow(X)
    p <- ncol(X)

    # raw covariance
    mu <- mean_Mest(X)
    cov_raw <- get_raw_cov(X,
                           mu = mu,
                           gr = gr,
                           diag = FALSE,
                           engine = engine)

    gr_expand <- expand.grid(gr, gr)   # dataframe type
    colnames(gr_expand) <- c("s","t")

    # bandwidth candidates
    if (is.null(bw_cand)) {
        a <- min(gr)
        b <- max(gr)
        bw_cand <- seq(min(diff(gr))*1.5, (b-a)/3, length.out = 10)
        # bw_cand <- 10^seq(-2, 0, length.out = 10) * (b - a)/3
    }

    # get index for each folds
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

    # K-fold cross-validation
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
                                     .packages = c("robfpca"),
                                     .errorhandling = "pass") %dopar% {
        bw <- bw_fold_mat$bw_cand[i]   # bandwidth candidate
        k <- bw_fold_mat$fold[i]   # fold for K-fold CV

        # data of kth fold
        X_train <- X[-folds[[k]], ]
        X_test <- X[folds[[k]], ]
        cov_train <- cov_raw[-which(cov_raw$ind %in% folds[[k]]), ]   # dataframe type
        cov_test <- cov_raw[which(cov_raw$ind %in% folds[[k]]), ]   # dataframe type

        # smoothed covariance without kth fold
        cov_sm <- cov_local_M(X_train,
                              bw = bw,
                              gr = gr,
                              raw_cov = cov_train,
                              diag = FALSE,
                              engine = engine)

        cov_sm <- cbind(gr_expand,
                        cov_sm = as.numeric(cov_sm))   # dataframe type

        df <- dplyr::left_join(cov_test,
                               cov_sm[cov_sm$s != cov_sm$t, ],   # remove diagonal parts
                               by = c("s","t"))

        # robust loss (Huber loss)
        a <- abs(df$cov - df$cov_sm)
        err_huber <- ifelse(a > delta, delta*(a - delta/2), a^2/2)
        err <- sum(err_huber)/K

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
            X_train <- X[-folds[[k]], ]
            X_test <- X[folds[[k]], ]
            cov_train <- cov_raw[-which(cov_raw$ind %in% folds[[k]]), ]   # dataframe type
            cov_test <- cov_raw[which(cov_raw$ind %in% folds[[k]]), ]   # dataframe type

            for (i in 1:length(bw_cand)) {
                # smoothed covariance without kth fold
                cov_sm <- cov_local_M(X_train,
                                      bw = bw_cand[i],
                                      gr = gr,
                                      raw_cov = cov_train,
                                      diag = FALSE,
                                      engine = engine)

                cov_sm <- cbind(gr_expand,
                                cov_sm = as.numeric(cov_sm))   # dataframe type

                df <- dplyr::left_join(cov_test,
                                       cov_sm[cov_sm$s != cov_sm$t, ],   # remove diagonal parts
                                       by = c("s","t"))

                # robust loss (Huber loss)
                a <- abs(df$cov - df$cov_sm)
                err_huber <- ifelse(a > delta, delta*(a - delta/2), a^2/2)
                cv_error[i] <- cv_error[i] + sum(err_huber)/K
            }
        }

        bw <- list(selected_bw = bw_cand[ which.min(cv_error) ],
                   cv.error = data.frame(bw = bw_cand,
                                         error = cv_error))
    }

    return(bw)
}


