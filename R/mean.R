########################################################
### Robust mean estimation for functional snippets
########################################################

### Local polynomial kernel smoothing with huber loss (mean estimator)
# Lt : a list of vectors or a vector containing time points for all curves
# Ly : a list of vectors or a vector containing observations for all curves
# newt : a vector containing time points to estimate
# method : "Huber"
# bw : bandwidth
# kernel : a kernel function for kernel smoothing ("epanechnikov", "gauss" are supported.)
# deg : a numeric scalar of polynomial degrees for the kernel smoother
# k2 : If method == "Huber", it uses for \rho function in Huber loss.
# cv : If cv == TRUE, it performs K-fold cross-validation for optimal bandwidth(bw).
# cv_optns : a option for K-fold cross-validation if cv == TRUE.
# loss : a loss function for kernel smoothing("L2" is squared loss, "Huber" is huber loss.)
#   For loss = "Huber", it uses `rlm()` in `MASS` and fits the robust regression with Huber loss.
#   So additional parameters of `rlm()` can be applied. (k2, maxit, ...)

#' Local polynomial kernel smoothing with huber loss (mean estimator)
#'
#' @param Lt a list of vectors containing time points for each curve
#' @param Ly a list of vectors containing observations for each curve
#' @param newt a vector containing time points to estimate
#' @param method "huber", "WRM" are supported
#' @param bw a bandwidth
#' @param kernel a kernel function for kernel smoothing ("epanechnikov", "gauss" are supported.)
#' @param deg a numeric scalar of polynomial degrees for the kernel smoother
#' @param k2 If method == "Huber", it uses for $\\rho$ function in Huber loss.
#' @param ncores If ncores > 1, it implements \code{foreach()} in \code{doParallel} for CV.
#' @param cv_delta_loss a loss function for K-fold cross-validation for delta in Huber function
#' @param cv_bw_loss a loss function for K-fold cross-validation for bandwidth
#' @param cv_K a number of folds for K-fold cross-validation
#' @param ... additional options
#'
#' @import dplyr
#'
#' @export
#' @useDynLib robfpca
meanfunc.rob <- function(Lt,
                         Ly,
                         newt = NULL,
                         method = c("L2","Huber","WRM","Bisquare"),
                         bw = NULL,
                         kernel = "epanechnikov",
                         deg = 1,
                         k2 = 1.345,
                         ncores = 1,
                         cv_delta_loss = "L1",
                         cv_bw_loss = "HUBER",
                         cv_K = 5,
                         # cv = FALSE,
                         # cv_optns = list(K = 5,
                         #                 ncores = 1,
                         #                 Loss = "Huber"),
                         ...) {
    method <- toupper(method)
    if (!(method %in% c("L2","HUBER","WRM","BISQUARE"))) {
        stop(paste0(method, " is not provided. Check method parameter."))
    }

    R <- NULL

    others <- list(...)
    domain <- c(0,0)
    domain[1] <- min(unlist(Lt))
    domain[2] <- max(unlist(Lt))
    domain <- c(domain[1]-0.01*(domain[2]-domain[1]),
                domain[2]+0.01*(domain[2]-domain[1]))

    if (!is.list(Lt)) {
        stop("Lt and Ly must be list type.")
    }

    # kernel <- tolower(get.optional.param('kernel',others,'epanechnikov'))
    # bw <- get.optional.param('bw',others,NULL)

    cv <- FALSE
    # 5-fold CV for delta in Huber function
    if ((method %in% c("HUBER","BISQUARE")) && (is.null(k2) | ncores > 1)) {
        print(paste0(cv_K, "-fold CV is performed for delta in Huber function."))
        delta_cv_obj <- delta.local_kern_smooth(Lt = Lt,
                                                Ly = Ly,
                                                method = method,
                                                kernel = kernel,
                                                deg = deg,
                                                # k2 = k2,
                                                bw = bw,
                                                cv_loss = cv_delta_loss,
                                                K = cv_K,
                                                ncores = ncores,
                                                ...)
        k2 <- delta_cv_obj$selected_delta
        cv <- TRUE
    }

    # 5-fold CV for bandwidth
    if (is.null(bw) | ncores > 1) {
        print(paste0(cv_K, "-fold CV is performed for bandwidth."))
        bw_cv_obj <- bw.local_kern_smooth(Lt = Lt,
                                          Ly = Ly,
                                          method = method,
                                          kernel = kernel,
                                          deg = deg,
                                          k2 = k2,
                                          cv_loss = cv_bw_loss,
                                          K = cv_K,
                                          ncores = ncores,
                                          ...)
        bw <- bw_cv_obj$selected_bw
        cv <- TRUE
    }

    n <- length(Lt)
    t <- unlist(Lt)
    y <- unlist(Ly)

    ord <- sort(t, index.return = T)$ix
    t <- t[ord]
    y <- y[ord]

    R <- list(bw = bw,
              t = t,
              y = y,
              n = n,
              domain = domain,
              method = method,
              kernel = kernel,
              deg = deg,
              k2 = k2,
              yend = c(NULL,NULL))
    class(R) <- 'meanfunc.rob'

    L0 <- domain[2]-domain[1]
    yend <- predict(R, c(domain[1]+L0/100, domain[2]-L0/100))
    R$yend <- yend

    if (!is.null(newt)) {
        R$fitted <- predict(R, newt)
    }

    if (cv == TRUE) {
        # R$cv_optns <- cv_optns
        R$cv_optns <- list(K = cv_K,
                           ncores = ncores,
                           delta_loss = cv_delta_loss,
                           bw_loss = cv_bw_loss)
    }

    return(R)
}




### Predict mean at new time points
#' Predict mean at new time points
#'
#' @param object an object from \code{meanfunc.rob()}
#' @param newt new time points to predict
#' @param ... not used
#'
#' @importFrom stats predict
#'
#' @export
predict.meanfunc.rob <- function(object, newt, ...) {
    meanfunc.obj <- object
    pred <- function(newt) {   # newt must be a vector
        idxl <- newt < meanfunc.obj$domain[1]
        idxu <- newt > meanfunc.obj$domain[2]
        idx <- (!idxl) & (!idxu)

        newt0 <- newt[idx]
        ord <- sort(newt0, index.return = T)$ix

        tmp <- rep(Inf, length(newt0))

        # print(meanfunc.obj$method)

        tmp[ord] <- local_kern_smooth(Lt = meanfunc.obj$t,
                                      Ly = meanfunc.obj$y,
                                      newt = newt0[ord],
                                      method = meanfunc.obj$method,
                                      bw = meanfunc.obj$bw,
                                      kernel = meanfunc.obj$kernel,
                                      # loss = "Huber",
                                      k2 = meanfunc.obj$k2)

        yhat <- rep(0, length(newt))
        yhat[idx] <- tmp
        yhat[idxl] <- meanfunc.obj$yend[1]
        yhat[idxu] <- meanfunc.obj$yend[2]

        return(yhat)
    }

    if (is.list(newt)) {   # Need to test this code for list type.
        mi <- lapply(newt, length)
        newt <- unlist(newt)
        fitted <- pred(newt)

        cm <- c(0, cumsum(mi))

        R <- sapply(1:length(mi), function(i) {
            res <- list()
            res[[1]] <- fitted[(cm[i]+1):cm[i+1]]
            res
        })

        return(R)
    } else if (is.vector(newt)) {
        # obtain unique time point => merge to newt
        df_newt <- data.frame(t = newt)
        newt_unique <- unique(newt)
        pred_unique <- pred(newt_unique)
        df_newt_unique <- data.frame(t = newt_unique,
                                     pred = pred_unique)
        pred_newt <- df_newt %>%
            left_join(df_newt_unique, by = "t") %>%
            dplyr::select(pred) %>%
            unlist()

        return(as.numeric(pred_newt))
        # return(pred(newt))
    }
    else stop('newt must be a vector or a list of vectors of real numbers')
}
