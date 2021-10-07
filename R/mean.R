########################################################
### Robust mean estimation for functional snippets
########################################################

#' Local polynomial kernel smoothing with huber loss (mean estimator)
#'
#' @param Lt a list of vectors containing time points for each curve
#' @param Ly a list of vectors containing observations for each curve
#' @param newt a vector containing time points to estimate
#' @param method "huber", "WRM" are supported
#' @param kernel a kernel function for kernel smoothing ("epanechnikov", "gauss" are supported.)
#' @param bw a bandwidth.
#' @param delta cut-off value for "huber"(Huber) or "bisquare"(Tukey's biweight function).
#' Default is 1.345 for "huber" and 4.685 for "bisquare" for 95\% ARE.
#' @param deg a numeric scalar of polynomial degrees for the kernel smoother
#' @param cv If TRUE, K-fold cross-validation is performed.
#' @param ncores a number of cores to implement \code{foreach()} in \code{doParallel} for K-fold cross-validation.
#' @param cv_bw_loss a loss function for K-fold cross-validation for bandwidth
#' @param cv_K a number of folds for K-fold cross-validation
#' @param ... additional options
#'
#' @return a \code{meanfunc.rob} object contatining follows:
#' \item{t}{a vector containing unlist(Lt)}
#' \item{y}{a vector containing unlist(Ly)}
#' \item{n}{a number of functional trajectories}
#' \item{method}{a method used for obtaining a mean function}
#' \item{kernel}{a kernel used for local linear smoothing}
#' \item{bw}{a bandwidth}
#' \item{delta}{a cuf-off value in M-type loss function. ("HUBER")}
#' \item{deg}{a degree of polnomials}
#' \item{domain}{a range of timepoints}
#' \item{yend}{yend}
#' \item{cv_obj}{an object of K-fold CV for bandwidth}
#' \item{cv_optns}{a list containing options of K-fold cross-validation. (\code{K} = \code{cv_K}, \code{ncore} = \code{ncores}, \code{delta_loss} = \code{cv_delta_loss}, \code{bw_loss} = \code{cv_bw_loss})}
#'
#' @importFrom dplyr %>% group_by summarise
#'
#' @export
#' @useDynLib robfpca
meanfunc.rob <- function(Lt,
                         Ly,
                         newt = NULL,
                         method = c("L2","huber","bisquare"),
                         kernel = "epanechnikov",
                         bw = NULL,
                         delta = NULL,
                         deg = 1,
                         cv = FALSE,
                         ncores = 1,
                         # cv_delta_loss = "L1",
                         cv_bw_loss = "HUBER",
                         cv_K = 5,
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

    # cv <- FALSE
    # # 5-fold CV for delta in Huber function
    # if ((method %in% c("HUBER","BISQUARE")) && is.null(delta)) {
    #     print(paste0("delta is not specified. ", cv_K, "-fold CV is performed for delta in Huber function."))
    #     delta_cv_obj <- delta.locpolysmooth(Lt = Lt,
    #                                         Ly = Ly,
    #                                         method = method,
    #                                         kernel = kernel,
    #                                         deg = deg,
    #                                         # delta = delta,
    #                                         bw = bw,
    #                                         cv_loss = cv_delta_loss,
    #                                         K = cv_K,
    #                                         ncores = ncores,
    #                                         ...)
    #     delta <- delta_cv_obj$selected_delta
    #     cv <- TRUE
    # }

    # 5-fold CV for bandwidth
    if (isTRUE(cv)) {
        print(paste0(cv_K, "-fold CV is performed for bandwidth."))
        bw_cv_obj <- bw.locpolysmooth(Lt = Lt,
                                      Ly = Ly,
                                      method = method,
                                      kernel = kernel,
                                      deg = deg,
                                      delta = delta,
                                      cv_loss = cv_bw_loss,
                                      K = cv_K,
                                      ncores = ncores,
                                      ...)
        bw <- bw_cv_obj$selected_bw
    }

    n <- length(Lt)
    t <- unlist(Lt)
    y <- unlist(Ly)

    ord <- sort(t, index.return = T)$ix
    t <- t[ord]
    y <- y[ord]

    R <- list(t = t,
              y = y,
              n = n,
              method = method,
              kernel = kernel,
              bw = bw,
              delta = delta,
              deg = deg,
              domain = domain,
              yend = c(NULL,NULL))
    class(R) <- 'meanfunc.rob'

    L0 <- domain[2]-domain[1]
    yend <- predict(R, c(domain[1]+L0/100, domain[2]-L0/100))
    R$yend <- yend

    if (!is.null(newt)) {
        R$fitted <- predict(R, newt)
    }

    if (cv == TRUE) {
        R$cv_obj <- bw_cv_obj
        R$cv_optns <- list(K = cv_K,
                           ncores = ncores,
                           # delta_loss = cv_delta_loss,
                           bw_loss = cv_bw_loss)
    }

    return(R)
}




### Predict mean at new time points
#' Predict mean at new time points
#'
#' @param object an \code{meanfunc.rob} object from \code{meanfunc.rob()}
#' @param newt a vector containing timepoints to predict
#' @param ... does not needed
#'
#' @return a vector containing a mean function corresponds to \code{newt}
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

        tmp[ord] <- locpolysmooth(Lt = meanfunc.obj$t,
                                  Ly = meanfunc.obj$y,
                                  newt = newt0[ord],
                                  method = meanfunc.obj$method,
                                  bw = meanfunc.obj$bw,
                                  kernel = meanfunc.obj$kernel,
                                  delta = meanfunc.obj$delta)

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
            dplyr::left_join(df_newt_unique, by = "t") %>%
            dplyr::select(pred) %>%
            unlist()

        return(as.numeric(pred_newt))
        # return(pred(newt))
    }
    else stop('newt must be a vector or a list of vectors of real numbers')
}
