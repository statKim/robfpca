#########################################################
### Robust covariance estimation for functional snippets
#########################################################

#' Covariance estimation via semiparametric approach
#'
#' @param Lt a list of vectors containing time points for each curve
#' @param Ly a list of vectors containing observations for each curve
#' @param newt a vector containing time points to estimate
#' @param method "huber", "WRM" are supported
#' @param mu Optional, \code{meanfunc.rob} object from \code{meanfunc.rob()}
#' @param sig2x Optional, \code{varfunc.rob} object from \code{varfunc.rob()}
#' @param sig2e Optional, a noise variance estimator obtained from \code{sigma2.rob()}
#' @param corf a correlation structure, it should be "function" class(x, y). Default is matern correlation.
#' @param kernel a kernel function for local polynomial smoothing ("epanechnikov", "gauss" are supported.)
#' @param bw a bandwidth for local polynomial smoothing.
#' @param delta cut-off value for "huber"(Huber) or "bisquare"(Tukey's biweight function).
#' Default is 1.345 for "huber" and 4.685 for "bisquare" for 95\% ARE.
#' @param deg a numeric scalar of polynomial degrees for local polynomial smoothing
#' @param cv If TRUE, K-fold cross-validation is performed for bandwidth.
#' @param ncores a number of cores to implement \code{foreach()} in \code{doParallel} for K-fold cross-validation.
#' @param cv_bw_loss a loss function for K-fold cross-validation
#' @param cv_K a number of folds for K-fold cross-validation
#' @param ... additional options
#'
#' @return a \code{covfunc.rob} object contatining as follows:
#' \item{sig2}{a noise variance estimate}
#' \item{theta}{a estiamted parameter in parametric correlation function}
#' \item{mu.hat}{a estimated mean function corresponded to unlist(Lt)}
#' \item{domain}{a range of timepoints}
#' \item{mu}{a \code{meanfunc.rob} object}
#' \item{sig2x}{a \code{varfunc.rob} object}
#' \item{rho}{a parametric correlation function}
#' \item{method}{a method used for obtaining a mean and variance function}
#'
#' @export
#' @useDynLib robfpca
covfunc.rob <- function(Lt,
                        Ly,
                        newt = NULL,
                        method = c("huber","bisquare"),
                        mu = NULL,
                        sig2x = NULL,
                        sig2e = NULL,
                        corf = NULL,
                        kernel = "epanechnikov",
                        bw = NULL,
                        delta = NULL,
                        deg = 1,
                        cv = FALSE,
                        ncores = 1,
                        cv_bw_loss = "HUBER",
                        cv_K = 5,
                        ...) {

    method <- toupper(method)
    if (!(method %in% c("HUBER","WRM","BISQUARE"))) {
        stop(paste0(method, " is not provided. Check method parameter."))
    }

    ### Mean estimation
    if (is.null(mu)) {
        mu <- meanfunc.rob(Lt = Lt,
                           Ly = Ly,
                           newt = NULL,
                           method = method,
                           kernel = kernel,
                           bw = bw,
                           delta = delta,
                           deg = deg,
                           cv = cv,
                           ncores = ncores,
                           cv_bw_loss = cv_bw_loss,
                           cv_K = cv_K)
        # mu <- meanfunc.rob(Lt, Ly, kernel = kernel, method = method)   # Huber option
    }
    # mu.hat <- predict(mu, unlist(Lt))   #  TOO SLOW => Need to improve speed
    n <- length(Lt)
    gr <- sort(unique(unlist(Lt)))
    mu.hat.unique <- predict(mu, gr)
    mu.hat <- lapply(1:n, function(i) {
        ind <- match(Lt[[i]], gr)
        if (length(ind) == length(Lt[[i]])) {
            return(mu.hat.unique[ind])
        } else {
            mui <- predict(mu, Lt[[i]])
            return(mui)
        }
    })
    # mu.hat <- unlist(mu.hat)

    # cat("Finish mean estimation! \n")


    ### Variance estimation
    if (is.null(sig2x)) {
        # sig2x <- varfunc(Lt,Ly,mu=mu,sig2=sig2e)
        sig2x <- varfunc.rob(Lt, Ly, mu = mu, sig2 = sig2e, kernel = kernel,
                             method = method, ...)   # Huber option
    }
    sig2e <- sig2x$sig2   # noise variance

    # var.hat <- predict(sig2x, unlist(Lt))   #  TOO SLOW => Need to improve speed
    var.hat.unique <- predict(sig2x, gr)
    var.hat <- lapply(1:n, function(i) {
        ind <- match(Lt[[i]], gr)
        if (length(ind) == length(Lt[[i]])) {
            return(var.hat.unique[ind])
        } else {
            vari <- predict(sig2x, Lt[[i]])
            return(vari)
        }
    })
    # var.hat <- unlist(var.hat)

    # cat("Finish variance estimation! \n")



    ### Covariance (off-diagonal parts) estimation
    if (is.null(corf)) {
        corf <- function(x, y, theta) {
            matern(x, y, nu = theta)
        }
        D <- 1
    }

    # range of timepoint
    domain <- range(unlist(Lt))

    th.est <- estimate.theta(Lt,
                             Ly,
                             D = D,
                             var.hat = var.hat,
                             mu.hat = mu.hat,
                             method = 'LS',
                             rho = corf,   # correlation function(theta,x,y)
                             weig = NULL,
                             pfunc = NULL,
                             theta.lb = NULL,
                             theta.ub = NULL,
                             theta0 = NULL,
                             domain = domain)$LS

    cov.obj <- list(sig2e = sig2e,
                    theta = th.est,
                    mu.hat = mu.hat,
                    domain = domain,
                    mu = mu,
                    sig2x = sig2x,
                    rho = function(x, y) { corf(x, y, th.est) },
                    method = method)
    class(cov.obj) <- 'covfunc.rob'

    if (!is.null(newt)) {
        cov.obj$fitted <- predict(cov.obj, newt)
    }
    # cat("Finish covariance estimation! \n")

    return(cov.obj)
}



#' Predict covariance at new time points
#'
#' @param object an object from \code{covfunc.rob()}
#' @param newt a vector containing timepoints to predict
#' @param ... does not needed
#'
#' @return a matrix containing a covariance function corresponds to \code{newt}
#'
#' @importFrom stats predict
#'
#' @export
predict.covfunc.rob <- function(object, newt, ...) {
    covobj <- object
    pred <- function(newt) {  # newt is a vector
        stopifnot(is.vector(newt))

        corr.est <- covobj$rho(newt, newt)
        var.hat <- predict(covobj$sig2x, newt)
        sig.hat <- sqrt(var.hat)
        cov.fitted <- corr.est * (sig.hat %*% t(sig.hat))
        return(cov.fitted)
    }

    if (is.list(newt)) {
        return(lapply(newt, pred))
    } else if (is.vector(newt)) {
        return(pred(newt))
    } else {
        stop('newt must be a vector or a list of vectors of real numbers')
    }
}


# ### estimate the window width of snippets
# estimate.delta <- function(Lt) {
#     if (is.list(Lt)) {
#         tmp <- lapply(Lt, function(v) max(v)-min(v))
#         return(max(unlist(tmp)))
#     } else if (is.vector(Lt)) {
#         return(max(Lt)-min(Lt))
#     } else {
#         stop('unsupported type of Lt')
#     }
# }





####################################################
### Functions from mcfda package
####################################################
# matern correlation
matern <- function(x,y=NULL,nu=1,rho=1)
{
    if(is.null(y)) y <- x

    G <- expand.grid(x,x)

    S <- apply(G,1,function(z){
        delta <- abs(z[1]-z[2])/rho
        if(delta == 0)
            1
        else
            (sqrt(2*nu)*delta)^nu * besselK(sqrt(2*nu)*delta,nu=nu) / (2^(nu-1)*gamma(nu))
    })

    C <- matrix(S,nrow=length(x),ncol=length(y),byrow=F)
    return(C)
}


# Estimate the parameters of the correlation structure
# @param weig a vector of length(Lt)
# @param D the number of parameters
# @param sig0.hat the estimate of variance of measure error
# @keywords internal
estimate.theta <- function(Lt,Ly,
                           D,
                           var.hat,
                           mu.hat,
                           method,
                           rho,
                           weig=NULL,
                           theta.lb=NULL,
                           theta.ub=NULL,
                           h=0.1,
                           domain=c(0,1),
                           pfunc=NULL,
                           theta0=NULL) {
    if(is.null(weig))   weig <- rep(1,length(Lt))
    if(is.null(pfunc)) pfunc <- function(th){return(0)}

    sig <- lapply(var.hat,function(x){
        x <- ifelse(x <=0,yes=1e-8,no=x)
        sqrt(x)
    })

    n <- length(Lt)
    result <- list()

    if ('LS' %in% method) {
        R <- sapply(1:n,function(i){
            # tobs <- Lt[[i]]
            # idx <- (tobs >= domain[1]+h) & (tobs <= domain[2]-h)
            # yy <- Ly[[i]]
            # yy <- yy[idx]
            # mu <- mu.hat[[i]]
            # mu <- mu[idx]
            resid <- Ly[[i]] - mu.hat[[i]]
            return(list(cbind(resid) %*% rbind(resid)))
        })

        TT <- Lt

        SS <- sapply(1:n,function(i){
            # tobs <- Lt[[i]]
            # idx <- (tobs >= domain[1]+h) & (tobs <= domain[2]-h)
            sig.t <- sig[[i]]
            # sig.t <- sig.t[idx]
            diag(sig.t)
        })
        # R <- sapply(1:n,function(i){
        #     tobs <- Lt[[i]]
        #     idx <- (tobs >= domain[1]+h) & (tobs <= domain[2]-h)
        #     yy <- Ly[[i]]
        #     yy <- yy[idx]
        #     mu <- mu.hat[[i]]
        #     mu <- mu[idx]
        #     resid <- yy-mu
        #     return(list(cbind(resid) %*% rbind(resid)))
        # })
        #
        # TT <- sapply(1:n,function(i){
        #     tobs <- Lt[[i]]
        #     idx <- (tobs >= domain[1]+h) & (tobs <= domain[2]-h)
        #     tobs[idx]
        # })
        #
        # SS <- sapply(1:n,function(i){
        #     tobs <- Lt[[i]]
        #     idx <- (tobs >= domain[1]+h) & (tobs <= domain[2]-h)
        #     sig.t <- sig[[i]]
        #     sig.t <- sig.t[idx]
        #     diag(sig.t)
        # })

        Q <- function(theta){
            if(any(theta < theta.lb) || any(theta > theta.ub)) return(1e300)
            v <- sapply(1:n,function(i){
                tobs <- TT[[i]]
                if(length(tobs)==1) return(0)

                rho.mat <- rho(tobs,tobs,theta)
                if(any(is.infinite(rho.mat)) || any(is.nan(rho.mat)))
                {
                    theta
                    tobs
                }
                S <- SS[[i]]
                # print(tobs)
                # print(dim(S))
                # print(dim(rho.mat))
                resid <- (S %*% rho.mat %*% S - R[[i]])^2
                diag(resid) <- 0 # remove points on diagonal
                return(sum(resid))
            })

            if(!is.finite(sum(v*weig)))    return(1e300)
            return(sum(v*weig) + pfunc(theta)) # regularization
        }

        # set up initial value for optimization
        if(is.null(theta.lb)) theta.lb <- rep(-Inf,D)
        if(is.null(theta.ub)) theta.ub <- rep(Inf,D)

        if(is.null(theta0))
            v <- sapply(1:D,function(d){
                lb <- theta.lb[d]
                ub <- theta.ub[d]
                if(lb > -Inf && ub < Inf) return(stats::runif(1)*(ub-lb)+lb)
                else if(lb > -Inf) return(stats::runif(1)+lb)
                else if(ub < Inf) return(ub-stats::runif(1))
                else return(stats::runif(1))
            })
        else v <- theta0
        res <- stats::optim(v,Q,lower=theta.lb,upper=theta.ub,method='L-BFGS-B')
        result$LS <- res$par
        result$Q <- Q(res$par)
    }

    # if('QMLE' %in% method)
    # {
    #     R <- sapply(1:n,function(i){
    #         tobs <- Lt[[i]]
    #         idx <- (tobs >= domain[1]+h) & (tobs <= domain[2]-h)
    #         yy <- Ly[[i]]
    #         yy <- yy[idx]
    #         mu <- mu.hat[[i]]
    #         mu <- mu[idx]
    #         resid <- yy-mu
    #         return(list(cbind(resid) %*% rbind(resid)))
    #     })
    #
    #     TT <- sapply(1:n,function(i){
    #         tobs <- Lt[[i]]
    #         idx <- (tobs >= domain[1]+h) & (tobs <= domain[2]-h)
    #         tobs[idx]
    #     })
    #
    #     SS <- sapply(1:n,function(i){
    #         tobs <- Lt[[i]]
    #         idx <- (tobs >= domain[1]+h) & (tobs <= domain[2]-h)
    #         sig.t <- sig[[i]]
    #         sig.t <- sig.t[idx]
    #         diag(1/sig.t)
    #     })
    #
    #
    #     # set up initial value for optimization
    #     if(is.null(theta.lb)) theta.lb <- rep(-Inf,D)
    #     if(is.null(theta.ub)) theta.ub <- rep(Inf,D)
    #
    #     Q <- function(theta){
    #         if(any(theta < theta.lb) || any(theta > theta.ub)) return(1e300)
    #         v <- sapply(1:n,function(i){
    #             tobs <- TT[[i]]
    #             m <- length(tobs)
    #             rho.mat <- rho(theta,tobs,tobs)
    #             if(any(is.infinite(rho.mat)) || any(is.nan(rho.mat)))
    #             {
    #                 theta
    #                 tobs
    #             }
    #             if(pracma::cond(rho.mat) > 1e2)
    #             {
    #                 rho.inv <- pracma::inv(rho.mat+0.01*diag(rep(1,m)))
    #             }
    #             else rho.inv <- pracma::inv(rho.mat)
    #             S <- SS[[i]]
    #             resid <- (t(R[[i]]) %*% S %*% rho.inv %*% S %*% R[[i]])
    #             resid <- resid + log(abs(det(rho.mat)))
    #             return(resid/2)
    #         })
    #         rs <- sum(v*weig)
    #         if(is.nan(rs) || is.infinite(rs))
    #         {
    #             rs <- 1e100
    #         }
    #         return(rs)
    #     }
    #
    #
    #     v <- sapply(1:D,function(d){
    #         lb <- theta.lb[d]
    #         ub <- theta.ub[d]
    #         if(lb > -Inf && ub < Inf) return(stats::runif(1)*(ub-lb)+lb)
    #         else if(lb > -Inf) return(stats::runif(1)+lb)
    #         else if(ub < Inf) return(ub-stats::runif(1))
    #         else return(stats::runif(1))
    #     })
    #     res <- stats::optim(v,Q,lower=theta.lb,upper=theta.ub,method='L-BFGS-B')
    #     result$QMLE <- res$par
    # }

    return(result)
}
