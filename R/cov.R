#########################################################
### Robust covariance estimation for functional snippets
#########################################################

#' Covariance estimation via semiparametric approach
#'
#' @param Lt a list of vectors containing time points for each curve
#' @param Ly a list of vectors containing observations for each curve
#' @param newt a vector containing time points to estimate
#' @param method "huber", "WRM" are supported
#' @param mu \code{meanfunc.rob} object from \code{meanfunc.rob()}
#' @param weig weight
#' @param ... additional options of \code{varfunc.rob()}
#'
#' @export
#' @useDynLib robfpca
covfunc.rob <- function(Lt,
                        Ly,
                        newt = NULL,
                        method = c("Huber","WRM","Bisquare"),
                        mu = NULL,
                        weig = NULL, ...) {

    method <- toupper(method)
    if (!(method %in% c("HUBER","WRM","BISQUARE"))) {
        stop(paste0(method, " is not provided. Check method parameter."))
    }

    return(covfunc.rob.huber(Lt, Ly, mu = mu, newt = newt, method = method, ...))
    # if (method == 'Huber') {
    #   return(cov.huber(Lt, Ly, mu = mu, newt = newt, ...))
    # } else {
    #   stop("method is not supported.")
    # }
}

covfunc.rob.huber <- function(Lt,
                              Ly,
                              newt = NULL,
                              method = "HUBER",
                              domain = NULL,
                              weig = NULL,
                              corf = NULL, # correlation function(theta,x,y)
                              mu = NULL,
                              sig2e = NULL,
                              sig2x = NULL,
                              pfunc = NULL,
                              theta0 = NULL,
                              lb = NULL,
                              ub = NULL,
                              D = NULL,
                              kernel = "epanechnikov", ...) {

    if (is.null(corf)) {
        corf <- function(x, y, theta) {
            matern(x, y, nu = theta)
        }
        D <- 1
    } else {
        if (is.null(theta0) && is.null(D)) {
            stop('The dimension D must be specified')
        }
    }

    if (is.null(mu)) {
        mu <- meanfunc.rob(Lt, Ly, kernel = kernel, method = method)   # Huber option
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
    mu.hat <- unlist(mu.hat)

    cat("Finish mean estimation! \n")

    # if (is.null(sig2e)) {
    #   sig2e <- sigma2.rob(Lt, Ly)
    # }

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
    var.hat <- unlist(var.hat)

    cat("Finish variance estimation! \n")

    if (is.null(domain)) {
        t.vec <- unlist(Lt)
        domain <- c(min(t.vec), max(t.vec))
    }

    th.est <- estimate.theta(Lt,
                             Ly,
                             D = D,
                             var.hat = var.hat,
                             mu.hat = mu.hat,
                             method = 'LS',
                             rho = corf,
                             weig = weig,
                             pfunc = pfunc,
                             theta.lb = lb,
                             theta.ub = ub,
                             theta0 = theta0,
                             domain = domain)$LS

    rslt <- list(sig2e = sig2e,
                 theta = th.est,
                 mu.hat = mu.hat,
                 domain = domain,
                 mu = mu,
                 sig2x = sig2x,
                 rho = function(x, y) { corf(x, y, th.est) },
                 method = method)
    class(rslt) <- 'covfunc.rob'

    if (!is.null(newt)) {
        rslt$fitted <- predict(rslt, newt)
    }
    cat("Finish covariance estimation! \n")

    return(rslt)
}


### Predict covariance at new time points
#' Predict covariance at new time points
#'
#' @param object an object from \code{covfunc.rob()}
#' @param newt new time points to predict
#' @param ... not used
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


### estimate the window width of snippets
estimate.delta <- function(Lt) {
    if (is.list(Lt)) {
        tmp <- lapply(Lt, function(v) max(v)-min(v))
        return(max(unlist(tmp)))
    } else if (is.vector(Lt)) {
        return(max(Lt)-min(Lt))
    } else {
        stop('unsupported type of Lt')
    }
}

