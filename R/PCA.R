#' Functional Principal Component Analysis (FPCA) via Conditional Expectation
#'
#' FPCA is performed via PACE (Principal Analysis via Conditional Expectation) proposed by Yao et al. (2005).
#'
#' @param Lt  a list of vectors or a vector containing time points for all curves
#' @param Ly  a list of vectors or a vector containing observations for all curves
#' @param mu a mean function estimated at work.grid
#' @param cov a covariance function estimated at work.grid
#' @param sig2 a noise variance
#' @param work.grid a work grid
#' @param K a number of FPCs. If K is NULL, K is selected by PVE.
#' @param PVE a proportion of variance explained
#'
#' @return a list contatining as follows:
#' \item{data}{a list containing Lt and Ly}
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
#'
#' @examples
#' ### Generate example data
#' set.seed(100)
#' x.list <- sim_delaigle(n = 100,
#'                        type = "partial",
#'                        out.prop = 0.2,
#'                        dist = "normal")
#' x <- list2matrix(x.list)
#'
#' ### Estimate robust covariance function
#' work.grid <- seq(0, 1, length.out = 51)
#' cov.obj <- cov_ogk(x,
#'                    type = "huber",
#'                    bw = 0.1)
#'
#' ### Functional principal component analysis
#' fpca.obj <- funPCA(Lt = x.list$Lt,
#'                    Ly = x.list$Ly,
#'                    mu = cov.obj$mean,
#'                    cov = cov.obj$cov,
#'                    work.grid = work.grid,
#'                    PVE = 0.95)
#' fpc.score <- fpca.obj$pc.score
#'
#' # ### Give same result in the above
#' # fpca.obj <- robfpca.partial(x,
#' #                             type = "huber",
#' #                             PVE = 0.95,
#' #                             bw = 0.1)
#' # fpc.score <- fpca.obj$pc.score
#'
#' @export
funPCA <- function(Lt,
                   Ly,
                   mu,
                   cov,
                   sig2 = 0,
                   work.grid = NULL,
                   K = NULL,
                   PVE = 0.99) {
    n <- length(Lt)   # number of observations

    # grid for numerical integration
    if (is.null(work.grid)) {
        time.range <- range(unlist(Lt))
        work.grid <- seq(time.range[1], time.range[2],
                         length.out = nrow(cov))
    }

    # eigen analysis
    eig.obj <- get_eigen(cov, work.grid)

    # fitted covariance which is transformed to positive semi-definite
    cov <- eig.obj$cov_psd

    # # noise variance
    # if (is.null(sig2)) {
    #     sig2 <- 0
    # }

    # number of PCs
    if (is.null(K)){
        K <- which(eig.obj$PVE >= PVE)[1]
        PVE <- eig.obj$PVE[K]   # PVE
    }

    ## estimate PC scores via conditional expectation
    # - for loop is as fast as sapply!
    pc_score <- matrix(NA, n, K)

    ## If there exist complete curves, compute by matrix mutliplication.
    # complete curves - calculate matrix multiplication
    ind_complete <- sapply(Lt, function(t) { identical(work.grid, t) })
    if (sum(ind_complete) > 0) {
        complete_curves <- list2rbind(Ly[ind_complete])   # combine complete curves
        pc_score[ind_complete, ] <- get_CE_score(work.grid,
                                                 complete_curves,
                                                 mu,
                                                 cov,
                                                 sig2,
                                                 eig.obj,
                                                 K,
                                                 work.grid)
        # index of incomplete curves
        ind_snippet <- (1:n)[!ind_complete]
    } else {
        # does not exist complete curves
        ind_snippet <- 1:n
    }

    # snippets or partially observed curves - calculate individually
    for (i in ind_snippet) {
        pc_score[i, ] <- get_CE_score(Lt[[i]],
                                      Ly[[i]],
                                      mu,
                                      cov,
                                      sig2,
                                      eig.obj,
                                      K,
                                      work.grid)
    }


    res <- list(
        data = list(Lt = Lt,
                    Ly = Ly),
        lambda = eig.obj$lambda[1:K],
        eig.fun = matrix(eig.obj$phi[, 1:K],
                         ncol = K),
        pc.score = pc_score,
        K = K,
        PVE = eig.obj$PVE[K],
        work.grid = work.grid,
        eig.obj = eig.obj,
        mu = mu,
        cov = cov,
        sig2 = sig2
    )

    class(res) <- "funPCA"

    return(res)
}


#' Reconstruction via functional PCA
#' newdata should be a list containing Lt and Ly.
#'
#' @param object a \code{funPCA} object from \code{funPCA()}
#' @param newdata a list containing \code{Lt} and \code{Ly}
#' @param K a number of PCs for reconstruction
#' @param ... Not used
#'
#' @method predict funPCA
#'
#' @export
predict.funPCA <- function(object, newdata = NULL, K = NULL, ...) {
    if (class(object) != "funPCA") {
        stop("Check the class of input object! class name 'funPCA' is only supported!")
    }

    if (is.null(K)) {
        K <- object$K
    }
    if (K > object$K) {
        stop(paste0("Selected number of PCs from funPCA object is less than K."))
    }

    if (is.null(newdata)) {
        pc.score <- matrix(object$pc.score[, 1:K],
                           ncol = K)
        n <- nrow(pc.score)
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

    mu <- matrix(rep(object$mu, n),
                 nrow = n, byrow = TRUE)
    eig.fun <- matrix(object$eig.fun[, 1:K],
                      ncol = K)
    pred <- mu + pc.score %*% t(eig.fun)   # reconstructed curves

    return(pred)
}


### Obtain reconstructed curve for missing parts and simple curve registration
pred_missing_curve <- function(x, pred, grid = NULL, align = FALSE, conti = TRUE) {
    num_grid <- length(x)
    if (is.null(grid)) {
        grid <- seq(0, 1, length.out = num_grid)
    }

    pred_missing <- rep(NA, num_grid)   # leave only missing parts

    # Wheter missing is outer only or is in the middle
    # is_snippets == TRUE, missing is outer side.
    # is_snippets == FALSE, missing is in middle.
    is_snippets <- (max( diff( which(!is.na(x)) ) ) == 1)
    if (is_snippets) {
        obs_range <- range(which(!is.na(x)))   # index range of observed periods

        if ((obs_range[1] > 1) & (obs_range[2] < num_grid)) {
            # start and end
            A_i <- obs_range[1] - 1
            B_i <- obs_range[2] + 1

            if (align == TRUE) {
                pred_missing[1:A_i] <- pred[1:A_i] - pred[A_i] + x[A_i + 1]
                pred_missing[B_i:num_grid] <- pred[B_i:num_grid] - pred[B_i] + x[B_i - 1]
            } else {
                pred_missing[1:A_i] <- pred[1:A_i]
                pred_missing[B_i:num_grid] <- pred[B_i:num_grid]
            }

            # add true endpoints
            if (conti == TRUE) {
                pred_missing[c(A_i+1, B_i-1)] <- x[c(A_i+1, B_i-1)]
            }
        } else if ((obs_range[1] > 1) | (obs_range[2] < num_grid)) {
            if (obs_range[1] > 1) {
                # start periods
                A_i <- obs_range[1] - 1

                if (align == TRUE) {
                    pred_missing[1:A_i] <- pred[1:A_i] - pred[A_i] + x[A_i + 1]
                } else {
                    pred_missing[1:A_i] <- pred[1:A_i]
                }

                # add true endpoints
                if (conti == TRUE) {
                    pred_missing[A_i+1] <- x[A_i+1]
                }
            } else if (obs_range[2] < num_grid) {
                # end periods
                B_i <- obs_range[2] + 1

                if (align == TRUE) {
                    pred_missing[B_i:num_grid] <- pred[B_i:num_grid] - pred[B_i] + x[B_i - 1]
                } else {
                    pred_missing[B_i:num_grid] <- pred[B_i:num_grid]
                }

                # add true endpoints
                if (conti == TRUE) {
                    pred_missing[B_i-1] <- x[B_i-1]
                }
            }
        }
    } else {
        # missing is in the middle.
        na_range <- range(which(is.na(x)))   # index range of missing parts

        # last observed point (different with above case!!)
        A_i <- na_range[1] - 1
        B_i <- na_range[2] + 1

        if (align == TRUE) {
            # align gradient(slope) and shift to the true scale
            slope <- (x[B_i] - x[A_i]) / (pred[B_i] - pred[A_i])
            intercept <- x[A_i] - slope * pred[A_i]
            pred_missing[A_i:B_i] <- slope*pred[A_i:B_i] + intercept

            # small missing values => Not align
            if (diff(na_range) >= 3) {
                # check convexity and difference of slope
                # if (+ => -) or (- => +), other transform is performed.
                # at least 7 obs are required to compare 3rd and 5th
                if (A_i - 3 > 0) {
                    # missing start part
                    first_deriv <- get_deriv(c(x[(A_i-3):A_i],
                                               pred_missing[(A_i+1):(A_i+3)]),
                                             grid[(A_i-3):(A_i+3)])
                    is_convex <- is.convex(c(x[(A_i-3):A_i],
                                             pred_missing[(A_i+1):(A_i+3)]),
                                           grid[(A_i-3):(A_i+3)])
                } else if (B_i + 3 <= num_grid) {
                    # missing end part
                    first_deriv <- get_deriv(c(pred_missing[(B_i-3):(B_i-1)],
                                               x[B_i:(B_i+3)]),
                                             grid[(B_i-3):(B_i+3)])
                    is_convex <- is.convex(c(pred_missing[(B_i-3):(B_i-1)],
                                             x[B_i:(B_i+3)]),
                                           grid[(B_i-3):(B_i+3)])
                }

                # proportional transform having linear gradient of slope
                # second derivative is changed or differnece of first derivative > 10
                if ( xor(is_convex[3], is_convex[5]) || abs(first_deriv[3] - first_deriv[5]) > 10 ) {
                    g <- function(t, y) {
                        numerator <- x[B_i] / pred[B_i] - x[A_i] / pred[A_i]
                        gamma_prime <- numerator / (grid[B_i] - grid[A_i])
                        const <- x[A_i] / pred[A_i] - gamma_prime*grid[A_i]
                        gamma <- gamma_prime * t + const

                        return(gamma * y)
                    }
                    pred_missing[A_i:B_i] <- g(grid[(A_i):(B_i)],
                                               pred[(A_i):(B_i)])
                }
            }
        } else {
            pred_missing[A_i:B_i] <- pred[A_i:B_i]
        }

        # remove true endpoints
        if (conti == FALSE) {
            pred_missing[c(A_i, B_i)] <- NA
        }
    }

    return(pred_missing)
}


### Get PC scores via conditional expectation
# - t: observed time points
# - y: matrix => fully observed curves
# - y: vector => snippets or partially observed curve
get_CE_score <- function(t, y, mu, cov, sig2, eig.obj, K, work.grid) {
    phi <- matrix(eig.obj$phi[, 1:K],
                  ncol = K)
    lambda <- eig.obj$lambda[1:K]
    Sigma_y <- cov + diag(sig2, nrow = nrow(cov))
    # Sigma_Y <- fpca.yao$smoothedCov
    # Sigma_Y <- eig.yao$phi %*% diag(eig.yao$lambda) %*% t(eig.yao$phi)


    # # convert phi and fittedCov to obsGrid.
    # mu <- ConvertSupport(work.grid, obs.grid, mu = mu)
    # phi <- ConvertSupport(work.grid, obs.grid, phi = phi)
    # Sigma_Y <- ConvertSupport(work.grid, obs.grid,
    #                           Cov = Sigma_y)
    # # get subsets at each observed grid for conditional expectation
    # mu_y_i <- approx(obs.grid, mu, t)$y
    # phi_y_i <- apply(phi, 2, function(eigvec){
    #   return(approx(obs.grid, eigvec, t)$y)
    # })
    # Sigma_y_i <- matrix(pracma::interp2(obs.grid,
    #                                     obs.grid,
    #                                     Sigma_y,
    #                                     rep(t, each = length(t)),
    #                                     rep(t, length(t))),
    #                     length(t),
    #                     length(t))

    # get CE scores via matrix multiplication for fully observed curves
    if (is.matrix(y)) {
        n <- nrow(y)   # number of complete curves

        # obtain PC score via conditional expectation
        y_mu <- y - matrix(rep(mu, n),
                           nrow = n,
                           byrow = TRUE)
        lamda_phi <- diag(lambda, ncol = K) %*% t(phi)
        # Sigma_y_mu <- solve(Sigma_y, t(y_mu))
        Sigma_y_mu <- MASS::ginv(Sigma_y) %*% t(y_mu)
        xi <- lamda_phi %*% Sigma_y_mu

        return( t(xi) )
    }

    # get subsets at each observed grid for conditional expectation
    if (!identical(work.grid, t)) {
        mu_y_i <- stats::approx(work.grid, mu, t)$y
        phi_y_i <- apply(phi, 2, function(eigvec){
            return(stats::approx(work.grid, eigvec, t)$y)
        })
        Sigma_y_i <- matrix(pracma::interp2(work.grid,
                                            work.grid,
                                            Sigma_y,
                                            rep(t, each = length(t)),
                                            rep(t, length(t))),
                            length(t),
                            length(t))
    } else {
        mu_y_i <- mu
        phi_y_i <- phi
        Sigma_y_i <- Sigma_y
    }

    # obtain PC score via conditional expectation
    lamda_phi <- diag(lambda, ncol = K) %*% t(phi_y_i)
    # Sigma_y_i_mu <- solve(Sigma_y_i, y - mu_y_i)
    Sigma_y_i_mu <- MASS::ginv(Sigma_y_i) %*% matrix(y - mu_y_i, ncol = 1)
    xi <- lamda_phi %*% Sigma_y_i_mu

    return(as.numeric(xi))
}


### Get PC scores using numerical integration
get_IN_score <- function(t, y, mu, cov, sig2, eig.obj, K, work.grid) {
    phi <- eig.obj$phi[, 1:K]
    lambda <- eig.obj$lambda[1:K]
    Sigma_y <- cov + diag(sig2, nrow = nrow(cov))

    # get subsets at each observed grid for numerical integration
    mu_y_i <- stats::approx(work.grid, mu, t)$y
    phi_y_i <- apply(phi, 2, function(eigvec){
        return(stats::approx(work.grid, eigvec, t)$y)
    })

    # obtain PC score using numerical integration
    h <- (y - mu_y_i) %*% phi_y_i
    xi <- fdapace::trapzRcpp(t, h)

    return(as.numeric(xi))
}

