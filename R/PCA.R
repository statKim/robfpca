
#' PCA for functional snippets via conditional expectation
#'
#' @param Lt  a list of vectors or a vector containing time points for all curves
#' @param Ly  a list of vectors or a vector containing observations for all curves
#' @param mu a mean function estimated at work.grid
#' @param cov a covariance function estimated at work.grid
#' @param sig2 a noise variance
#' @param work.grid a work grid
#' @param K a number of PCs. If K is NULL, K is selected by PVE.
#' @param PVE a proportion of variance explained
#'
#' @export
#' @useDynLib robfpca
funPCA <- function(Lt, Ly, mu, cov, sig2 = NULL, work.grid, K = NULL, PVE = 0.99) {
    n <- length(Lt)   # number of observations

    # eigen analysis
    eig.obj <- get_eigen(cov, work.grid)

    # noise variance
    if (is.null(sig2)) {
        sig2 <- 0
    }

    # number of PCs
    if (is.null(K)){
        K <- which(eig.obj$PVE > PVE)[1]
        PVE <- eig.obj$PVE[K]   # PVE
    }

    ## estimate PC scores via conditional expectation
    # - for loop is as fast as sapply!
    PC_score <- matrix(NA, n, K)

    # complete curves - calculate matrix multiplication
    ind_complete <- sapply(Lt, function(t) { identical(work.grid, t) })
    complete_curves <- list2rbind(Ly[ind_complete])   # combine complete curves
    PC_score[ind_complete, ] <- get_CE_score(work.grid, complete_curves,
                                             mu, cov,
                                             sig2, eig.obj, K, work.grid)

    # snippets or partially observed curves - calculate individually
    ind_snippet <- (1:n)[!ind_complete]
    for (i in ind_snippet) {
        PC_score[i, ] <- get_CE_score(Lt[[i]], Ly[[i]],
                                      mu, cov,
                                      sig2, eig.obj, K, work.grid)
    }
    # obs.grid <- sort(unique(unlist(Lt)))
    # PC_score <- sapply(1:n, function(i) {
    #   # xi <- get_CE_score(Lt[[i]], Ly[[i]], mu, cov, sig2, eig.obj, K, work.grid, obs.grid)
    #   xi <- get_CE_score(Lt[[i]], Ly[[i]], mu, cov, sig2, eig.obj, K, work.grid)
    #   return(xi)
    # })

    res <- list(
        lambda = eig.obj$lambda[1:K],
        eig.fun = eig.obj$phi[, 1:K],
        pc.score = PC_score,
        # pc.score = t(PC_score),
        K = K,
        PVE = PVE,
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
#' @param newdata a new data
#' @param K a number of PCs for reconstruction
#' @param ... Not used
#'
#' @importFrom stats predict
#'
#' @export
predict.funPCA <- function(object, newdata = NULL, K = NULL, ...) {
    funPCA.obj <- object
    if (is.null(K)) {
        K <- funPCA.obj$K
    }
    if (K > funPCA.obj$K) {
        stop(paste0("Selected number of PCs from funPCA object is less than K."))
    }

    if (is.null(newdata)) {
        pc.score <- funPCA.obj$pc.score[, 1:K]
        n <- nrow(pc.score)
    } else {
        Lt <- newdata$Lt
        Ly <- newdata$Ly
        n <- length(Lt)

        pc.score <- matrix(NA, n, K)
        for (i in 1:n) {
            pc.score[i, ] <- get_CE_score(Lt[[i]], Ly[[i]],
                                          funPCA.obj$mu, funPCA.obj$cov, funPCA.obj$sig2,
                                          funPCA.obj$eig.obj, K, funPCA.obj$work.grid)
        }
    }

    mu <- matrix(rep(funPCA.obj$mu, n),
                 nrow = n, byrow = TRUE)
    eig.fun <- funPCA.obj$eig.fun[, 1:K]
    pred <- mu + pc.score %*% t(eig.fun)   # reconstructed curves

    return(pred)
}


#' Obtain reconstructed curve for missing parts and simple curve registration
#'
#' @param x a vector containing a partially observed curve (NA for not observed grids)
#' @param pred a vector containing a reconstructed curve for whole grid points
#' @param grid a vector containing grid points (if NULL, it uses equally spaced grids between (0, 1))
#' @param align If TRUE, a simple curve registration are performed for missing parts.
#' @param conti include last observed value, thus looks continuous prediction.
#'
#' @export
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

            # check convexity and difference of slope
            # if (+ => -) or (- => +), other transform is performed.
            # at least 7 obs are required to compare 3rd and 5th
            if (A_i-3 > 0) {
                # missing start part
                first_deriv <- get_deriv(c(x[(A_i-3):A_i],
                                           pred_missing[(A_i+1):(A_i+3)]),
                                         grid[(A_i-3):(A_i+3)])
                is_convex <- is.convex(c(x[(A_i-3):A_i],
                                         pred_missing[(A_i+1):(A_i+3)]),
                                       grid[(A_i-3):(A_i+3)])
            } else {
                # missing end part
                first_deriv <- get_deriv(c(pred_missing[(B_i-3):(B_i-1)],
                                           x[B_i:(B_i+3)]),
                                         grid[(B_i-3):(B_i+3)])
                is_convex <- is.convex(c(pred_missing[(B_i-3):(B_i-1)],
                                         x[B_i:(B_i+3)]),
                                       grid[(B_i-3):(B_i+3)])
            }

            # proportional transform having linear gradient of slope
            # || prod(first_deriv[c(3, 5)]) < 0
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

            # # check first derivative
            # # if (+ => -) or (- => +), other transform is performed.
            # first_deriv <- get_deriv(c(x[(A_i-1):A_i], pred_missing[A_i+1]),
            #                          grid[(A_i-1):(A_i+1)])
            # first_deriv <- sign(first_deriv)   # obtain the sign
            #
            # # proportional transform having linear gradient of slope
            # if (prod(first_deriv[1], first_deriv[3]) < 0) {
            #   # print(paste(ind, ":", prod(first_deriv[1], first_deriv[3])))
            #   g <- function(t, y) {
            #     numerator <- x[B_i] / pred[B_i] - x[A_i] / pred[A_i]
            #     gamma_prime <- numerator / (grid[B_i] - grid[A_i])
            #     const <- x[A_i] / pred[A_i] - gamma_prime*grid[A_i]
            #     gamma <- gamma_prime * t + const
            #
            #     return(gamma * y)
            #   }
            #   pred_missing[A_i:B_i] <- g(grid[(A_i):(B_i)],
            #                              pred[(A_i):(B_i)])
            # }
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
    phi <- eig.obj$phi[, 1:K]
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
        lamda_phi <- diag(lambda) %*% t(phi)
        Sigma_y_mu <- solve(Sigma_y, t(y_mu))
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
    lamda_phi <- diag(lambda) %*% t(phi_y_i)
    Sigma_y_i_mu <- solve(Sigma_y_i, y - mu_y_i)
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

