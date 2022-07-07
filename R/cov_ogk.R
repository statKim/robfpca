#' Robust covariance function esimation for partially observed functional data
#'
#' The location and scale functions are computed via pointwise M-estimator, and the covariance function is obtained via robust pairwise computation based on Orthogonalized Gnanadesikan-Kettenring (OGK) estimation.
#' Additionally, bivariate Nadaraya-Watson smoothing is applied for smoothed covariance surfaces.
#'
#' @param X  a n x p matrix. It allows NA.
#' @param type  the option for robust dispersion estimator. "huber", "bisquare", and "tdist" are supported.
#' @param MM the option for M-scale estimator in GK identity. If it is FALSE, the same method using \code{type} is used, that is the iterative algorithm. The closed form solution using method of moments can be used when \code{MM == TRUE}. Defalut is TRUE.
#' @param smooth If it is TRUE, bivariate Nadaraya-Watson smoothing is performed using \code{fields::smooth2d()}. Default is TRUE.
#' @param grid a vector containing the observed timepoints
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
#' \item{mean}{the vector containing the robust mean function.}
#' \item{cov}{a matrix containing robust covariance function.}
#' \item{bw}{a bandwidth of the bivariate smoothing selected from K-fold cross-validation}
#' \item{cv.obj}{cv.obj from bandwidth selection}
#'
#' @examples
#' set.seed(100)
#' x.list <- sim_delaigle(n = 100,
#'                        type = "partial",
#'                        out.prop = 0.2,
#'                        dist = "normal")
#' x <- list2matrix(x.list)
#' cov.obj <- cov_ogk(x,
#'                    type = "huber",
#'                    bw = 0.1)
#' mu.ogk.sm <- cov.obj$mean
#' cov.ogk.sm <- cov.obj$cov
#'
#' @references
#' \cite{Park, Y., Kim, H., & Lim, Y. (2022+). Functional principal component analysis for partially observed elliptical process, Under review.}
#'
#' \cite{Maronna, R. A., & Zamar, R. H. (2002). Robust estimates of location and dispersion for high-dimensional datasets. Technometrics, 44(4), 307-317.}
#'
#' @import stats
#' @importFrom MASS fitdistr
#' @importFrom RobStatTM locScaleM
#'
#' @export
### OGK covariance estimation
# df : degrees of freedom for type = "tdist"
cov_ogk <- function(X,
                    type = c("huber","bisquare","tdist"),
                    MM = TRUE,
                    smooth = TRUE,
                    grid = NULL,
                    bw = NULL,
                    cv = FALSE,
                    df = 3,
                    cv_optns = list(bw_cand = NULL,
                                    K = 5,
                                    ncores = 1)) {
  n <- nrow(X)   # number of curves
  p <- ncol(X)   # number of timepoints

  # Step 2. correlation matrix
  # - Using Rcpp function (It gives same result with the below "cov_gk()".)
  # obj.gk <- cor_gk_cpp(X,
  #                      type = type,
  #                      MM = MM,
  #                      df = df)
  obj.gk <- cov_gk(X,
                   type = type,
                   MM = MM,
                   cor = TRUE,
                   smooth = FALSE,
                   psd = FALSE,
                   df = df)
  U <- obj.gk$cov
  rob.disp <- obj.gk$disp

  # Step 1. scaling
  Y <- sweep(X, 2, rob.disp, "/")
  Y[is.na(Y)] <- 0
  # Y <- sweep(X, 2, rob.disp, "/") %>%
  #   tidyr::replace_na(0)

  # Step 3. spectral decomposition
  eig <- eigen(U)
  E <- eig$vectors

  # Step 4. PC
  Z <- Y %*% E

  # Step 5. compute location and dispersion of Z
  nu <- rep(NA, p)
  Gamma <- matrix(0, p, p)
  for (i in 1:p) {
    if (type %in% c("huber","bisquare")) {
      tmp <- RobStatTM::locScaleM(Z[, i],
                                  psi = type,
                                  eff = 0.95,
                                  maxit = 50,
                                  tol = 1e-04,
                                  na.rm = TRUE)
      nu[i] <- tmp$mu
      Gamma[i, i] <- tmp$disper^2
    } else if (type == "tdist") {
      tmp <- MASS::fitdistr(Z[, i],
                            densfun = "t",
                            df = df)
      nu[i] <- tmp$estimate[1]
      Gamma[i, i] <- tmp$estimate[2]^2
    }
    # Gamma[i, i] <- sd_trim(Z[, i], trim = 0.2)^2
  }

  # Step 6. Transform back to X
  D <- diag(rob.disp)
  A <- D %*% E
  rob.cov <- A %*% Gamma %*% t(A)
  rob.mean <- as.numeric( A %*% matrix(nu, ncol = 1) )

  # # Step 7. Re-weighting
  # # Hard rejection using Beta-quantile chi-squared dist
  # if (reweight == TRUE) {
  #   z.scale <- sqrt(diag(Gamma))
  #   d <- sweep(Z, 2, nu, "-") %>%
  #     sweep(2, z.scale, "/")
  #   d <- rowSums(d^2)   # mahalanobis distance
  #   d0 <- qchisq(beta, p)*median(d) / qchisq(0.5, p)   # cut-off of weight
  #   W <- ifelse(d <= d0, 1, 0)   # weight
  #
  #   # re-weighted mean
  #   X0 <- X %>%
  #     tidyr::replace_na(0)
  #   rob.mean <- as.numeric( matrix(W, nrow = 1) %*% X0 / sum(W) )
  #
  #   # re-weighted covariance
  #   Xmu0 <- sweep(X, 2, rob.mean, "-") %>%
  #     sweep(1, W, "*") %>%
  #     tidyr::replace_na(0)
  #   rob.cov <- t(Xmu0) %*% Xmu0 / sum(W)
  # }


  # # subtract noise variance
  # if (noise.var == TRUE) {
  #   noise.var <- noise_var_gk(X,
  #                             cov = rob.cov)
  # } else {
  #   noise.var <- 0
  # }
  # diag(rob.cov) <- diag(rob.cov) - noise.var


  # 2-dimensional smoothing - does not need to adjust noise variance
  if (smooth == TRUE) {
    # Obtain work.grid
    if (is.null(grid)) {
      work.grid <- seq(0, 1, length.out = p)
    } else {
      if (length(grid) != p) {
        stop("length(grid) does not equal to ncol(X).")
      }
      work.grid <- grid
    }

    # K-fold cross-validation for bandwidth selection
    if (cv == TRUE) {
      cv.obj <- cv.cov_ogk(X,
                           type = type,
                           MM = MM,
                           bw_cand = cv_optns$bw_cand,
                           K = cv_optns$K,
                           ncores = cv_optns$ncores)
      bw <- cv.obj$selected_bw
    } else {
      if (is.null(bw)) {
        bw <- 5 * max(diff(work.grid))
      }
      cv.obj <- NULL
    }

    # Bivariate Nadaraya-Watson smoothing using optimal or given bandwidth
    rob.cov <- fields::smooth.2d(as.numeric(rob.cov),
                                 x = expand.grid(work.grid, work.grid),
                                 surface = F,
                                 theta = bw,
                                 nrow = p,
                                 ncol = p)
    # knots <- min(p/2, 35)   # Remark 3 from Xiao(2013)
    # cov.sm.obj <- refund::fbps(rob.cov,  knots = knots,list(x = gr, z = gr))
    # rob.cov <- cov.sm.obj$Yhat
  } else {
    # Does not perform bivariate smoothing
    bw <- NULL
    cv.obj <- NULL
  }

  return(list(mean = rob.mean,
              cov = rob.cov,
              # noise.var = noise.var,
              bw = bw,
              cv.obj = cv.obj))
}



### Function kappa in equation (3.4)
get_kappa <- function(v) {
  # x <- ifelse(abs(v) <= m, v^2, 2*m*(abs(v) - m/2))   # 2*huber loss

  # 2*Hampel loss
  # cut-off values are obtained from Sinova(2018)
  v <- abs(v)
  v2 <- v[!is.na(v)]
  a <- median(v2)
  b <- quantile(v2, 0.75)
  c <- quantile(v2, 0.85)
  x <- ifelse(v < a, v^2,
              ifelse(v < b, 2*a*v - a^2,
                     ifelse(v < c, a*(v-c)^2/(b-c) + a*(b+c-a), a*(b+c-a))))

  return(x)
}

### Robust loss - scale estimator using Method of moments (See eq (3.4))
get_sigma2_rob <- function(v) {
  x <- get_kappa(v)
  return( mean(x, na.rm = T) )
}


### Gnanadesikan-Kettenring (GK) covariance estimation
# df : degrees of freedom for type = "tdist"
cov_gk <- function(X,
                   type = "huber",
                   MM = FALSE,
                   cor = FALSE,
                   smooth = FALSE,
                   psd = FALSE,
                   # noise.var = FALSE,
                   df = 3) {
  p <- ncol(X)

  # Corss-sectional robust location estimator
  rob.mean <- rep(0, p)
  rob.disp <- rep(0, p)
  for (i in 1:p) {
    if (type %in% c("huber","bisquare")) {
      obj <- RobStatTM::locScaleM(X[, i], psi = type, na.rm = TRUE)
      rob.mean[i] <- obj$mu
      rob.disp[i] <- obj$disper
    } else if (type == "tdist") {
      obj <- MASS::fitdistr(X[which(!is.na(X[, i])), i],
                            densfun = "t",
                            df = df)
      rob.mean[i] <- obj$estimate[1]
      rob.disp[i] <- obj$estimate[2]
    }
    # }  else if (type == "trim") {
    #   rob.mean[i] <- mean(X[which(!is.na(X[, i])), i],
    #                       trim = 0.1)
    #   rob.disp[i] <- chemometrics::sd_trim(X[which(!is.na(X[, i])), i],
    #                                        trim = 0.1)
    # }
  }

  # Compute GK correlation
  cov.gk <- matrix(NA, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i < j) {
        # index of not missing
        ind_not_NA <- which(!is.na(X[, i] + X[, j]))

        if (type %in% c("huber","bisquare")) {   # M-estimator of dispersion
          if (cor == TRUE) {
            # # Scaling to obtain correlation matrix
            # disp <- apply(X[, c(i,j)], 2, function(col){
            #   locScaleM(col[ind_not_NA], psi = type)$disper
            # })
            # z1 <- X[, i]/disp[1] + X[, j]/disp[2]
            # z2 <- X[, i]/disp[1] - X[, j]/disp[2]

            # Scaling to obtain correlation matrix
            obj1 <- RobStatTM::locScaleM(X[ind_not_NA, i], psi = type)
            obj2 <- RobStatTM::locScaleM(X[ind_not_NA, j], psi = type)
            z1 <- (X[, i] - obj1$mu)/obj1$disper + (X[, j] - obj2$mu)/obj2$disper
            z2 <- (X[, i] - obj1$mu)/obj1$disper - (X[, j] - obj2$mu)/obj2$disper

            # # Scaling to obtain correlation matrix
            # z1 <- (X[, i] - rob.mean[i])/rob.disp[i] + (X[, j] - rob.mean[j])/rob.disp[j]
            # z2 <- (X[, i] - rob.mean[i])/rob.disp[i] - (X[, j] - rob.mean[j])/rob.disp[j]
          } else {
            # Not scaling
            z1 <- X[, i] + X[, j]
            z2 <- X[, i] - X[, j]
          }

          if (MM == TRUE) {
            # Closed form of robust loss scale estimator using Method of moments
            # See eq (3.4) in our paper
            z1.disp <- sqrt(get_sigma2_rob(z1))
            if (sd(z2[ind_not_NA]) < 10^(-10)) {
              z2.disp <- 0
            } else {
              z2.disp <- sqrt(get_sigma2_rob(z2))
            }
          } else {
            # Default, using "RobStatTM" package
            z1.disp <- RobStatTM::locScaleM(z1[ind_not_NA],
                                            psi = type)$disper
            if (sd(z2[ind_not_NA]) < 10^(-10)) {
              z2.disp <- 0
            } else {
              z2.disp <- RobStatTM::locScaleM(z2[ind_not_NA],
                                              psi = type)$disper
            }
          }

        } else if (type == "tdist") {   # MLE of t-distribution
          if (cor == TRUE) {
            # Scaling to obtain correlation matrix
            obj1 <- MASS::fitdistr(X[ind_not_NA, i],
                                   densfun = "t",
                                   df = df)
            obj2 <- MASS::fitdistr(X[ind_not_NA, j],
                                   densfun = "t",
                                   df = df)
            z1 <- (X[, i] - obj1$estimate[1])/obj1$estimate[2] + (X[, j] - obj2$estimate[1])/obj2$estimate[2]
            z2 <- (X[, i] - obj1$estimate[1])/obj1$estimate[2] - (X[, j] - obj2$estimate[1])/obj2$estimate[2]
            # disp <- apply(X[, c(i,j)], 2, function(col){
            #   MASS::fitdistr(col[ind_not_NA],
            #                  densfun = "t",
            #                  df = df)$estimate[2]
            # })
            # z1 <- X[, i]/disp[1] + X[, j]/disp[2]
            # z2 <- X[, i]/disp[1] - X[, j]/disp[2]
          } else {
            # Not scaling
            z1 <- X[, i] + X[, j]
            z2 <- X[, i] - X[, j]
          }

          if (MM == TRUE) {
            # Closed form of robust loss scale estimator using Method of moments
            # See eq (3.4) in our paper
            z1.disp <- sqrt(get_sigma2_rob(z1))
            if (sd(z2[ind_not_NA]) < 10^(-10)) {
              z2.disp <- 0
            } else {
              z2.disp <- sqrt(get_sigma2_rob(z2))
            }
          } else {
            # Default, using t-MLE
            z1.disp <- MASS::fitdistr(z1[ind_not_NA],
                                      densfun = "t",
                                      df = df)$estimate[2]
            if (sd(z2[ind_not_NA]) < 10^(-10)) {
              z2.disp <- 0
            } else {
              z2.disp <- MASS::fitdistr(z2[ind_not_NA],
                                        densfun = "t",
                                        df = df)$estimate[2]
            }
          }

        # } else if (type == "trim") {   # Trimmed standard deviation
        #   if (cor == TRUE) {
        #     # Scaling to obtain correlation matrix
        #     disp <- apply(X[, c(i,j)], 2, function(col){
        #       chemometrics::sd_trim(col[ind_not_NA], trim = 0.1)
        #     })
        #     z1 <- X[, i]/disp[1] + X[, j]/disp[2]
        #     z2 <- X[, i]/disp[1] - X[, j]/disp[2]
        #   } else {
        #     # Not scaling
        #     z1 <- X[, i] + X[, j]
        #     z2 <- X[, i] - X[, j]
        #   }
        #
        #   z1.disp <- chemometrics::sd_trim(z1[ind_not_NA], trim = 0.1)
        #   if (sd(z2[ind_not_NA]) < 10^(-10)) {
        #     z2.disp <- 0
        #   } else {
        #     z2.disp <- chemometrics::sd_trim(z2[ind_not_NA], trim = 0.1)
        #   }

        }

        # cov.gk[i, j] <- 0.25*(z1.disp^2 - z2.disp^2)
        cov.gk[i, j] <- (z1.disp^2 - z2.disp^2) / (z1.disp^2 + z2.disp^2)


      } else {
        cov.gk[i, j] <- cov.gk[j, i]
      }

    }
  }
  if (cor == TRUE) {
    diag(cov.gk) <- 1
  } else {
    diag(cov.gk) <- rob.disp^2
  }
  # range(cov.gk)


  rob.cov <- cov.gk

  # # subtract noise variance
  # if (noise.var == TRUE) {
  #   noise.var <- noise_var_gk(X,
  #                             cov = rob.cov)
  # } else {
  #   noise.var <- 0
  # }
  # diag(rob.cov) <- diag(rob.cov) - noise.var

  # # 2-dimensional smoothing - does not need to adjust noise variance
  # if (smooth == T) {
  #   p <- nrow(rob.cov)
  #   knots <- min(p/2, 15)   # Remark 3 from Xiao(2013)
  #   gr <- seq(0, 1, length.out = p)
  #   cov.sm.obj <- refund::fbps(rob.cov,
  #                              knots = knots,
  #                              list(x = gr,
  #                                   z = gr))
  #   rob.cov <- cov.sm.obj$Yhat
  # }

  # make positive-semi-definite
  if (isTRUE(psd)) {
    eig <- eigen(rob.cov)

    # if complex eigenvalues exists, get the real parts only.
    if (is.complex(eig$values)) {
      idx <- which(abs(Im(eig$values)) < 1e-6)
      eig$values <- Re(eig$values[idx])
      eig$vectors <- Re(eig$vectors[, idx])
    }

    k <- which(eig$values > 0)
    lambda <- eig$values[k]
    phi <- matrix(eig$vectors[, k],
                  ncol = length(k))

    rob.cov <- phi %*% diag(lambda, ncol = length(k)) %*% t(phi)

    rob.cov <- (rob.cov + t(rob.cov)) / 2
    # if (length(k) > 1) {
    #     rob.cov <- eig$vectors[, k] %*% diag(eig$values[k]) %*% t(eig$vectors[, k])
    # } else {
    #     rob.cov <- eig$values[k] * (eig$vectors[, k] %*% t(eig$vectors[, k]))
    # }
  }

  return(list(mean = rob.mean,
              cov = rob.cov,
              disp = rob.disp))
}





#' K-fold cross-validation for Bivariate Nadaraya-Watson smoothing for Covariance
#'
#' K-fold cross-validation for Bivariate Nadaraya-Watson smoothing for Covariance
#'
#' @param X  a n x p matrix with or without NA.
#' @param type  the option for robust dispersion estimator. "huber", "bisquare", and "tdist" are supported.
#' @param MM the option for M-scale estimator in GK identity. Default is same method using \code{type} which is iterative algorithm. The closed form solution using method of moments can be used when \code{MM == TRUE}.
#' @param bw_cand  a vector contains the candidates of bandwidths for bivariate smoothing.
#' @param K the number of folds for K-fold cross validation.
#' @param ncores the number of cores on \code{foreach} for parallel computing.
#'
#' @return a list contatining as follows:
#' \item{selected_bw}{the optimal bandwidth selected from the robust K-fold cross-validation.}
#' \item{cv.error}{a matrix containing CV error per bandwidth candidates.}
#'
#' @examples
#' set.seed(100)
#' x.list <- sim_delaigle(n = 100,
#'                        type = "partial",
#'                        out.prop = 0.2,
#'                        dist = "normal")
#' x <- list2matrix(x.list)
#'
#' # 5-fold CV for bivariate smoothing
#' # You can improve the computaion speed by setting "ncores" option.
#' bw_cand = seq(0.01, 0.1, length.out = 10)
#' cov.cv.obj <- cv.cov_ogk(x,
#'                          type = 'huber',
#'                          bw_cand = bw_cand,
#'                          K = 5,
#'                          ncores = 1)
#' print(cov.cv.obj$selected_bw)
#'
#' # Robust smoothed covariance using selected bw from cross-validation
#' cov.obj <- cov_ogk(x,
#'                    type = "huber",
#'                    bw = cov.cv.obj$selected_bw)
#' mu.ogk.sm <- cov.obj$mean
#' cov.ogk.sm <- cov.obj$cov
#' noise.ogk.sm <- cov.obj$noise.var
#'
#' @importFrom foreach %dopar% foreach
#' @importFrom dplyr %>% group_by summarise
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#'
#' @export
### K-fold cross-validation for Bivariate Nadaraya-Watson smoothing for Covariance
### - Not exactly observation-wise cross-validation
### - It is conducted for element-wise covariance
cv.cov_ogk <- function(X,
                       type = c("huber","bisquare","tdist"),
                       MM = TRUE,
                       bw_cand = NULL,
                       K = 5,
                       ncores = 1) {

  if (is.list(X)) {
    gr <- sort(unique(unlist(X$Lt)))
    x <- list2matrix(X)
  } else {
    gr <- seq(0, 1, length.out = ncol(X))
  }

  n <- nrow(X)
  p <- ncol(X)

  # bandwidth candidates
  if (is.null(bw_cand)) {
    bw_min <- max(diff(gr)) / 2
    bw_max <- (max(gr) - min(gr)) / 3
    if (bw_min >= bw_max) {
      bw_max <- (max(gr) - min(gr)) / 2
    }
    bw_cand <- seq(bw_min, bw_max, length.out = 10)
  }

  # obtain the raw covariance (Not smoothed)
  cov_hat <- cov_ogk(X,
                     type = type,
                     MM = MM,
                     smooth = FALSE)
  cov_hat <- cov_hat$cov

  # element-wise covariances
  st <- expand.grid(gr, gr)
  cov_st <- as.numeric(cov_hat)

  # # remove diagonal parts from raw covariance (See Yao et al.(2005))
  # ind <- which(st[, 1] == st[, 2], arr.ind = T)
  # st <- st[-ind, ]
  # cov_st <- cov_st[-ind]

  # get index for each folds
  n <- nrow(st)
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




# ### Yao version
# noise_var_gk <- function(x, gr = NULL, cov = NULL) {
#   m <- ncol(x)
#
#   if (is.null(cov)) {
#     cov <- cov_gk(x,
#                   smooth = FALSE,
#                   noise.var = FALSE)$cov
#   }
#
#   if (is.null(gr)) {
#     gr <- seq(0, 1, length.out = m)
#   }
#   h <- max(diff(gr))
#
#   # 1D smoothing
#   var_y <- diag(cov)
#   var_y <- stats::smooth.spline(gr, var_y)$y
#
#   df <- data.frame(v = as.numeric(cov),
#                    s = rep(gr, m),
#                    t = rep(gr, each = m))
#
#   # 2D smoothing
#   var_x <- rep(NA, m)
#   for (i in 1:m) {
#     idx <- which((abs(df$s - gr[i]) <= h + .Machine$double.eps) &
#                    (abs(df$t - gr[i]) <= h + .Machine$double.eps) &
#                    (df$s != df$t))
#     var_x[i] <- mean(df$v[idx])
#   }
#   diag(cov) <- var_x
#   cov.sm.obj <- refund::fbps(cov,
#                              knots = m/2,   # recommendation of Xiao(2013)
#                              list(x = gr,
#                                   z = gr))
#   # cov.sm.obj <- refund::fbps(cov, list(x = gr,
#   #                                      z = gr))
#   rob.var <- cov.sm.obj$Yhat
#
#   int_inf <- min(gr) + (max(gr) - min(gr)) / 4
#   int_sup <- max(gr) - (max(gr) - min(gr)) / 4
#   idx <- which(gr > int_inf & gr < int_sup)
#   noise_var <- 2 * fdapace::trapzRcpp(gr[idx], var_y[idx] - diag(rob.var)[idx])
#
#   if (noise_var < 0) {
#     noise_var <- 1e-6
#     warning("noise variance is estimated negative")
#   }
#
#   return(noise_var)
# }
