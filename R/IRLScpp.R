#' Local polynomial kernel smoothing with huber loss (mean estimator)
#'
#' Robust local polynomial smoothing via Huber function using iteratively re-weighted least sqruares (IRLS) algorithm
#'
#' @param Lt  a list of vectors or a vector containing time points for all curves
#' @param Ly  a list of vectors or a vector containing observations for all curves
#' @param newt a vector containing time points to estimate
#' @param method "Huber" or "WRM" or "Bisquare"
#' @param bw bandwidth
#' @param kernel a kernel function for kernel smoothing ("epan", "gauss" are supported.)
#' @param deg a degree of polynomial
#' @param ncores number of cores for k-fold cross-validation
#' @param k2 cut-off value for Huber function
#' For loss = "Huber", it uses \code{rlm()} in \code{MASS} and fits the robust regression with Huber loss.
#' So additional parameters of \code{rlm()} can be applied. (k2, maxit, ...)
#' @param ... additional parameters
#'
#' @return a smoothed vector
#' @export
local_kern_smooth <- function(Lt, Ly, newt = NULL, method = c("L2","HUBER","WRM","BISQUARE"),
                              bw = NULL, deg = 1, ncores = 1,
                              kernel = "epanechnikov", k2 = 1.345, ...) {
  method <- toupper(method)
  if (!(method %in% c("L2","HUBER","WRM","BISQUARE"))) {
    stop(paste0(method, " is not provided. Check method parameter."))
  }

  # If `bw` is not defined, 5-fold CV is performed.
  if (is.null(bw)) {
    if (!(is.list(Lt) & is.list(Ly))) {
      stop("Lt or Ly are not list type. If bw is NULL, 5-fold CV are performed but it is needed list type.")
    }
    bw <- bw.local_kern_smooth(Lt = Lt, Ly = Ly, method = method, kernel = kernel,
                               ncores = ncores, k2 = k2)
  }

  if (is.list(Lt) | is.list(Ly)) {
    Lt <- unlist(Lt)
    Ly <- unlist(Ly)
  }

  if (is.null(newt)) {
    newt <- Lt
  }
  if (is.list(newt)) {
    newt <- unlist(newt)
  }

  if (method %in% c("HUBER","BISQUARE")) {   # proposed Huber loss
    mu_hat <- locpolysmooth(Lt = Lt,
                            Ly = Ly,
                            newt = newt,
                            kernel = kernel,
                            bw = bw,
                            k = k2,
                            deg = deg)
  } else if (method == "L2") {   # squared loss
    w <- 1/length(Lt)
    mu_hat <- sapply(newt, function(t) {
      tmp <- (Lt - t) / bw

      if (kernel == "epanechnikov") {
        kern <- (3/4) * (1 - tmp^2)   # Epanechnikov kernel
      } else if (kernel == "gauss") {
        kern <- 1/sqrt(2*pi) * exp(-1/2 * tmp^2)   # gaussian kernel
      }

      idx <- which(kern > 0)   # non-negative values
      W <- w * kern[idx] / bw   # weight vector for w_i * K_h(x)
      # W <- diag(w * kern[idx]) / bw   # weight vector for w_i * K_h(x)
      X <- matrix(1, length(idx), deg+1)
      for (d in 1:deg) {
        X[, d+1] <- (Lt[idx] - t)^d
      }
      Y <- Ly[idx]
      # print(paste(t, ":", length(kern[idx])))   # number of non-negative weights

      # Weighted least squares
      beta <- solve(t(X) %*% diag(W) %*% X) %*% t(X) %*% diag(W) %*% Y

      return(beta[1, ])
    })
  } else if (method == "WRM") {   # robfilter package
    # robfilter
    if (kernel == "epanechnikov") {
      kern <- 2
    } else if (kernel == "gauss") {
      kern <- 3
    }
    wrm.obj <- robfilter::wrm.smooth(Lt,
                                     Ly,
                                     h = bw,
                                     xgrid = newt,
                                     weight = kern)
    mu_hat <- wrm.obj$level

    # # C++
    # wrm.obj <- wrm_smooth(x = Lt,
    #                       y = Ly,
    #                       h = bw,
    #                       xgrid = newt,
    #                       kernel = kernel)
    # mu_hat <- wrm.obj$mu
  }

  return( as.numeric(mu_hat) )
}



#' K-fold cross validation to find optimal bandwidth for local polynomial kernel smoother
#'
#' @param Lt a list of vectors containing time points for each curve
#' @param Ly a list of vectors containing observations for each curve
#' @param method "Huber" or "WRM" or "Bisquare"
#' @param kernel a kernel function for kernel smoothing ("epan", "gauss" are supported.)
#' @param cv_loss "Huber" or "L1" or "L2"
#' @param k2 cut-off value for Huber function
#' @param bw_cand user defined bandwidth candidates for CV
#' @param K the number of folds
#' @param ncores If ncores > 1, it implements \code{foreach()} in \code{doParallel} for CV.
#' @param ... parameters are same with \code{local_kern_smooth()}.
#'
#' @return bandwidth
#' @import foreach
#' @import dplyr
#'
#' @export
bw.local_kern_smooth <- function(Lt, Ly, method = "HUBER", kernel = "epanechnikov",
                                 cv_loss = "HUBER", k2 = 1.345, bw_cand = NULL,
                                 K = 5,  ncores = 1, ...) {
  cv_loss <- toupper(cv_loss)
  if (!(cv_loss %in% c("HUBER","L1","L2"))) {
    stop(paste0(cv_loss, " is not provided. Check cv_loss parameter."))
  }

  if (!(is.list(Lt) & is.list(Ly))) {
    stop("Lt and Ly should be only a list type.")
  }

  if (is.null(bw_cand)) {
    a <- min(unlist(Lt))
    b <- max(unlist(Lt))
    bw_cand <- 10^seq(-2, 0, length.out = 10) * (b - a)/3

    if (kernel == "epanechnikov") {
      bw_cand <- 10^seq(-1, 0, length.out = 10) * (b - a)/3
    }
  }

  # get index for each folds
  folds <- list()
  n <- length(Lt)   # the number of curves
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
      ncores <- parallel::detectCores() - 3
      warning(paste0("ncores is too large. We now use", ncores, " cores."))
    }
    # ncores <- detectCores() - 3
    cl <-parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    # matrix of bw_cand and fold
    bw_fold_mat <- data.frame(bw_cand = rep(bw_cand, each = K),
                              fold = rep(1:K, length(bw_cand)))

    cv_error <- foreach::foreach(i = 1:nrow(bw_fold_mat), .combine = "c",
                                 .export = c("local_kern_smooth"),
                                 # .noexport = c("locpolysmooth","IRLScpp","get_positive_elements"),
                                 .packages = c("robfilter"),
                                 .errorhandling = "pass") %dopar% {

      bw <- bw_fold_mat$bw_cand[i]   # bandwidth candidate
      k <- bw_fold_mat$fold[i]   # fold for K-fold CV

      # data of kth fold
      Lt_train <- Lt[ -folds[[k]] ]
      Ly_train <- Ly[ -folds[[k]] ]
      Lt_test <- Lt[ folds[[k]] ]
      Ly_test <- Ly[ folds[[k]] ]

      y_hat <- local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
                                 bw = bw, kernel = kernel, k2 = k2, ...)
      # y_hat <- tryCatch({
      #   local_kern_smooth_cpp(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
      #                         bw = bw, kernel = kernel, k2 = k2, ...)
      # }, error = function(e) {
      #   return(NA)
      # })
      # # if error occurs in kernel smoothing, return Inf
      # if (is.na(y_hat)) {
      #   return(Inf)
      # }

      y <- unlist(Ly_test)
      if (cv_loss == "L2") {   # squared errors
        err <- sum((y - y_hat)^2)
      } else if (cv_loss == "HUBER") {   # Huber errors
        a <- abs(y - y_hat)
        err_huber <- ifelse(a > k2, k2*(a - k2/2), a^2/2)
        err <- sum(err_huber)
      } else if (cv_loss == "L1") {   # absolute errors
        err <- sum(abs(y - y_hat))
      }

      return(err)
    }
    parallel::stopCluster(cl)

    bw_fold_mat$cv_error <- cv_error
    cv_obj <- bw_fold_mat %>%
      group_by(bw_cand) %>%
      summarise(cv_error = sum(cv_error))

    bw <- list(selected_bw = cv_obj$bw_cand[ which.min(cv_obj$cv_error) ],
               cv.error = as.data.frame(cv_obj))
  } else {
    cv_error <- rep(0, length(bw_cand))
    for (k in 1:K) {
      Lt_train <- Lt[ -folds[[k]] ]
      Ly_train <- Ly[ -folds[[k]] ]
      Lt_test <- Lt[ folds[[k]] ]
      Ly_test <- Ly[ folds[[k]] ]

      for (i in 1:length(bw_cand)) {
        y_hat <- local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
                                   bw = bw_cand[i], kernel = kernel, k2 = k2, ...)
        # y_hat <- tryCatch({
        #   local_kern_smooth_cpp(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
        #                         bw = bw_cand[i], kernel = kernel, k2 = k2, ...)
        # }, error = function(e) {
        #   return(NA)
        # })
        # # if error occurs in kernel smoothing, return Inf
        # if (is.na(y_hat)) {
        #   return(Inf)
        # }
        # if (i == 1 && k == 1) {
        #   print(y_hat)
        # }

        y <- unlist(Ly_test)
        if (cv_loss == "L2") {
          cv_error[i] <- cv_error[i] + sum((y - y_hat)^2)   # squared errors
        } else if (cv_loss == "HUBER") {
          a <- abs(y - y_hat)
          err_huber <- ifelse(a > k2, k2*(a - k2/2), a^2/2)
          cv_error[i] <- cv_error[i] + sum(err_huber)
        } else if (cv_loss == "L1") {   # absolute errors
          cv_error[i] <- cv_error[i] + sum(abs(y - y_hat))
        }
      }
    }

    bw <- list(selected_bw = bw_cand[ which.min(cv_error) ],
               cv.error = data.frame(bw = bw_cand,
                                     error = cv_error))
  }

  return(bw)
}



#' K-fold cross validation to find optimal delta for Huber loss function
#'
#' @param Lt a list of vectors containing time points for each curve
#' @param Ly a list of vectors containing observations for each curve
#' @param method "Huber" or "WRM" or "Bisquare"
#' @param kernel a kernel function for kernel smoothing ("epan", "gauss" are supported.)
#' @param cv_loss "L1" or "L2"
#' @param bw bandwidth
#' @param delta_cand user defined delta candidates for CV
#' @param K the number of folds
#' @param ncores If ncores > 1, it implements \code{foreach()} in \code{doParallel} for CV.
#' @param ... parameters are same with \code{local_kern_smooth()}.
#'
#' @return delta
#' @import foreach
#' @import dplyr
#' @import doParallel
#' @import parallel
#'
#' @export
delta.local_kern_smooth <- function(Lt, Ly, method = "HUBER", kernel = "epanechnikov",
                                    cv_loss = "L1", bw = NULL, delta_cand = NULL,
                                    K = 5, ncores = 1, ...) {
  cv_loss <- toupper(cv_loss)
  if (!(cv_loss %in% c("L1","L2"))) {
    stop(paste0(cv_loss, " is not provided. Check cv_loss parameter."))
  }

  if (!(is.list(Lt) & is.list(Ly))) {
    stop("Lt and Ly should be only a list type.")
  }

  if (is.null(delta_cand)) {
    a <- min(unlist(Lt))
    b <- max(unlist(Lt))
    delta_cand <- union(1.345,
                        10^seq(-3, 2, length.out = 10) * (b - a)/3)
  }

  # fixed bandwidth
  if (is.null(bw)) {
    a <- min(unlist(Lt))
    b <- max(unlist(Lt))
    bw <- (b - a) / 5
  }

  # get index for each folds
  folds <- list()
  n <- length(Lt)   # the number of curves
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
      ncores <- parallel::detectCores() - 3
      warning(paste0("ncores is too large. We now use", ncores, " cores."))
    }
    # ncores <- detectCores() - 3
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    # matrix of delta_cand and fold
    delta_fold_mat <- data.frame(delta_cand = rep(delta_cand, each = K),
                                 fold = rep(1:K, length(delta_cand)))

    cv_error <- foreach::foreach(i = 1:nrow(delta_fold_mat), .combine = "c",
                                 # .export = c("local_kern_smooth"),
                                 .packages = c("robfilter","robfpca"),
                                 .errorhandling = "pass") %dopar% {
      delta <- delta_fold_mat$delta_cand[i]   # bandwidth candidate
      k <- delta_fold_mat$fold[i]   # fold for K-fold CV

      # data of kth fold
      Lt_train <- Lt[ -folds[[k]] ]
      Ly_train <- Ly[ -folds[[k]] ]
      Lt_test <- Lt[ folds[[k]] ]
      Ly_test <- Ly[ folds[[k]] ]

      y_hat <- robfpca::local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
                                          bw = bw, k2 = delta, kernel = kernel, ...)
      # y_hat <- tryCatch({
      #   local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
      #                     bw = bw, kernel = kernel, k2 = delta, ...)
      # }, error = function(e) {
      #   return(NA)
      # })
      # # if error occurs in kernel smoothing, return Inf
      # if (is.na(y_hat)) {
      #   return(Inf)
      # }

      y <- unlist(Ly_test)
      if (cv_loss == "L2") {   # squared errors
        err <- sum((y - y_hat)^2)
        # } else if (cv_loss == "HUBER") {   # Huber errors
        #   a <- abs(y - y_hat)
        #   err_huber <- ifelse(a > k2, k2*(a - k2/2), a^2/2)
        #   err <- sum(err_huber)
      } else if (cv_loss == "L1") {   # absolute errors
        err <- sum(abs(y - y_hat))
      }

      return(err)
    }
    parallel::stopCluster(cl)

    delta_fold_mat$cv_error <- cv_error
    cv_obj <- delta_fold_mat %>%
      group_by(delta_cand) %>%
      summarise(cv_error = sum(cv_error))

    delta <- list(selected_delta = cv_obj$delta_cand[ which.min(cv_obj$cv_error) ],
                  cv.error = as.data.frame(cv_obj))
  } else {
    cv_error <- rep(0, length(delta_cand))
    for (k in 1:K) {
      Lt_train <- Lt[ -folds[[k]] ]
      Ly_train <- Ly[ -folds[[k]] ]
      Lt_test <- Lt[ folds[[k]] ]
      Ly_test <- Ly[ folds[[k]] ]

      for (i in 1:length(delta_cand)) {
        y_hat <- local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, method = method,
                                   bw = bw, k2 = delta_cand[i], kernel = kernel, ...)
        y <- unlist(Ly_test)
        if (cv_loss == "L2") {
          cv_error[i] <- cv_error[i] + sum((y - y_hat)^2)   # squared errors
        } else if (cv_loss == "L1") {
          cv_error[i] <- cv_error[i] + sum(abs(y - y_hat))   # LAD
          # } else if (cv_loss == "HUBER") {
          #   a <- abs(y - y_hat)
          #   err_huber <- ifelse(a > k2, k2*(a - k2/2), a^2/2)
          #   cv_error[i] <- cv_error[i] + sum(err_huber)
        }
      }
    }

    delta <- list(selected_delta = delta_cand[ which.min(cv_error) ],
                  cv.error = data.frame(delta = delta_cand,
                                        error = cv_error))
  }

  return(delta)
}

