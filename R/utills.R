##########################################
### Utill functions
##########################################

#' Calculate Intergrated Squared Errors (ISE)
#'
#' @param x vectors or matrices to compare (y-axis)
#' @param x_hat vectors or matrices to compare (y-axis)
#' @param grid corresponding observed grid (x-axis)
#'
#' @return a value of ISE
#'
#' @export
get_ise <- function(x, x_hat, grid) {
    z <- (x - x_hat)^2

    # fdapace package
    if (is.matrix(z)) {   # 2-dim matrix (covariance matrix)
        row.ise <- apply(z, 1, function(row){
            fdapace::trapzRcpp(grid, row)
        })
        ise <- fdapace::trapzRcpp(grid, row.ise)
    } else {   # 1-dim vector
        ise <- fdapace::trapzRcpp(grid, z)
    }

    # # pracma package
    # if (is.matrix(z)) {   # 2-dim matrix (covariance matrix)
    #   row.ise <- apply(z, 1, function(row){
    #     trapz(grid, row)
    #   })
    #   ise <- trapz(grid, row.ise)
    # } else {   # 1-dim vector
    #   ise <- trapz(grid, z)
    # }

    return(ise)
}


#' Calculate discrete derivatives
#'
#' @param f_x a vector of f(x)
#' @param x a vector of x
#'
#' @return a vector containing derivatives of \code{f_x}
#'
#' @export
get_deriv <- function(f_x, x) {
    f_prime <- pracma::gradient(f_x, x)

    return(f_prime)
}


#' Check whether vector is convex
#'
#' @param f_x a vector of f(x)
#' @param x a vector of x
#'
#' @return a logical vector indicating \code{f_x} is convex on each \code{x}
#'
#' @export
is.convex <- function(f_x, x) {
    first_deriv <- get_deriv(f_x, x)
    second_deriv <- get_deriv(first_deriv, x)

    ind <- stats::median(1:length(x))

    return( second_deriv > 0 )
}


#' Get design points
#'
#' @param Lt a list of vectors containing time points for each curve
#' @param work.grid a vector containing sorted unique grids
#'
#' @return a matrix containing 2 columns (each row is a point which is a pair of observed timepoints)
#'
#' @export
get_design_index <- function(Lt, work.grid) {
    N <- length(work.grid)   # length of unique grid points
    design_mat <- matrix(0, N, N)
    for (t in Lt){
        ind <- which((work.grid <= max(t)) & (work.grid >= min(t)))
        design_mat[ind, ind] <- 1   # multiple indexing for matrix
    }

    # make the design points to (x, y) form
    x <- as.vector(design_mat)
    y <- as.vector(t(design_mat))

    t1 <- rep(1:N, times = N)
    t2 <- rep(1:N, each = N)

    res <- cbind(t1[x != 0],
                 t2[y != 0])

    return(res)
}



#' Rbind the list containing matrix having same number of columns
#'
#' @param x a list containing matrices having same number of columns
#'
#' @return a matrix concatenating all matrices row-wise in a list
#'
#' @export
list2rbind <- function(x) {
    if (!is.list(x)) {
        stop("Input data is not list type.")
    }

    x_is_numeric <- sapply(x, is.numeric)
    if (FALSE %in% x_is_numeric) {
        stop("At least one object in list is not numeric type.")
    }

    # x_dim <- sapply(x, dim)
    # x_dim_uniq <- apply(x_dim, 1, unique)
    # if (length(x_dim_uniq) != 2) {
    #   stop("The dimension of matrix for each list is not same.")
    # }

    M <- length(x)
    A <- x[[1]]
    for (m in 2:M) {
        A <- rbind(A,
                   x[[m]])
    }
    return(A)
}



#' Convert a matrix to a list
#'
#' @param X a n x p matrix containing functional trajectories.
#' n is the number of curves, and p is the number of timepoints.
#' @param grid a vector containing observed timepoints.
#' Default is setting to equal grids between 0 and 1.
#'
#' @return a list containing 2 list (Lt, Ly).
#'
#' @export
matrix2list <- function(X, grid = NULL) {
    n <- nrow(X)

    if (is.null(grid)) {
        grid <- seq(0, 1, length.out = ncol(X))
    }

    Ly <- list()
    Lt <- list()

    for (i in 1:n) {
        t <- grid
        y <- X[i, ]

        NA_ind <- which(is.na(y))
        if (length(NA_ind) > 0) {
            t <- t[-NA_ind]
            y <- y[-NA_ind]
        }
        Ly[[i]] <- as.numeric(y)
        Lt[[i]] <- t
    }

    x.2 <- list(Ly = Ly,
                Lt = Lt)

    return(x.2)
}


#' Convert a list to a matrix
#'
#' @param X a list containing 2 list (Lt, Ly).
#'
#' @return a n x p matrix containing functional trajectories.
#' n is the number of curves, and p is the number of timepoints.
#'
#' @export
list2matrix <- function(X) {
    Lt <- X$Lt
    Ly <- X$Ly
    id <- lapply(1:length(Lt), function(i){ rep(i, length(Lt[[i]]))  })

    # spread data
    x <- data.frame(id = unlist(id),
                    y = unlist(Ly),
                    t = unlist(Lt)) %>%
        tidyr::spread(key = "t", value = "y")
    x <- x[, -1] %>%
        as.matrix

    return(x)
}


#' Plot functional trajectories
#'
#' @param x a list containing 2 list (Lt, Ly) or n x p matrix.
#' @param Lt a list containing obeserved time grids of each trajectory. (If x is NULL)
#' @param Ly a list containing obeserved values of each trajectory. (If x is NULL)
#' @param ... Options of \code{plot()}.
#'
#' @export
plot_curve <- function(x = NULL, Lt = NULL, Ly = NULL, ...) {
    if (is.null(x)) {
        if (!(is.list(Lt) & is.list(Ly))) {
            stop("Lt and Ly should be only list.")
        }

        Lt_len <- sapply(Lt, length)
        Ly_len <- sapply(Ly, length)

        if (!all.equal(Lt_len, Ly_len)) {
            stop("Length of each trajectory of Lt and Ly are not same.")
        }

        x <- list(Lt = Lt,
                  Ly = Ly)
        x <- list2matrix(x)
    }

    if (is.list(x)) {
        gr <- sort(unique(unlist(Lt)))
        graphics::matplot(gr, t(x), ...)
    } else if (is.matrix(x)) {
        graphics::matplot(t(x), ...)
    }
}
