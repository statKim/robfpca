##########################################
### Utill functions
##########################################


#' Get function name in the global environment
#'
#' @export
fun2char <- function() {
    env <- ls(envir = .GlobalEnv)
    ind <- sapply(env, function(x) { is.function(get(x)) })
    return(env[ind])
}


#' Calculate Intergrated Squared Errors (ISE)
#'
#' @param x vectors or matrices to compare (y-axis)
#' @param x_hat vectors or matrices to compare (y-axis)
#' @param grid corresponding observed grid (x-axis)
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


#' Calculate raw derivatives
#'
#' @param f_x a vector of f(x)
#' @param x a vector of x
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
#'
#' @export
get_design_index <- function(Lt) {
    obs_grid <- sort(unique(unlist(Lt)))
    N <- length(obs_grid)   # length of unique grid points
    design_mat <- matrix(0, N, N)

    for (t in Lt){
        ind <- which(obs_grid %in% t)
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



