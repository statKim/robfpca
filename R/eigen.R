
#' Get eigen analysis results
#'
#' @param cov a matrix of covariance matrix
#' @param grid a vector containing observed time points. (same with dimension of \code{cov})
#'
#' @export
get_eigen <- function(cov, grid) {
    eig <- eigen(cov, symmetric = T)
    positive_ind <- which(eig$values > 0)
    lambda <- eig$values[positive_ind]
    phi <- eig$vectors[, positive_ind]
    PVE <- cumsum(lambda) / sum(lambda)

    # normalization since discretization technique - from "fdapace" package
    lambda <- lambda * (grid[2] - grid[1])
    phi <- apply(phi, 2, function(x) {
        x <- x / sqrt(fdapace::trapzRcpp(grid, x^2))
        if ( 0 <= sum(x * 1:length(grid)) )
            return(x)
        else
            return(-x)
    })

    return(list(lambda = lambda,
                phi = phi,
                PVE = PVE))
}



#' Change the sign of eigenvectors to target eigenvectors
#'
#' @param eig_vec a vector (n x 1) or matrix (n x q) containing estimated eigenvectors
#' @param target a vector (n x 1) or matrix (n x q) containing true (or want to convert) eigenvectors
#'
#' @export
check_eigen_sign <- function(eig_vec, target) {
    if (is.matrix(eig_vec) & is.matrix(target)) {
        ## if inputs are eigenvector matrices
        eig_vec_dim <- dim(eig_vec)
        target_dim <- dim(target)
        if (!isTRUE(all.equal(eig_vec_dim, target_dim))) {
            stop("The dimensions of 2 eigenvector matrices are not equal.")
        }

        for (i in 1:eig_vec_dim[2]) {
            sse_pos <- sum((eig_vec[, i] - target[, i])^2)
            sse_neg <- sum(((-eig_vec[, i]) - target[, i])^2)

            if (sse_pos > sse_neg) {
                eig_vec[, i] <- -eig_vec[, i]
            }
        }
    } else if (is.numeric(eig_vec) & is.numeric(target)) {
        ## if inputs are eigenvectors
        if (length(eig_vec) != length(target)) {
            stop("The dimensions of 2 eigenvector matrices are not equal.")
        }

        sse_pos <- sum((eig_vec - target)^2)
        sse_neg <- sum(((-eig_vec) - target)^2)

        if (sse_pos > sse_neg) {
            eig_vec <- -eig_vec
        }
    } else {
        stop("Inputs are not n x q matrices or q-dim vectors. Check the input objects.")
    }

    return(eig_vec)
}
