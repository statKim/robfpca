##############################################################
### Robust noise variance estimation for functional snippets
##############################################################

#' Noise variance estimator
#'
#' @param t a list of vectors containing time points for each curve
#' @param y a list of vectors containing observations for each curve
#' @param h a bandwidth
#'
#' @importFrom stats median
#'
#' @export
#' @useDynLib robfpca
sigma2.rob <- function(t, y, h = NULL) {
    if (is.list(y)) {   # irregular data
        n <- length(t)

        t.min <- min(unlist(t))
        t.max <- max(unlist(t))
        if (is.null(h)) {   # Not recommended
            h <- select.sig2.bw(t,y)
        } else if (2*h >= t.max-t.min) {
            stop('h is too large')
        }


        A0A1 <- sapply(1:n, function(i){
            tobs <- t[[i]]
            y <- y[[i]]
            m <- length(tobs)
            A0 <- 0
            A1 <- 0
            B <- 0

            for (j in 1:m) {
                for (k in 1:m) {
                    if ( (k!=j) && (abs(tobs[j]-tobs[k]) < h) ) {
                        A0 <- A0 + (y[j]^2 + y[k]^2)/2
                        A1 <- A1 + y[j]*y[k]
                        B <- B + 1
                    }
                }
            }
            return(c(A0/(m*(m-1)),
                     A1/(m*(m-1)),
                     B/(m*(m-1))))
        })

        # mean((A0A1[1, ] - A0A1[2, ])/A0A1[3, ])
        # mean( mean(A0A1[1, ] - A0A1[2, ]) / mean(A0A1[3, ]) )
        sig2 <- median(A0A1[1, ] - A0A1[2, ]) / median(A0A1[3, ])

    } else {
        stop('unsupported data type')
    }

    return(sig2)
}


### select optimal bandwidth
select.sig2.rob.bw <- function(Lt, Ly, ss = NULL) {
    n <- length(Lt)

    t.min <- min(unlist(Lt))
    t.max <- max(unlist(Lt))

    delta <- max(sapply(Lt, function(ts){ max(ts)-min(ts) }))
    m <- mean(sapply(Lt, length))
    M <- n * m^2
    #h <- delta * (M^(-1/5)) / 7

    # calculate sum of squares - Not recommended (Already it was calculated once!)
    if (is.null(ss)) {
        # mu.hat <- predict(meanfunc(Lt,Ly),Ly)
        mu.hat <- predict(meanfunc.rob(Lt, Ly), Lt)
        ss <- lapply(1:length(Lt),function(i){
            rr <- Ly[[i]] - mu.hat[[i]]
            rr^2
        })
    }

    vn <- sqrt(mean(unlist(ss)))

    h <- 0.29 * delta * vn * (M^(-1/5))

    max.it <- 1000
    it <- 0
    while (it < max.it) {  # h0 two small
        it <- it + 1
        # print(paste(it, ":"))

        cnt <- sapply(1:n, function(i) {
            tobs <- Lt[[i]]
            # y <- Ly[[i]]
            m <- length(tobs)
            v1 <- 0

            if (m < 2) {
                return(0)
            }
            for (j in 1:m) {
                for (k in 1:m) {
                    if ( (k!=j) && (abs(tobs[j]-tobs[k]) < h) ) {
                        v1 <- v1 + 1
                    }
                }
            }
            return(v1)
        })

        cnt <- sum(cnt)
        if (cnt >= min(50, 0.1*n*m*(m-1))) {
            break
        } else {
            h <- h*1.01
        }
    }

    return(h)
}
