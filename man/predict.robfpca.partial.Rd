% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robfpca.partial.R
\name{predict.robfpca.partial}
\alias{predict.robfpca.partial}
\title{Predict FPC scores, reconstruction and completion for a new data}
\usage{
\method{predict}{robfpca.partial}(
  object,
  type = c("score", "reconstr", "comp"),
  newdata = NULL,
  K = NULL,
  ...
)
}
\arguments{
\item{object}{a \code{robfpca.partial} object from \code{robfpca.partial()}}

\item{type}{"score" gives FPC scores, "reconstr" gives reconstruction of each curves, and "comp" gives completion of each curves.}

\item{newdata}{a n x p matrix containing n curves observed at p timepoints}

\item{K}{a number of FPCs}

\item{...}{Not used}
}
\description{
Predict FPC scores, reconstruction and completion for a new data
}
\examples{
### Generate example data
set.seed(100)
x.list <- sim_delaigle(n = 100,
                       type = "partial",
                       out.prop = 0.2,
                       out.type = 1,
                       dist = "normal")
x <- list2matrix(x.list)
matplot(t(x), type = "l")

### Robust FPCA for partially observed functional data
### Given bandwidth = 0.1
fpca.obj <- robfpca.partial(x,
                            type = "huber",
                            PVE = 0.95,
                            bw = 0.1)
fpc.score <- fpca.obj$pc.score

new_data <- x[1:5, ]   # example of new data

### Predict FPC score
pred_score <- predict(fpca.obj, type = "score", newdata = new_data)
pred_score

### Reconstruction
pred_reconstr <- predict(fpca.obj, type = "reconstr", newdata = new_data)
pred_reconstr
par(mfrow = c(1, 2))
matplot(t(new_data), type = "l",
        xlab = "t", ylab = "", main = "Observed curves")
matplot(t(pred_reconstr), type = "l",
        xlab = "t", ylab = "", main = "Reconstructed curves")

### Completion
pred_comp <- predict(fpca.obj, type = "comp", newdata = new_data)
pred_comp
matplot(t(new_data), type = "l",
        xlab = "t", ylab = "", main = "Completion")
matlines(t(pred_comp), type = "l", lty = 1, lwd = 2)

}
