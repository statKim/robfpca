#include <RcppEigen.h>
#include <Rcpp.h>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

// [[Rcpp::depends(RcppEigen)]]


typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;

using Rcpp::List;
using Rcpp::as;

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision

using std::exp;
using std::isnan;


// Get data with positive kernel weights for local regression
// - Y : observed data (unlist(Ly))
// - X : observed grid (unlist(Lt))
// - W : weight from kernel function
// [[Rcpp::export]]
Rcpp::List get_positive_elements(Eigen::VectorXd Y,
                                 Eigen::MatrixXd X,
                                 Eigen::VectorXd W) {
  int n = W.size();

  Rcpp::NumericVector Y_pos(n);
  Rcpp::NumericMatrix X_pos(n, X.cols());
  Rcpp::NumericVector W_pos(n);
  Rcpp::NumericVector row_na(X.cols(), R_NaN);
  Rcpp::NumericVector row_x(X.cols());
  for (int i = 0; i < n; i++) {
    if (W(i) > 0) {
      Y_pos(i) = Y(i);
      row_x = Rcpp::wrap(X.row(i));
      X_pos.row(i) = row_x;
      W_pos(i) = W(i);
    } else {
      Y_pos(i) = R_NaN;
      X_pos.row(i) = row_na;
      W_pos(i) = R_NaN;
    }
  }

  Rcpp::NumericVector Y_rm_neg = na_omit(Y_pos);
  Rcpp::NumericVector W_rm_neg = na_omit(W_pos);

  int n_pos = Y_rm_neg.length();

  Rcpp::NumericMatrix X_rm_neg(n_pos, X.cols());
  Rcpp::NumericVector col_x(n_pos);
  for (int j = 0; j < X.cols(); j++) {
    col_x = na_omit(X_pos.column(j));
    X_rm_neg.column(j) = col_x;
  }

  return Rcpp::List::create(Rcpp::_["Y"] = Y_rm_neg,
                            Rcpp::_["X"] = X_rm_neg,
                            Rcpp::_["W"] = W_rm_neg,
                            Rcpp::_["n_pos"] = n_pos);
}


// Get psi function of robust loss function
// [[Rcpp::export]]
Eigen::VectorXd get_psi(Eigen::VectorXd tmp,
                        int n,
                        double k,
                        std::string method = "HUBER") {
  Eigen::VectorXd psi(n);

  if (method == "HUBER") {
    // psi function of Huber loss
    for (int j = 0; j < n; j++) {
      if (abs(tmp(j)) >= k) {
        if (tmp(j) < 0) {
          psi(j) = -1 * k;
        } else {
          psi(j) = k;
        }
      } else {
        psi(j) = tmp(j);
      }
      psi(j) = 2*psi(j);   // exact derivative
    }
  } else if (method == "BISQUARE") {
    // psi function of Tukey's Biweight loss
    for (int j = 0; j < n; j++) {
      if (abs(tmp(j)) >= k) {
        psi(j) = 0;
      } else {
        psi(j) = pow(1 - pow((tmp(j) / k), 2), 2)*tmp(j);
      }
      psi(j) = 6*psi(j) / pow(k, 2);   // exact derivative
    }
  }

  return psi;
}



// Get M-scale estimator using Bisquare
// [[Rcpp::export]]
double scale_M(Rcpp::NumericVector resid,
               double s,
               std::string method = "BISQUARE",
               const int maxit = 50,
               const double tol = 0.0001) {
  int n = resid.size();
  Rcpp::NumericVector W(n);
  int iter;
  Rcpp::NumericVector s_M(maxit);
  s_M[0] = s;

  if (method == "HUBER") {

  } else if (method == "BISQUARE") {
    double delta = 0.5;
    double c = 1.56;
    // iterate until scale M-estimator is converged
    for (int i = 0; i < maxit; i++) {
      iter = i+1;
      // Weight using rho function of Tukey's Biweight with k = 1
      Rcpp::NumericVector resid_s = resid / s_M[i];
      Rcpp::NumericVector rho = pmin(1,
                                     3*pow(resid_s/c, 2)-3*pow(resid_s/c, 4) + pow(resid_s/c, 6));
      // W = pmin(1/pow(resid_s, 2),
      //          3/pow(c, 2)-3*pow(resid_s, 2)/pow(c, 4) + pow(resid_s, 4)/pow(c, 6));
      W = rho / pow(resid_s, 2);
      // for (int j = 0; j < n; j++) {
      //   // compute weight
      //   if (resid_s(j) == 0) {
      //     W[j] = 1;
      //   } else {
      //     // rho function of Tukey's Biweight loss
      //     if (abs(resid_s(j)) > k) {
      //       // rho(j) = 1;
      //       rho = 1;
      //     } else {
      //       // rho(j) = 1 - pow(1 - pow((resid_s(j) / k), 2), 3);
      //       rho = 1 - pow(1 - pow((resid_s(j) / k), 2), 3);
      //     }
      //     W[j] = rho / pow(resid_s(j), 2);
      //   }
      // }

      s_M[i+1] = sqrt(mean(W * pow(resid, 2)) / delta);

      // Check convergence
      if ((iter > 2) & (abs(s_M[i+1]/s_M[i] - 1) < tol)) {
        break;
      }

      // If fail to convege, return warning message.
      if (iter == maxit) {
        Rcpp::warning("M-scale is not converged.");
        // Rcpp::Rcout << "Converge fail in the IRLS algorithm for maxit=" << maxit;
      }
    }
  } else {
    Rcpp::stop("Not supported option.");
  }

  // double s_M = sqrt(sum(W * pow(resid, 2)) / (n*delta));

  return s_M[iter];
}



// Iteratively re-weighted least squares (IRLS) for robust regression (M-estimation)
// - Y : observed data (unlist(Ly))
// - X : observed grid (unlist(Lt))
// - weight_ : weight vector (Default = NULL)
// - method : "HUBER" or "BISQUARE" are supported
// - maxit : maximun iteration numbers
// - weight : additional weight (for kernel regression)
// - tol : tolerence rate
// - k : cut-off value for Huber or Tukey's biweight function
// [[Rcpp::export]]
Rcpp::List IRLScpp(const Eigen::VectorXd Y,
                   const Eigen::MatrixXd X,
                   Rcpp::Nullable<Rcpp::NumericVector> weight_ = R_NilValue,
                   std::string method = "HUBER",
                   const int maxit = 50,
                   const double tol = 0.0001,
                   const double k = 1.345) {
  int n = Y.size();   // number of observations (unlist(Lt))

  // initial value for beta (LSE estimator)
  Eigen::MatrixXd beta(maxit+1, X.cols());   // container of beta for iterations
  beta.setZero();   // initialize to 0
  Eigen::LDLT<Eigen::MatrixXd> ldlt_XTX(X.transpose() * X);
  beta.row(0) = ldlt_XTX.solve(X.transpose() * Y);   // LSE estimator

  // initial scale estimator
  Eigen::VectorXd Y_hat(n);   // Y_hat
  Eigen::VectorXd beta_hat(X.cols());   // beta_hat
  beta_hat = beta.row(0);
  Y_hat = X * beta_hat;   // prediction
  Rcpp::NumericVector resid = Rcpp::wrap(Y - Y_hat);   // residuals
  Rcpp::NumericVector resid_med = abs(resid - median(resid));   // MAD
  double s = median(resid_med) * 1.4826;  // re-scaled MAD by MAD*1.4826

  // if weight is NULL, set 1
  Eigen::VectorXd weight(n);
  if (weight_.isNotNull()) {
    weight = Rcpp::as<Eigen::VectorXd>(weight_);
  } else {
    weight.setOnes();
  }

  // variables for iterations
  Eigen::VectorXd tmp(n);
  Eigen::VectorXd psi(n);
  Eigen::VectorXd w(n);
  Eigen::MatrixXd XTWX(X.cols(), X.cols());   // X'WX
  Eigen::LDLT<Eigen::MatrixXd> ldlt_XTWX(XTWX);   // LDLT obj for X'WX
  int iter = 0;

  // iterate until beta converged
  for (int i = 0; i < maxit; i++) {
    iter = i+1;
    // beta_hat = beta.row(i);
    // Y_hat = X * beta_hat;

    // // update scale estimator
    // resid = Rcpp::wrap(Y - Y_hat);
    // resid_med = abs(resid - median(resid));
    // s = median(resid_med) * 1.4826;  // re-scaled MAD by MAD*1.4826

    tmp = (Y - Y_hat).array() / s;

    // psi function of robust loss function
    psi = get_psi(tmp, n, k, method);

    // weight matrix for WLS
    w = (weight.array() * psi.array()).array() / (Y - Y_hat).array();

    // Weighted Least Squares using current estimates
    XTWX = X.transpose() * w.asDiagonal() * X;
    ldlt_XTWX.compute(XTWX);
    beta.row(i+1) = ldlt_XTWX.solve(X.transpose() * w.asDiagonal() * Y);
    beta_hat = beta.row(i+1);

    // Check convergence (See Marronna(2006) page 105)
    Rcpp::NumericVector resid_before = resid;   // residuals of before iteration
    Y_hat = X * beta_hat;
    resid = Rcpp::wrap(Y - Y_hat);   // current residuals
    s = scale_M(resid, s, "BISQUARE");   // update M-scale estimator
    if ((iter > 2) & (max(abs(resid_before - resid))/s < tol)) {
    // if ((beta.row(i+1) - beta.row(i)).cwiseAbs().sum() < tol) {
      break;
    }

    // If fail to convege, return warning message.
    if (iter == maxit) {
      Rcpp::warning("IRLS algorithm is failed to converge for maxit = %i iterations.", maxit);
      // Rcpp::Rcout << "Converge fail in the IRLS algorithm for maxit=" << maxit;
    }
  }

  return Rcpp::List::create(Rcpp::_["beta"] = beta_hat,
                            Rcpp::_["iter"] = iter,
                            Rcpp::_["scale.est"] = s);
}


// Get kernel weights for local smoothing
// [[Rcpp::export]]
Eigen::VectorXd get_kernel_weight(Eigen::VectorXd tmp,
                                  std::string kernel = "epanechnikov") {
  int n = tmp.size();
  Eigen::VectorXd kern(n);
  if (kernel == "epanechnikov") {   // Epanechnikov kernel
    kern = (3./4.) * (1 - tmp.array().pow(2));
  } else if (kernel == "gauss") {   // Gaussian kernel
    kern = (1./sqrt(2.*M_PI)) * (-0.5 * tmp.array().pow(2)).exp();
  }

  return kern;
}


// Local polynomial kernel smoothing with huber loss (mean estimator)
// It is used for local_kern_smooth() in R function
// - Lt : a list of vectors or a vector containing time points for all curves
// - Ly : a list of vectors or a vector containing observations for all curves
// - newt : a vector containing time points to estimate
// - method : "HUBER" or "BISQUARE" are supported
// - kernel : a kernel function for kernel smoothing ("epan", "gauss" are supported.)
// - bw : bandwidth
// - k : cut-off value for Huber or Tukey's biweight function
// - deg : degree of polynomial
// - maxit : maximum iteration for IRLS algorithm
// - tol : tolerance rate for stop iterations
// [[Rcpp::export]]
Eigen::VectorXd locpolysmooth_cpp(Eigen::VectorXd Lt,
                                  Eigen::VectorXd Ly,
                                  Eigen::VectorXd newt,
                                  std::string method = "HUBER",
                                  std::string kernel = "epanechnikov",
                                  const double bw = 0.1,
                                  const double k = 1.345,
                                  const int deg = 1,
                                  int maxit = 50,
                                  const double tol = 0.0001) {
  int n_newt = newt.size();   // number of grid which is predicted
  int n = Lt.size();   // number of Lt
  double weig = 1. / n;   // 1/length(Lt)


  Eigen::VectorXd beta_hat(deg+1);
  Eigen::VectorXd mu_hat(n_newt);
  for (int t = 0; t < n_newt; t++) {
    Eigen::VectorXd tmp(n);
    tmp = (Lt.array() - newt(t)) / bw;

    // Kernel weights
    Eigen::VectorXd kern(n);
    kern = get_kernel_weight(tmp, kernel);
    // if (kernel == "epanechnikov") {   // Epanechnikov kernel
    //   kern = (3./4.) * (1 - tmp.array().pow(2));
    // } else if (kernel == "gauss") {   // Gaussian kernel
    //   kern = (1./sqrt(2.*M_PI)) * (-0.5 * tmp.array().pow(2)).exp();
    // }

    // X matrix
    Eigen::MatrixXd X(n, deg+1);
    X.setOnes();
    for (int d = 0; d < deg; d++) {
      X.col(d+1) = (Lt.array() - newt(t)).pow(d+1);
    }

    // obtain data which has non-negative weights
    Rcpp::List pos_idx_obj = get_positive_elements(Ly, X, kern);
    int n_pos = pos_idx_obj["n_pos"];
    Eigen::MatrixXd X_sub = pos_idx_obj["X"];
    Eigen::VectorXd Y_sub = pos_idx_obj["Y"];
    Eigen::VectorXd kern_sub = pos_idx_obj["W"];

    // weight vector for w_i * K_h(x)
    Rcpp::NumericVector W(n_pos);
    W = (weig * kern_sub.array()) / bw;

    // Rcpp::Rcout << W << "\n";
    // Rcpp::Rcout << n << "\t" << W.size() << "\t" << kern.size() << "\n";

    // M-type Robust regression
    Rcpp::List fit = IRLScpp(Y_sub, X_sub, W, method, maxit, tol, k);
    beta_hat = fit["beta"];
    mu_hat(t) = beta_hat(0);
  }

  return mu_hat;
}



//
// /*** R
// library(mvtnorm)
// library(robfpca)
// source("R/sim_delaigle.R")
// # source("R/sim_Lin_Wang(2020).R")
// # source("load_source.R")
//
// num_sim <- 30   # number of simulations
// out_prop <- 0   # proportion of outliers
// data_type <- "snippet"   # type of functional data
// # kernel <- "epanechnikov"   # kernel function for local smoothing
// kernel <- "gauss"   # kernel function for local smoothing
// bw <- 0.1   # fixed bandwidth
// n_cores <- 12   # number of threads for parallel computing
//
// set.seed(1)
// n <- 100
// n.grid <- 51
// # x.2 <- sim_kraus(n = n, out.prop = 0.2, out.type = 1, grid.length = n.grid)
// x.2 <- sim_delaigle(n = n,
//                     model = 2,
//                     type = data_type,
//                     out.prop = out_prop,
//                     out.type = 1)
// df <- data.frame(
//   id = factor(unlist(sapply(1:length(x.2$Lt),
//                             function(id) {
//                               rep(id, length(x.2$Lt[[id]]))
//                             })
//   )),
//   y = unlist(x.2$Ly),
//   t = unlist(x.2$Lt)
// )
//
// work.grid <- seq(0, 1, length.out = n.grid)
//
//
//
// Lt <- x.2$Lt
// Ly <- x.2$Ly
// method <- "HUBER"
// deg <- 1
// K <- 5
// bw_cand = 10^seq(-1, 0, length.out = 5)/3
// cv_loss<-"HUBER"
// delta <- 1.345
//
//
// # get index for each folds
// folds <- list()
// n <- length(Lt)   # the number of curves
// fold_num <- n %/% K   # the number of curves for each folds
// fold_sort <- sample(1:n, n)
// for (k in 1:K) {
//   ind <- (fold_num*(k-1)+1):(fold_num*k)
//   if (k == K) {
//     ind <- (fold_num*(k-1)+1):n
//   }
//   folds[[k]] <- fold_sort[ind]
// }
// k <- 1
//
//
// for (k in 1:K) {
//   print(k)
//   Lt_train <- Lt[ -folds[[k]] ]
//   Ly_train <- Ly[ -folds[[k]] ]
//   Lt_test <- Lt[ folds[[k]] ]
//   Ly_test <- Ly[ folds[[k]] ]
//
//   for (i in 1:length(bw_cand)) {
//     # print(i)
//     # y_hat <- local_kern_smooth(Lt = Lt_train,
//     #                            Ly = Ly_train,
//     #                            newt = Lt_test,
//     #                            method = method,
//     #                            bw = bw_cand[i],
//     #                            kernel = kernel,
//     #                            delta = 1.345)
//
//     mu_hat <- locpolysmooth(Lt = unlist(Lt_train),
//                             Ly = unlist(Ly_train),
//                             newt = work.grid,
//                             kernel = kernel,
//                             bw = bw_cand[i],
//                             k = delta,
//                             deg = deg)
//   }
// }
// mu_hat
// */
