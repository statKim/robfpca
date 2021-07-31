#include <Rcpp.h>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// smoothing spline using a R function "stats::smooth.spline"
// [[Rcpp::export]]
Rcpp::NumericVector smooth_spline_cpp(Rcpp::NumericVector x) {
    // Obtain environment containing function
    Rcpp::Environment base("package:stats");

    // Make function callable from C++
    Rcpp::Function sm_spline_r = base["smooth.spline"];

    // Call the function and receive its list output
    Rcpp::List res = sm_spline_r(Rcpp::_["x"] = x);

    // Return test object in list structure
    return res["y"];
}

// Huber's location M-estimator
// [[Rcpp::export]]
double huber_cpp(NumericVector y,
                 double k = 1.5,
                 double tol = 1e-6) {
    y = y[!is_na(y)];
    int n = y.size();
    double mu = median(y);   // median

    if (n == 0) {
        Rcpp::stop("There is no available observation to compute Huber's location M-estimator.");
    }

    // MAD
    NumericVector s_ = abs(y - mu);
    double s = 1.4826 * median(s_);
    if (s == 0) {
        Rcpp::stop("Estimated MAD is 0 in Huber's location M-estimator.");
        // Rcout << "Estimated MAD is 0." << "\n";
        return 0;
    }

    // Obtain Huber's M-estimator
    NumericVector yy;
    double mu1;
    while(1) {   // infinite loop
        yy = pmin(pmax(mu - k * s, y), mu + k * s);
        mu1 = sum(yy) / n;

        if (abs(mu - mu1) < tol*s) {
            break;
        }

        mu = mu1;
    }

    return mu;
}


// Obtain column-wise M-estimator
// [[Rcpp::export]]
NumericVector mean_Mest_cpp(NumericMatrix X,
                            bool smooth = false) {
    int p = X.ncol();
    NumericVector mu(p);
    for (int i = 0; i < p; i++) {
        mu[i] = huber_cpp(X.column(i));
    }

    // smoothing spline
    if (smooth == true) {
        mu = smooth_spline_cpp(mu);
    }

    return mu;
}


// Obtain marginal M-estimator for covariance
// [[Rcpp::export]]
NumericMatrix cov_Mest_cpp(NumericMatrix X,
                           bool smooth = false) {
    int n = X.nrow();
    int p = X.ncol();
    NumericMatrix rob_var(p, p);

    LogicalVector NA_ind(n);
    IntegerVector ind_full = seq(0, n-1);
    IntegerVector ind_st;
    IntegerVector ind;
    NumericVector X_i;
    NumericVector mu;
    NumericMatrix X_sub;
    NumericVector A;

    if (sum(is_na(X)) == 0) {
        // Complete curves
        mu = mean_Mest_cpp(X, smooth);
        for (int i = 0; i < X.nrow(); i++) {
            X(i, _) = X(i, _) - mu;   // X-mu
        }
        A = NumericVector(n);
        for (int s = 0; s < p; s++) {
            for (int t = 0; t < p; t++) {
                if (s <= t) {
                    for (int i = 0; i < X.nrow(); i++) {
                        A[i] = X(i, s) * X(i, t);
                    }

                    rob_var(s, t) = huber_cpp(A);
                } else {
                    rob_var(s, t) = rob_var(t, s);
                }
            }
        }
    } else {
        // Partially observed curves having missing
        for (int s = 0; s < p; s++) {
            for (int t = 0; t < p; t++) {
                if (s <= t) {
                    // Find the index of curve which have NA at time s or t
                    ind_st = IntegerVector::create(s, t);
                    for (int i = 0; i < n; i++) {
                        X_i = X.row(i);
                        X_i = X_i[ind_st];
                        NA_ind[i] = sum(is_na(X_i)) > 0;
                    }

                    if (sum(NA_ind) > 0) {
                        ind = ind_full[!NA_ind];
                    } else {
                        ind = ind_full;
                    }

                    X_sub = NumericMatrix(ind.size(), p);
                    for (int i = 0; i < ind.size(); i++) {
                        X_sub(i, _) = X(ind[i], _);
                    }

                    mu = mean_Mest_cpp(X_sub, smooth);

                    A = NumericVector(ind.size());
                    for (int i = 0; i < X_sub.nrow(); i++) {
                        X_sub(i, _) = X_sub(i, _) - mu;
                        A[i] = X_sub(i, s) * X_sub(i, t);
                    }

                    rob_var(s, t) = huber_cpp(A);
                } else {
                    rob_var(s, t) = rob_var(t, s);
                }
            }
        }
    }

    // Rcout << rob_var << "\n";

    return rob_var;
}


//
// /*** R
// # x <- matrix(rnorm(10000), nrow = 10)
// system.time({
//   cov1 <- cov_Mest_cpp(x)
// })
// system.time({
//   cov2 <- var.rob.missfd(x, make.pos.semidef=F)
// })
// all.equal(cov1, cov2)
//
// # y <- matrix(rnorm(10000), nrow = 10)
// # system.time({
// #   mu1 <- mean_Mest_cpp(x)
// # })
// #
// # system.time({
// #   mu2 <- as.numeric(mean.rob.missfd(x))
// # })
// #
// # all.equal(mu1, mu2)
// # # # mu <- apply(x, 2, function(t){ MASS::huber(t)$mu })
// */


///////////////////////////////////
/// Local M-estimator
//////////////////////////////////

// expand.grid
// [[Rcpp::export]]
NumericMatrix expand_grid_cpp(NumericVector x,
                              NumericVector y) {
    int n = x.length();

    NumericMatrix res(pow(n, 2), 2);
    res.column(0) = rep_len(x, pow(n, 2));
    res.column(1) = rep_each(y, n);

    return res;
}


// 1-dimensional Local M-estimator
// [[Rcpp::export]]
NumericVector local_M_1D(NumericVector x,
                         NumericVector t,
                         NumericVector new_t,
                         double h) {
  int p = new_t.length();
  NumericVector mu(p);
  for (int i = 0; i < p; i++) {
    LogicalVector i_neighbor = abs(t - new_t[i]) < h;
    NumericVector X_i_neighbor = x[i_neighbor];
    mu[i] = huber_cpp(X_i_neighbor);
  }

  return mu;
}



// for loop of "get_raw_cov" in R function
// [[Rcpp::export]]
NumericMatrix get_raw_cov_cpp(NumericMatrix X,
                              NumericVector mu,
                              NumericVector gr,
                              bool diag = false) {
    int n = X.rows();
    int p = X.cols();
    // IntegerVector ind_vec = seq(0, p-1);
    LogicalVector ind(p);

    NumericMatrix raw_cov;
    for (int i = 0; i < n; i++) {
        // IntegerVector ind = ind_vec[!is_na(X.row(i))];
        ind = !is_na(X.row(i));
        NumericVector gr_i = gr[ind];
        NumericVector X_i = X.row(i);
        NumericVector mu_i = mu[ind];
        X_i = X_i[ind];
        X_i = X_i - mu_i;

        NumericMatrix exp_grid = expand_grid_cpp(gr_i, gr_i);
        NumericMatrix exp_X_i = expand_grid_cpp(X_i, X_i);
        NumericVector exp_grid_1 = exp_grid.column(0);
        NumericVector exp_grid_2 = exp_grid.column(1);

        if (diag == false) {
            // eliminate diagonal parts
            IntegerVector idx = seq(0, exp_X_i.nrow()-1);
            idx = idx[exp_grid_1 != exp_grid_2];
            NumericVector cov_i = exp_X_i.column(0) * exp_X_i.column(1);
            cov_i = cov_i[idx];
            NumericVector exp_grid_1_sub = exp_grid_1[idx];
            NumericVector exp_grid_2_sub = exp_grid_2[idx];

            // We use "cbind" in Rcpp, so we should transpose it.
            NumericMatrix raw_cov_i(4, cov_i.length());
            raw_cov_i(0, _) = cov_i;
            raw_cov_i(1, _) = exp_grid_1_sub;
            raw_cov_i(2, _) = exp_grid_2_sub;
            raw_cov_i(3, _) = rep(i+1, cov_i.length());   // index of curve

            if (i == 0) {
                raw_cov = raw_cov_i;
            } else {
                raw_cov = cbind(raw_cov, raw_cov_i);
            }
        } else {
            NumericVector cov_i = exp_X_i.column(0) * exp_X_i.column(1);

            // We use "cbind" in Rcpp, so we should transpose it.
            NumericMatrix raw_cov_i(4, cov_i.length());
            raw_cov_i(0, _) = cov_i;
            raw_cov_i(1, _) = exp_grid_1;
            raw_cov_i(2, _) = exp_grid_2;
            raw_cov_i(3, _) = rep(i, cov_i.length());

            if (i == 0) {
                raw_cov = raw_cov_i;
            } else {
                raw_cov = cbind(raw_cov, raw_cov_i);
            }
        }
    }

    // transpose because we cbind the raw_cov_i
    raw_cov = transpose(raw_cov);

    return raw_cov;
}


// 2-dimensional Local M-estimator
// for loop of "cov_local_M" in R function
// [[Rcpp::export]]
NumericMatrix local_M_2D(NumericVector raw_cov,
                         NumericVector s,
                         NumericVector t,
                         NumericVector gr,
                         double h = 0.02) {
    int p = gr.length();
    LogicalVector i_neighbor(s.length());
    LogicalVector j_neighbor(s.length());
    LogicalVector ind(s.length());
    NumericMatrix cov_mat(p, p);
    for (int i = 0; i < p; i++) {
        i_neighbor = abs(s - gr[i]) < h;

        for (int j = 0; j < p; j++) {
            if (i <= j) {
                j_neighbor = abs(t - gr[j]) < h;
                ind = i_neighbor * j_neighbor;
                NumericVector raw_cov_sub = raw_cov[ind];

                cov_mat(i, j) = huber_cpp(raw_cov_sub);
            } else {
                cov_mat(i, j) = cov_mat(j, i);
            }
        }
    }

    return cov_mat;
}




//
// /*** R
// x <- matrix(c(1,NA,NA,10,2,3,4,NA), nrow = 2)
// get_raw_cov_cpp(x, rep(0, 4), seq(0, 1, length.out = 4), FALSE)
// # set.seed(100)
// # x <- matrix(rnorm(510), ncol = 51)
// # gr <- seq(0, 1, length.out = 51)
// # x.2 <- sim_delaigle(n = ,
// #                     model = 2,
// #                     type = data_type,
// #                     out.prop = out_prop,
// #                     out.type = out_type)
// # x <- list2matrix(x.2)
//
// # system.time({
// #   raw_cov <- get_raw_cov_cpp(x,
// #                              rep(0, 51),
// #                              gr, TRUE)
// #   a1 <- cov_local_M_cpp(raw_cov[, 1],
// #                         raw_cov[, 2],
// #                         raw_cov[, 3],
// #                         gr,
// #                         0.1)
// # })
// # system.time({
// #   a2 <- test(x, 0.1, gr)
// # })
//
// # all.equal(a1, a2)
//
// # GA::persp3D(gr, gr, a1,
// #             theta = -70, phi = 30, expand = 1)
// # GA::persp3D(gr, gr, a2,
// #             theta = -70, phi = 30, expand = 1)
// # system.time({
// #   cov_mat <- cov_local_M_cpp(raw_cov = raw_cov,
// #                              s = s,
// #                              t = t,
// #                              gr = gr,
// #                              h = 0.05)
// #
// # })
//
//
// */
