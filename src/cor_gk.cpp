#include <Rcpp.h>

using namespace Rcpp;

// MASS::fitdistr()
// [[Rcpp::export]]
Rcpp::NumericVector fitdistr_cpp(Rcpp::NumericVector x,
                                 Rcpp::String densfun,
                                 int df = 3) {
    // Obtain environment containing function
    Rcpp::Environment env("package:MASS");

    // Make function callable from C++
    Rcpp::Function fitdistr_r = env["fitdistr"];

    // Call the function and receive its list output
    Rcpp::List res = fitdistr_r(Rcpp::_["x"] = x,
                                Rcpp::_["densfun"] = densfun,
                                Rcpp::_["df"] = df);

    // Return test object in list structure
    return res["estimate"];
}


// RobStatTM::locScaleM()
// [[Rcpp::export]]
Rcpp::List locScaleM_cpp(Rcpp::NumericVector x,
                         Rcpp::String psi,
                         bool na_rm = 1) {
    // Obtain environment containing function
    Rcpp::Environment env("package:RobStatTM");

    // Make function callable from C++
    Rcpp::Function locScaleM_r = env["locScaleM"];

    // Call the function and receive its list output
    Rcpp::List res = locScaleM_r(Rcpp::_["x"] = x,
                                 Rcpp::_["psi"] = psi,
                                 Rcpp::_["na.rm"] = na_rm);

    // Return test object in list structure
    return res;
}


// stats::quantile()
// [[Rcpp::export]]
double quantile_cpp(Rcpp::NumericVector x,
                    double probs) {
    // Obtain environment containing function
    Rcpp::Environment env("package:stats");

    // Make function callable from C++
    Rcpp::Function ftn_r = env["quantile"];

    // Call the function and receive its list output
    Rcpp::List res = ftn_r(Rcpp::_["x"] = x,
                           Rcpp::_["probs"] = probs,
                           Rcpp::_["names"] = 0);

    // Return test object in list structure
    return res[0];
}


// Robust scale estimator using Method of moments with kappa (See eq (3.4))
// [[Rcpp::export]]
double get_sigma2_rob_cpp(Rcpp::NumericVector v) {
    // 2*Hampel loss
    // cut-off values are obtained from Sinova(2018)
    v = abs(v);
    Rcpp::NumericVector v2 = v[!is_na(v)];
    double a = median(v2);
    double b = quantile_cpp(v2, 0.75);
    double c = quantile_cpp(v2, 0.85);

    for (int i = 0; i < v2.length(); i++) {
        double v_i = v2[i];
        if (v_i < a) {
            v_i = pow(v_i, 2);
        } else if (v_i < b) {
            v_i = 2*a*v_i - pow(a, 2);
        } else if (v_i < c) {
            v_i = a * pow(v_i-c, 2)/(b-c) + a*(b+c-a);
        } else {
            v_i = a*(b+c-a);
        }
        v2[i] = v_i;
    }

    double res = mean(v2);

    return res;
}


// Gnanadesikan-Kettenring (GK) correlation estimation
// [[Rcpp::export]]
Rcpp::List cor_gk_cpp(Rcpp::NumericMatrix X,
                      Rcpp::String type = "huber",
                      bool MM = 1,
                      int df = 3) {
    int p = X.cols();

    // Corss-sectional robust location estimator
    Rcpp::NumericVector rob_mean(p);
    Rcpp::NumericVector rob_disp(p);
    for (int i = 0; i < p; i++) {
        Rcpp::NumericVector X_t = X.column(i);
        X_t = X_t[!is_na(X_t)];
        if (type == "tdist") {
            // "tdist"
            Rcpp::NumericVector obj = fitdistr_cpp(X_t, "t", df);
            rob_mean[i] = obj[0];
            rob_disp[i] = obj[1];
        } else {
            // "huber" or "bisquare"
            Rcpp::List obj = locScaleM_cpp(X_t, type);
            rob_mean[i] = obj["mu"];
            rob_disp[i] = obj["disper"];
        }
    }

    // Compute GK correlation
    Rcpp::NumericMatrix rob_cov(p, p);
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            if (i < j) {
                Rcpp::NumericVector X_i = X.column(i);
                Rcpp::NumericVector X_j = X.column(j);
                Rcpp::LogicalVector ind_not_NA = !is_na(X_i + X_j);

                // Rcout << ind_not_NA << "\n";

                double z1_disp = 0;
                double z2_disp = 0;

                if (type == "tdist") {
                    // "tdist"

                    // Scaling to obtain correlation matrix
                    Rcpp::NumericVector obj_1 = fitdistr_cpp(X_i[ind_not_NA], "t", df);
                    Rcpp::NumericVector obj_2 = fitdistr_cpp(X_j[ind_not_NA], "t", df);

                    Rcpp::NumericVector z1 = (X_i - obj_1[0])/obj_1[1] + (X_j - obj_2[0])/obj_2[1];
                    Rcpp::NumericVector z2 = (X_i - obj_1[0])/obj_1[1] - (X_j - obj_2[0])/obj_2[1];

                    // double disp_1 = fitdistr_cpp(X_i[ind_not_NA], "t", df)[1];
                    // double disp_2 = fitdistr_cpp(X_j[ind_not_NA], "t", df)[1];
                    //
                    // // Rcout << disp_1 << "\t" << disp_2 << "\n";
                    //
                    // Rcpp::NumericVector z1 = X_i/disp_1 + X_j/disp_2;
                    // Rcpp::NumericVector z2 = X_i/disp_1 - X_j/disp_2;

                    // Closed form of robust loss scale estimator using Method of moments
                    // See eq (3.4) in our paper
                    if (MM == 1) {
                        z1_disp = sqrt(get_sigma2_rob_cpp(z1));

                        Rcpp::NumericVector z2_not_NA = z2[ind_not_NA];
                        if (sd(z2_not_NA) >= pow(10, -10)) {
                            z2_disp = sqrt(get_sigma2_rob_cpp(z2));
                        }
                    } else {
                        // Using t-MLE
                        z1_disp = fitdistr_cpp(z1[ind_not_NA], "t", df)[1];

                        Rcpp::NumericVector z2_not_NA = z2[ind_not_NA];
                        if (sd(z2_not_NA) >= pow(10, -10)) {
                            z2_disp = fitdistr_cpp(z2[ind_not_NA], "t", df)[1];
                        }
                    }
                } else {
                    // "huber" or "bisquare"

                    // Scaling to obtain correlation matrix
                    Rcpp::List obj_1 = locScaleM_cpp(X_i[ind_not_NA], type);
                    Rcpp::List obj_2 = locScaleM_cpp(X_j[ind_not_NA], type);

                    // Rcout << disp_1 << "\t" << disp_2 << "\n";

                    Rcpp::NumericVector z1 = (X_i - obj_1["mu"])/obj_1["disper"] + (X_j - obj_2["mu"])/obj_2["disper"];
                    Rcpp::NumericVector z2 = (X_i - obj_1["mu"])/obj_1["disper"] - (X_j - obj_2["mu"])/obj_2["disper"];

                    // Closed form of robust loss scale estimator using Method of moments
                    // See eq (3.4) in our paper
                    if (MM == 1) {
                        z1_disp = sqrt(get_sigma2_rob_cpp(z1));

                        Rcpp::NumericVector z2_not_NA = z2[ind_not_NA];
                        if (sd(z2_not_NA) >= pow(10, -10)) {
                            z2_disp = sqrt(get_sigma2_rob_cpp(z2));
                        }
                    } else {
                        // Using "RobStatTM" package
                        z1_disp = locScaleM_cpp(z1[ind_not_NA], type)["disper"];

                        Rcpp::NumericVector z2_not_NA = z2[ind_not_NA];
                        if (sd(z2_not_NA) >= pow(10, -10)) {
                            z2_disp = locScaleM_cpp(z2[ind_not_NA], type)["disper"];
                        }
                    }
                }

                rob_cov(i, j) = (pow(z1_disp, 2) - pow(z2_disp, 2)) / (pow(z1_disp, 2) + pow(z2_disp, 2));
            } else if (i == j) {
                rob_cov(i, i) = 1;
            } else {
                rob_cov(i, j) = rob_cov(j, i);
            }
        }
    }

    return Rcpp::List::create(_["mean"] = rob_mean,
                              _["cov"] = rob_cov,
                              _["disp"] = rob_disp);
}




//
//
// /*** R
// library(MASS)
// library(RobStatTM)
// # x2 <- rt(250, df = 9)
// #
// # fitdistr(x2, "t", df = 9)
// # fitdistr_cpp(x2, "t", df = 9)
// #
// # x2[c(1, 100, 30, 50)] <- NA
// # locScaleM(x2, psi = "huber", na.rm = T)
// # locScaleM_cpp(x2, psi = "huber", na_rm = T)
// #
// # # robfpca:::get_sigma2_rob(c(1,2,3,NA,NA,6,7,8,NA,10))
// # # get_sigma2_rob_cpp(c(1,2,3,NA,NA,6,7,8,NA,10))
// # robfpca:::get_sigma2_rob(x2)
// # get_sigma2_rob_cpp(x2)
//
// set.seed(1000)
// X <- matrix(rnorm(100), 25, 4)
// type <- "tdist"
// # cor_gk_cpp(X, type = type, MM = T)$cov
// # robfpca:::cov_gk(X, type = type, MM = T, smooth = F, cor = T)$cov
//
// # x2 <- X[, 4]
// # robfpca:::get_sigma2_rob(x2)
// # get_sigma2_rob_cpp(x2)
// library(robfpca)
// x.list <- sim_delaigle(n = 100,
//                        type = "partial",
//                        out.prop = 0.2,
//                        dist = "normal")
// x <- list2matrix(x.list)
// all.equal(
//     cor_gk_cpp(x, type = type, MM = T)$cov,
//     robfpca:::cov_gk(x, type = type, MM = T, smooth = F, cor = T, psd = F)$cov
// )
// */
