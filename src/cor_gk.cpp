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


// stats::mad()
// [[Rcpp::export]]
double mad(Rcpp::NumericVector x) {
    // remove NAs
    x = x[!is_na(x)];

    x = abs(x - median(x));
    double res = 1.4826 * median(x);

    return res;
}


// // Too slow....
// // RobStatTM::rho()
// // [[Rcpp::export]]
// Rcpp::NumericVector rho_cpp(Rcpp::NumericVector u,
//                             Rcpp::String family = "bisquare",
//                             double cc = 1.547645) {
//     // Obtain environment containing function
//     Rcpp::Environment base("package:RobStatTM");
//
//     // Make function callable from C++
//     Rcpp::Function ftn = base["rho"];
//
//     // Call the function and receive its list output
//     Rcpp::NumericVector res = ftn(Rcpp::_["u"] = u,
//                                   Rcpp::_["family"] = family,
//                                   Rcpp::_["cc"] = cc);
//
//     // Return test object in list structure
//     return res;
// }
//
// // RobStatTM::mscale()
// // [[Rcpp::export]]
// double mscale_cpp(Rcpp::NumericVector u,
//                   double delta = 0.5,
//                   double tuning_chi = 1.547645,
//                   Rcpp::String family = "bisquare",
//                   int maxit = 100,
//                   double tol = 1e-6,
//                   double tolerancezero = 2.220446e-16) {
//     u = abs(u);
//     double s0 = median(u) / 0.6745;
//     if (s0 < tolerancezero) {
//         return 0;
//     }
//     double err = tol + 1;
//     int it = 0;
//     double s1;
//     Rcpp::NumericVector sub(u.length());
//     while ((err > tol) & (it < maxit)) {
//         it = it + 1;
//         s1 = mean(rho_cpp(u/s0, family, tuning_chi));
//         s1 = pow(s0, 2) * s1 / delta;
//         s1 = sqrt(s1);
//         err = abs(s1 - s0) / s0;
//         s0 = s1;
//     }
//
//     return s0;
// }


// RobStatTM::mscale()
// [[Rcpp::export]]
Rcpp::NumericVector mscale_cpp(Rcpp::NumericVector u,
                               double delta = 0.5,
                               double tuning_chi = 1.547645,
                               Rcpp::String family = "bisquare") {
    // Obtain environment containing function
    Rcpp::Environment base("package:RobStatTM");

    // Make function callable from C++
    Rcpp::Function ftn = base["mscale"];

    // Call the function and receive its list output
    Rcpp::NumericVector res = ftn(Rcpp::_["u"] = u,
                                  Rcpp::_["delta"] = delta,
                                  Rcpp::_["tuning.chi"] = tuning_chi,
                                  Rcpp::_["family"] = family);

    // Return test object in list structure
    return res;
}



// [[Rcpp::export]]
Rcpp::NumericVector wfun(Rcpp::NumericVector x,
                         Rcpp::String psi) {
    Rcpp::NumericVector w(x.length());
    if (psi == "huber") {
        Rcpp::NumericVector w1(x.length());
        w1 = ifelse(abs(x) <= 1, 1, w1);
        Rcpp::NumericVector w2 = 1 / (abs(x) + 1.e-20);
        w = w1 + ifelse(abs(x) > 1, w2, 0);
    } else if (psi == "bisquare") {
        w = pow((1 - pow(x, 2)), 2);
        w = ifelse(abs(x) <= 1, w, 0);
    }

    return w;
}

// [[Rcpp::export]]
Rcpp::NumericVector psif(Rcpp::NumericVector x,
                         Rcpp::String psi) {
    return x * wfun(x, psi);
}

// [[Rcpp::export]]
Rcpp::NumericVector psipri(Rcpp::NumericVector x,
                           Rcpp::String psi) {
    Rcpp::NumericVector w(x.length());
    if (psi == "huber") {
        w = ifelse(abs(x) <= 1, 1, w);
    } else if (psi == "bisquare") {
        w = pow(1 - pow(x, 2), 2) - 4 * pow(x, 2) * (1 - pow(x, 2));
        w = ifelse(abs(x) < 1, w, 0);
    }

    return w;
}



// RobStatTM::locScaleM()
// [[Rcpp::export]]
Rcpp::List loc_scale_M(Rcpp::NumericVector x,
                       Rcpp::String psi = "huber",
                       double eff = 0.95,
                       int maxit = 50,
                       double tol = 1.e-4) {
    // remove NAs
    x = x[!is_na(x)];
    int n = x.length();

    double k;   // cut-off value of robust loss
    if (psi == "huber") {
        if (eff == 0.95) {
            k = 1.34;
        } else if (eff == 0.90) {
            k = 0.981;
        } else if (eff == 0.85) {
            k = 0.732;
        }
    } else if (psi == "bisquare") {
        if (eff == 0.95) {
            k = 4.685;
        } else if (eff == 0.90) {
            k = 3.88;
        } else if (eff == 0.85) {
            k = 3.44;
        }
    }

    double mu0 = median(x);
    double sig0 = mad(x);

    // intialization
    double dife = 1e10;
    int iter = 0;
    Rcpp::NumericVector resi(n);
    Rcpp::NumericVector ww(n);
    double mu = 0;
    while (dife > tol & iter < maxit) {
        iter = iter + 1;
        resi = (x - mu0) / sig0;
        ww = wfun(resi/k, psi);
        mu = sum(ww * x) / sum(ww);
        dife = abs(mu - mu0) / sig0;
        mu0 = mu;
    }


    Rcpp::NumericVector rek = resi / k;
    Rcpp::NumericVector pp = psif(rek, psi) * k;
    double a = mean(pow(pp, 2));
    double b = mean(psipri(rek, psi));
    double sigmu = pow(sig0, 2) * a / (n * pow(b, 2));
    sigmu = sqrt(sigmu);
    Rcpp::NumericVector scat = mscale_cpp(x - mu, 0.5, 1.56, "bisquare");

    Rcpp::List resu = Rcpp::List::create(_["mu"] = mu,
                                         _["std.mu"] = sigmu,
                                         _["disper"] = scat);

    return resu;
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

    // Rcpp::NumericVector v3(v2.length());
    // v3 = ifelse(v2 < a, pow(v2, 2),
    //             ifelse(v2 < b, 2*a*v2 - pow(a, 2),
    //                    ifelse(v2 < c, a * pow(v2-c, 2)/(b-c) + a*(b+c-a), a*(b+c-a))));
    // double res = mean(v3);
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
        if (type == "tdist") {   // "tdist"
            Rcpp::NumericVector obj = fitdistr_cpp(X_t, "t", df);
            rob_mean[i] = obj[0];
            rob_disp[i] = obj[1];
        } else {   // "huber" or "bisquare"
            Rcpp::List obj = loc_scale_M(X_t, type);
            // Rcpp::List obj = locScaleM_cpp(X_t, type);
            rob_mean[i] = obj["mu"];
            rob_disp[i] = obj["disper"];
        }
    }

    // Compute GK correlation
    Rcpp::NumericMatrix rob_cov(p, p);
    for (int i = 0; i < p; i++) {
        Rcpp::NumericVector X_i = X.column(i);

        for (int j = 0; j < p; j++) {
            if (i < j) {
                Rcpp::NumericVector X_j = X.column(j);
                Rcpp::LogicalVector ind_not_NA = !is_na(X_i + X_j);

                // Rcout << ind_not_NA << "\n";

                double z1_disp = 0;
                double z2_disp = 0;

                if (type == "tdist") {   // "tdist"
                    // Scaling to obtain correlation matrix
                    Rcpp::NumericVector obj_1 = fitdistr_cpp(X_i[ind_not_NA], "t", df);
                    Rcpp::NumericVector obj_2 = fitdistr_cpp(X_j[ind_not_NA], "t", df);

                    Rcpp::NumericVector z1 = (X_i - obj_1[0])/obj_1[1] + (X_j - obj_2[0])/obj_2[1];
                    Rcpp::NumericVector z2 = (X_i - obj_1[0])/obj_1[1] - (X_j - obj_2[0])/obj_2[1];

                    // Closed form of robust loss scale estimator using Method of moments
                    // See eq (3.4) in our paper
                    if (MM == 1) {
                        z1_disp = sqrt(get_sigma2_rob_cpp(z1));

                        Rcpp::NumericVector z2_not_NA = z2[ind_not_NA];
                        if (sd(z2_not_NA) >= 1e-10) {
                            z2_disp = sqrt(get_sigma2_rob_cpp(z2));
                        }
                    } else {
                        // Using t-MLE
                        z1_disp = fitdistr_cpp(z1[ind_not_NA], "t", df)[1];

                        Rcpp::NumericVector z2_not_NA = z2[ind_not_NA];
                        if (sd(z2_not_NA) >= 1e-10) {
                            z2_disp = fitdistr_cpp(z2[ind_not_NA], "t", df)[1];
                        }
                    }
                } else {   // "huber" or "bisquare"
                    // Scaling to obtain correlation matrix
                    Rcpp::List obj_1 = loc_scale_M(X_i[ind_not_NA], type);
                    Rcpp::List obj_2 = loc_scale_M(X_j[ind_not_NA], type);
                    // Rcpp::List obj_1 = locScaleM_cpp(X_i[ind_not_NA], type);
                    // Rcpp::List obj_2 = locScaleM_cpp(X_j[ind_not_NA], type);

                    Rcpp::NumericVector z1 = (X_i - obj_1["mu"])/obj_1["disper"] + (X_j - obj_2["mu"])/obj_2["disper"];
                    Rcpp::NumericVector z2 = (X_i - obj_1["mu"])/obj_1["disper"] - (X_j - obj_2["mu"])/obj_2["disper"];

                    // Closed form of robust loss scale estimator using Method of moments
                    // See eq (3.4) in our paper
                    if (MM == 1) {
                        z1_disp = sqrt(get_sigma2_rob_cpp(z1));

                        Rcpp::NumericVector z2_not_NA = z2[ind_not_NA];
                        if (sd(z2_not_NA) >= 1e-10) {
                            z2_disp = sqrt(get_sigma2_rob_cpp(z2));
                        }
                    } else {
                        // Using "RobStatTM" package
                        z1_disp = loc_scale_M(z1[ind_not_NA], type)["disper"];
                        // z1_disp = locScaleM_cpp(z1[ind_not_NA], type)["disper"];

                        Rcpp::NumericVector z2_not_NA = z2[ind_not_NA];
                        if (sd(z2_not_NA) >= 1e-10) {
                            z2_disp = loc_scale_M(z2[ind_not_NA], type)["disper"];
                            // z2_disp = locScaleM_cpp(z2[ind_not_NA], type)["disper"];
                        }
                    }
                }

                rob_cov(i, j) = (pow(z1_disp, 2) - pow(z2_disp, 2)) / (pow(z1_disp, 2) + pow(z2_disp, 2));
            } else if (i == j) {
                // diagonal parts
                rob_cov(i, i) = 1;
            } else {
                // upper-triangular parts
                rob_cov(i, j) = rob_cov(j, i);
            }
        }
    }

    return Rcpp::List::create(Rcpp::_["mean"] = rob_mean,
                              Rcpp::_["cov"] = rob_cov,
                              Rcpp::_["disp"] = rob_disp);
}





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
// type <- "huber"
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
// # all.equal(
// #     cor_gk_cpp(x, type = type, MM = T)$cov,
// #     robfpca:::cov_gk(x, type = type, MM = T, smooth = F, cor = T, psd = F)$cov
// # )
//
// system.time({
//     cor_gk_cpp(x, type = type, MM = T)
// })
// system.time({
//     robfpca:::cov_gk(x, type = type, MM = T, smooth = F, cor = T, psd = F)
// })
// system.time({
//     robfpca:::cor_gk_cpp(x, type = type, MM = T)
// })
//
// # system.time({
// #     RobStatTM::locScaleM(x[!is.na(x[, 1]), 1], psi = "huber")
// # })
// # system.time({
// #     loc_scale_M(x[!is.na(x[, 1]), 1], psi = "huber")
// # })
// # mscale_cpp(x[!is.na(x[, 1]), 1], 0.5, 1.56, "bisquare")
// */
