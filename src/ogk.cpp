#include <Rcpp.h>

using namespace Rcpp;

// MASS::fitdistr()
// [[Rcpp::export]]
Rcpp::List fitdistr_cpp(Rcpp::NumericVector x,
                        Rcpp::String densfun,
                        int df) {
    // Obtain environment containing function
    Rcpp::Environment env("package:MASS");

    // Make function callable from C++
    Rcpp::Function fitdistr_r = env["fitdistr"];

    // Call the function and receive its list output
    Rcpp::List res = fitdistr_r(Rcpp::_["x"] = x,
                                Rcpp::_["densfun"] = densfun,
                                Rcpp::_["df"] = df);

    // Return test object in list structure
    return res;
}


// RobStatTM::locScaleM()
// [[Rcpp::export]]
Rcpp::List locScaleM_cpp(Rcpp::NumericVector x,
                         Rcpp::String psi,
                         bool na_rm) {
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
Rcpp::List quantile_cpp(Rcpp::NumericVector x,
                    double probs) {
    // Obtain environment containing function
    Rcpp::Environment env("package:stats");

    // Make function callable from C++
    Rcpp::Function ftn_r = env["quantile"];

    // Call the function and receive its list output
    Rcpp::List res = ftn_r(Rcpp::_["x"] = x,
                           Rcpp::_["probs"] = probs);

    // Return test object in list structure
    return res;
}



// // [[Rcpp::export]]
// double get_sigma2_rob_cpp(Rcpp::NumericVector v) {
//     // 2*Hampel loss
//     // cut-off values are obtained from Sinova(2018)
//     v = abs(v);
//     Rcpp::NumericVector v2 = v[!is_na(v)];
//     double a = median(v2);
//     b <- quantile(v2, 0.75)
//     c <- quantile(v2, 0.85)
//     x <- ifelse(v < a, v^2,
//                 ifelse(v < b, 2*a*v - a^2,
//                        ifelse(v < c, a*(v-c)^2/(b-c) + a*(b+c-a), a*(b+c-a))))
//     double res = mean(x, na.rm = T)
//
//     return res;
// }






//
// /*** R
// x2 <- rt(250, df = 9)
//
// fitdistr(x2, "t", df = 9)
// fitdistr_cpp(x2, "t", df = 9)
//
// x2[c(1, 100, 30, 50)] <- NA
// locScaleM(x2, psi = "huber", na.rm = T)
// locScaleM_cpp(x2, psi = "huber", na_rm = T)
//
// quantile_cpp(1:10, 0.5)
// */
