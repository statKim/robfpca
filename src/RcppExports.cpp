// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// get_positive_elements
Rcpp::List get_positive_elements(Eigen::VectorXd Y, Eigen::MatrixXd X, Eigen::VectorXd W);
RcppExport SEXP _robfpca_get_positive_elements(SEXP YSEXP, SEXP XSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(get_positive_elements(Y, X, W));
    return rcpp_result_gen;
END_RCPP
}
// IRLScpp
Rcpp::List IRLScpp(const Eigen::VectorXd Y, const Eigen::MatrixXd X, Rcpp::Nullable<Rcpp::NumericVector> weight_, const int maxit, const double tol, const double k);
RcppExport SEXP _robfpca_IRLScpp(SEXP YSEXP, SEXP XSEXP, SEXP weight_SEXP, SEXP maxitSEXP, SEXP tolSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type weight_(weight_SEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(IRLScpp(Y, X, weight_, maxit, tol, k));
    return rcpp_result_gen;
END_RCPP
}
// locpolysmooth
Eigen::VectorXd locpolysmooth(Eigen::VectorXd Lt, Eigen::VectorXd Ly, Eigen::VectorXd newt, std::string kernel, const double bw, const double k, const int deg);
RcppExport SEXP _robfpca_locpolysmooth(SEXP LtSEXP, SEXP LySEXP, SEXP newtSEXP, SEXP kernelSEXP, SEXP bwSEXP, SEXP kSEXP, SEXP degSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type Lt(LtSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type Ly(LySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type newt(newtSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< const double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< const double >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type deg(degSEXP);
    rcpp_result_gen = Rcpp::wrap(locpolysmooth(Lt, Ly, newt, kernel, bw, k, deg));
    return rcpp_result_gen;
END_RCPP
}
// order_
IntegerVector order_(NumericVector x);
RcppExport SEXP _robfpca_order_(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(order_(x));
    return rcpp_result_gen;
END_RCPP
}
// weighted_quantile_cpp
double weighted_quantile_cpp(NumericVector dat, NumericVector weights, double p);
RcppExport SEXP _robfpca_weighted_quantile_cpp(SEXP datSEXP, SEXP weightsSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type dat(datSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(weighted_quantile_cpp(dat, weights, p));
    return rcpp_result_gen;
END_RCPP
}
// weighted_quantile_top_cpp
double weighted_quantile_top_cpp(NumericVector dat, NumericVector weights, double p);
RcppExport SEXP _robfpca_weighted_quantile_top_cpp(SEXP datSEXP, SEXP weightsSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type dat(datSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(weighted_quantile_top_cpp(dat, weights, p));
    return rcpp_result_gen;
END_RCPP
}
// wrm_fit
NumericVector wrm_fit(NumericVector xdat, NumericVector ydat, double x0, NumericVector weight_vec);
RcppExport SEXP _robfpca_wrm_fit(SEXP xdatSEXP, SEXP ydatSEXP, SEXP x0SEXP, SEXP weight_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xdat(xdatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ydat(ydatSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight_vec(weight_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(wrm_fit(xdat, ydat, x0, weight_vec));
    return rcpp_result_gen;
END_RCPP
}
// get_kernel
NumericVector get_kernel(NumericVector x_i, double x, double h, String method);
RcppExport SEXP _robfpca_get_kernel(SEXP x_iSEXP, SEXP xSEXP, SEXP hSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x_i(x_iSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(get_kernel(x_i, x, h, method));
    return rcpp_result_gen;
END_RCPP
}
// wrm_smooth_cpp
List wrm_smooth_cpp(NumericVector x, NumericVector y, double h, NumericVector xgrid, String kernel);
RcppExport SEXP _robfpca_wrm_smooth_cpp(SEXP xSEXP, SEXP ySEXP, SEXP hSEXP, SEXP xgridSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< String >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(wrm_smooth_cpp(x, y, h, xgrid, kernel));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_robfpca_get_positive_elements", (DL_FUNC) &_robfpca_get_positive_elements, 3},
    {"_robfpca_IRLScpp", (DL_FUNC) &_robfpca_IRLScpp, 6},
    {"_robfpca_locpolysmooth", (DL_FUNC) &_robfpca_locpolysmooth, 7},
    {"_robfpca_order_", (DL_FUNC) &_robfpca_order_, 1},
    {"_robfpca_weighted_quantile_cpp", (DL_FUNC) &_robfpca_weighted_quantile_cpp, 3},
    {"_robfpca_weighted_quantile_top_cpp", (DL_FUNC) &_robfpca_weighted_quantile_top_cpp, 3},
    {"_robfpca_wrm_fit", (DL_FUNC) &_robfpca_wrm_fit, 4},
    {"_robfpca_get_kernel", (DL_FUNC) &_robfpca_get_kernel, 4},
    {"_robfpca_wrm_smooth_cpp", (DL_FUNC) &_robfpca_wrm_smooth_cpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_robfpca(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
