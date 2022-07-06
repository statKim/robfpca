// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fitdistr_cpp
Rcpp::List fitdistr_cpp(Rcpp::NumericVector x, Rcpp::String densfun, int df);
RcppExport SEXP _robfpca_fitdistr_cpp(SEXP xSEXP, SEXP densfunSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type densfun(densfunSEXP);
    Rcpp::traits::input_parameter< int >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(fitdistr_cpp(x, densfun, df));
    return rcpp_result_gen;
END_RCPP
}
// locScaleM_cpp
Rcpp::List locScaleM_cpp(Rcpp::NumericVector x, Rcpp::String psi, bool na_rm);
RcppExport SEXP _robfpca_locScaleM_cpp(SEXP xSEXP, SEXP psiSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(locScaleM_cpp(x, psi, na_rm));
    return rcpp_result_gen;
END_RCPP
}
// quantile_cpp
Rcpp::List quantile_cpp(Rcpp::NumericVector x, double probs);
RcppExport SEXP _robfpca_quantile_cpp(SEXP xSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(quantile_cpp(x, probs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_robfpca_fitdistr_cpp", (DL_FUNC) &_robfpca_fitdistr_cpp, 3},
    {"_robfpca_locScaleM_cpp", (DL_FUNC) &_robfpca_locScaleM_cpp, 3},
    {"_robfpca_quantile_cpp", (DL_FUNC) &_robfpca_quantile_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_robfpca(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
