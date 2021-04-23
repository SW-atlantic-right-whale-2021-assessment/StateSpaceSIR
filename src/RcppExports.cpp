// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// GENERALIZED_LOGISTIC
List GENERALIZED_LOGISTIC(double r_max, double K, double N1, double z, double start_yr, double num_Yrs, NumericVector catches, NumericVector proc_error, double MVP);
RcppExport SEXP _StateSpaceSIR_GENERALIZED_LOGISTIC(SEXP r_maxSEXP, SEXP KSEXP, SEXP N1SEXP, SEXP zSEXP, SEXP start_yrSEXP, SEXP num_YrsSEXP, SEXP catchesSEXP, SEXP proc_errorSEXP, SEXP MVPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type r_max(r_maxSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type N1(N1SEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type start_yr(start_yrSEXP);
    Rcpp::traits::input_parameter< double >::type num_Yrs(num_YrsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type catches(catchesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type proc_error(proc_errorSEXP);
    Rcpp::traits::input_parameter< double >::type MVP(MVPSEXP);
    rcpp_result_gen = Rcpp::wrap(GENERALIZED_LOGISTIC(r_max, K, N1, z, start_yr, num_Yrs, catches, proc_error, MVP));
    return rcpp_result_gen;
END_RCPP
}
// dlnorm_zerb
NumericVector dlnorm_zerb(NumericVector x, NumericVector meanlog, NumericVector sdlog, bool return_log);
RcppExport SEXP _StateSpaceSIR_dlnorm_zerb(SEXP xSEXP, SEXP meanlogSEXP, SEXP sdlogSEXP, SEXP return_logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type meanlog(meanlogSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sdlog(sdlogSEXP);
    Rcpp::traits::input_parameter< bool >::type return_log(return_logSEXP);
    rcpp_result_gen = Rcpp::wrap(dlnorm_zerb(x, meanlog, sdlog, return_log));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StateSpaceSIR_GENERALIZED_LOGISTIC", (DL_FUNC) &_StateSpaceSIR_GENERALIZED_LOGISTIC, 9},
    {"_StateSpaceSIR_dlnorm_zerb", (DL_FUNC) &_StateSpaceSIR_dlnorm_zerb, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_StateSpaceSIR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
