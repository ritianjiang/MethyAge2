// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// dftest
RcppExport SEXP dftest(SEXP df, SEXP vname);
RcppExport SEXP _MethyAge2_dftest(SEXP dfSEXP, SEXP vnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type df(dfSEXP);
    Rcpp::traits::input_parameter< SEXP >::type vname(vnameSEXP);
    rcpp_result_gen = Rcpp::wrap(dftest(df, vname));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP _MethyAge2_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MethyAge2_dftest", (DL_FUNC) &_MethyAge2_dftest, 2},
    {"_MethyAge2_rcpp_hello", (DL_FUNC) &_MethyAge2_rcpp_hello, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_MethyAge2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
