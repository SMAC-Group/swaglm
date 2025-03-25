// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fast_toeplitz_matrix_from_vector_cpp
arma::mat fast_toeplitz_matrix_from_vector_cpp(const arma::vec& v);
RcppExport SEXP _swaglm_fast_toeplitz_matrix_from_vector_cpp(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_toeplitz_matrix_from_vector_cpp(v));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_swaglm_fast_toeplitz_matrix_from_vector_cpp", (DL_FUNC) &_swaglm_fast_toeplitz_matrix_from_vector_cpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_swaglm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
