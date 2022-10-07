// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_clip_dcd_optimizer
Rcpp::List cpp_clip_dcd_optimizer(arma::mat H, arma::mat q, arma::mat lb, arma::mat ub, double eps, unsigned int max_steps, arma::mat u);
RcppExport SEXP _manysvms_cpp_clip_dcd_optimizer(SEXP HSEXP, SEXP qSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP epsSEXP, SEXP max_stepsSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type q(qSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_steps(max_stepsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_clip_dcd_optimizer(H, q, lb, ub, eps, max_steps, u));
    return rcpp_result_gen;
END_RCPP
}
// cpp_rbf_kernel
SEXP cpp_rbf_kernel(arma::mat x1, arma::mat x2, double gamma);
RcppExport SEXP _manysvms_cpp_rbf_kernel(SEXP x1SEXP, SEXP x2SEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_rbf_kernel(x1, x2, gamma));
    return rcpp_result_gen;
END_RCPP
}
// cpp_poly_kernel
SEXP cpp_poly_kernel(arma::mat x1, arma::mat x2, double gamma, int degree, double coef0);
RcppExport SEXP _manysvms_cpp_poly_kernel(SEXP x1SEXP, SEXP x2SEXP, SEXP gammaSEXP, SEXP degreeSEXP, SEXP coef0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< double >::type coef0(coef0SEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_poly_kernel(x1, x2, gamma, degree, coef0));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_manysvms_cpp_clip_dcd_optimizer", (DL_FUNC) &_manysvms_cpp_clip_dcd_optimizer, 7},
    {"_manysvms_cpp_rbf_kernel", (DL_FUNC) &_manysvms_cpp_rbf_kernel, 3},
    {"_manysvms_cpp_poly_kernel", (DL_FUNC) &_manysvms_cpp_poly_kernel, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_manysvms(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
