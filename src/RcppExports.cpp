// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// psi
NumericVector psi(const NumericVector& x, const NumericVector& y);
RcppExport SEXP _classCleaner_psi(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(psi(x, y));
    return rcpp_result_gen;
END_RCPP
}
// sim_by_class
arma::mat sim_by_class(arma::uword n, const arma::uvec& Nk, const arma::colvec& s, double tau, const arma::mat& rho);
RcppExport SEXP _classCleaner_sim_by_class(SEXP nSEXP, SEXP NkSEXP, SEXP sSEXP, SEXP tauSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type Nk(NkSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_by_class(n, Nk, s, tau, rho));
    return rcpp_result_gen;
END_RCPP
}
// sim_by_instance
arma::mat sim_by_instance(arma::uword n, const arma::uvec& Nk, const arma::vec& s, const arma::mat& rho);
RcppExport SEXP _classCleaner_sim_by_instance(SEXP nSEXP, SEXP NkSEXP, SEXP sSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type Nk(NkSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_by_instance(n, Nk, s, rho));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_classCleaner_psi", (DL_FUNC) &_classCleaner_psi, 2},
    {"_classCleaner_sim_by_class", (DL_FUNC) &_classCleaner_sim_by_class, 5},
    {"_classCleaner_sim_by_instance", (DL_FUNC) &_classCleaner_sim_by_instance, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_classCleaner(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
