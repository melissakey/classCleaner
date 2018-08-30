#ifndef GET_SIMULATED_MEANS_H
#define GET_SIMULATED_MEANS_H

#include <Rcpp.h>

double generate_mean(double omega, double omega0, int N, int K, int i, const Rcpp::IntegerVector& s_indices, const Rcpp::NumericMatrix& D, const Rcpp::IntegerVector k_indices, std::string prior);
Rcpp::NumericVector get_simulated_means(double omega, double omega0, int K, const Rcpp::CharacterVector& simulated_ks,  std::string k, int i, const Rcpp::CharacterVector& assignment, const Rcpp::IntegerVector& k_indices, const Rcpp::NumericMatrix& D, std::string prior);

# endif