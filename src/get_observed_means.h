#ifndef GET_OBSERVED_MEANS_H
#define GET_OBSERVED_MEANS_H

#include <Rcpp.h>

Rcpp::NumericVector get_observed_means(int N, double omega, double omega0, const Rcpp::IntegerVector& k_indices, const Rcpp::NumericMatrix& D);
  
# endif