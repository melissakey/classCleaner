#ifndef GET_ESTIMATED_PROB_H
#define GET_ESTIMATED_PROB_H

#include <Rcpp.h>

Rcpp::NumericVector get_estimated_prob(double omega, double omega0, const Rcpp::IntegerVector& k_indices, const Rcpp::IntegerVector& class_table, std::string k, const Rcpp::CharacterVector& assignment,  const Rcpp::NumericMatrix& D, int B, const Rcpp::NumericVector& observed_means, std::string prior);

# endif