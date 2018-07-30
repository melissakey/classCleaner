#ifndef BAYESIANPROBABILITY_H
#define BAYESIANPROBABILITY_H

#include <Rcpp.h>
#include "misc.h"
#include "samplers.h"

double generate_mean(int N, int N1, int i, const Rcpp::IntegerVector& sim_vector, const Rcpp::NumericMatrix& D, const Rcpp::IntegerVector& k_indices, bool bypass, std::string prior);
Rcpp::NumericVector sample_means(const Rcpp::CharacterVector& P_ks, std::string k, int i, const Rcpp::CharacterVector& assignment, const Rcpp::IntegerVector& k_indices, const Rcpp::NumericMatrix& D, std::string prior);
Rcpp::NumericVector bayesian_prob(const Rcpp::CharacterVector& assignment, const Rcpp::NumericMatrix& D, std::string k, int B, std::string prior);

#endif
