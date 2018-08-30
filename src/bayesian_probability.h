#ifndef BAYESIAN_PROBABILITY_H
#define BAYESIAN_PROBABILITY_H

#include <Rcpp.h>

Rcpp::NumericVector bayesian_prob(double omega, double omega0, const Rcpp::CharacterVector& assignment, const Rcpp::NumericMatrix& D, std::string k, int B, std::string prior);
#endif
