#ifndef EMPIRICALPROBABILITY_H
#define EMPIRICALPROBABILITY_H

#include <Rcpp.h>
#include "misc.h"

Rcpp::NumericVector empirical_prob(const Rcpp::CharacterVector& assignment, const Rcpp::NumericMatrix& D, std::string k);
#endif