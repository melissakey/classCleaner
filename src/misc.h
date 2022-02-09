#ifndef MISC_H
#define MISC_H

#include <Rcpp.h>
Rcpp::IntegerVector find_string_locs(Rcpp::CharacterVector input_vector, std::string id);

#endif