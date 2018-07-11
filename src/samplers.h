#ifndef SAMPLERS_H
#define SAMPLERS_H

#include <Rcpp.h>

int sample_F0(int N, int i, const Rcpp::IntegerVector::iterator& row_indices_start, const Rcpp::IntegerVector::iterator& row_indices_end);
int sample_Fi(int N, int i, const Rcpp::IntegerVector& k_indices);
int sample_Fik(int N, int i, const Rcpp::IntegerVector::iterator& sim_vector_start, const Rcpp::IntegerVector::iterator& sim_vector_end, const Rcpp::IntegerVector& k_indices);


# endif