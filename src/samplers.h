#ifndef SAMPLERS_H
#define SAMPLERS_H

#include <Rcpp.h>

Rcpp::IntegerMatrix sample_F0(int n, int Nk, const Rcpp::IntegerVector& o_indices);
Rcpp::IntegerVector sample_Fk(int N, int n, int i, double omega0, const Rcpp::IntegerVector& k_indices, const Rcpp::IntegerVector& s_indices);
Rcpp::IntegerVector sample_Fi(int N, int K, int n, int i, double omega0, const Rcpp::IntegerVector& k_indices);
Rcpp::IntegerVector sample_Fik(int N, int K, int i, double omega, double omega0, const Rcpp::IntegerVector& k_indices, const Rcpp::IntegerVector& s_indices, std::string prior);

# endif