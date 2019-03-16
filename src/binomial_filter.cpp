// #include <Rcpp.h>
// #include "misc.h"
// using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// // [[Rcpp::export]]
// DataFrame binomial_filter(const CharacterVector& assignment, const NumericMatrix& D, std::string k, const NumericVector& cutoff) {
//   IntegerVector k_ind = find_string_locs(assignment, k);
//   IntegerVector::iterator k_p = k_ind.begin();
//   NumericVector::const_iterator cutoff_p = cutoff.begin();
//   int N_k = k_ind.size();
//   int N_o = assignment.size() - N_k;
//   
//   // NumericVector D1(N_k * (N_k - 1) / 2);
//   // NumericVector D2(N_k * N_o);
//   IntegerMatrix Zi(N_k, cutoff.size());
// 
//   
//   for(int j = 0; k_p < k_ind.end(); k_p++, j++) {
//     for(int i = 0; i < *k_ind.end(); i++) {
//       Zi(_,j) = Z(_,j) + ()
//     }
//   }
//   
// }
