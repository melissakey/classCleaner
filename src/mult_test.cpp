#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


//' @title test function
//' 
//' @description this is a trivial example.
//' 
//' @param vec1 a numeric vector
//' @details a trivial example to figure out why matrix multiplication is giving errors.
//' 
//' @export
// [[Rcpp::export]]
arma::mat test(const arma::vec& vec1){
  return vec1 * vec1.t();
}

/***R

test(1:5)

*/