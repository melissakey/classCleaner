#include <Rcpp.h>
using namespace Rcpp;

// Calculate prior mean
// [[Rcpp::export]]
NumericVector get_observed_means(int N, double omega, double omega0, const IntegerVector& k_indices, const NumericMatrix& D) {

  int N1 = k_indices.size();
  
  
  // Identify indices outside selected class
  IntegerVector all_indices = seq(0, N - 1);
  IntegerVector out_indices = setdiff(all_indices, k_indices);
  
  // Define iterators for each indice
  IntegerVector::const_iterator it = k_indices.begin();
  IntegerVector::iterator jt = out_indices.begin();
  
 
  // ************************************************  
  // 
  //  Step 1 - calculate prior mean for each i
  // 
  // ************************************************
  
  // start by calculating the mean of all non-within-class distances for each columns
  NumericVector observed_means(N1);
  NumericVector::iterator dt = observed_means.begin();
  
  for(; it != k_indices.end(); it++, dt++) {
    for(jt = out_indices.begin(); jt != out_indices.end(); jt++) {
      *dt += D(*jt, *it);
      // Rcout << observed_means << "\n";
    }
    *dt = *dt / out_indices.size();
  }
  // get the prior of this by calculating the overall mean
  double sum_means = sum(observed_means);
  // the prior mean is a weighted sum of the overall mean and the 
  NumericVector prior_mean = (observed_means * (N1 - 1) + omega0 * (sum_means - observed_means)/ (N1 - 1)) / (N1 - 1 + omega0);
  // Rcout << "sum_means = " << sum_means << "\n";
  // Rcout << "prior_mean = " << prior_mean << "\n";
  // // ************************************************  
  // 
  //  Step 2 - calculate the observed within-class mean distance
  // 
  // ************************************************  
  
  NumericVector muk(N1);
  for(int i = 0; i < N1; ++i) {
    // calculate observed mean
    muk(i) = 0;
    it = k_indices.begin();
    for(; it != k_indices.end(); ++it){
      muk(i) += D(*it, k_indices(i));
    }
    muk(i) = (muk(i) + omega * prior_mean(i)) / (N1 - 1 + omega);
  }
  
  return muk;
}
/***R
m <- 1 - diag(6)
m[lower.tri(m)] <- runif(15)
m <- m + t(m)

N <- nrow(m)
omega = 1
omega0 = 1
k_indices = 1:3

get_observed_means(N, omega, omega0, k_indices, m)

*/