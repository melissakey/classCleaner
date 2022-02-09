#include <Rcpp.h>
using namespace Rcpp;

Rcpp::IntegerVector find_string_locs(Rcpp::CharacterVector input_vector, std::string id){
  Rcpp::IntegerVector id_locs; // Setup storage for found IDs

  for(int i =0; i < input_vector.size(); i++) // Loop through input
    if(input_vector[i] == id) // check if input matches target
      id_locs.push_back(i);

    return id_locs; // send locations to R (c++ index shift!)
}

Rcpp::NumericVector quantile(const Rcpp::NumericVector x, const NumericVector q) {
  NumericVector y = clone(x);
  NumericVector ans(q.size());

  NumericVector::iterator a_it = ans.begin();
  NumericVector::const_iterator q_it = q.begin();
  
  std::sort(y.begin(), y.end());
  NumericVector k = floor((x.size() - 1) * q + 1) - 1;
  NumericVector::iterator k_it = k.begin();

  for(; a_it < ans.end(); a_it++, q_it++, k_it++){
    *a_it = y(*k_it) * *q_it + y(*k_it + 1) * (1 - *q_it);
  }
  
  return ans;
}

/*** R
x <- runif(100)

Cquantile(x, 0.05)
quantile(x, 0.05)
*/