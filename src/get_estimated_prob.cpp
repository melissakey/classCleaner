#include <Rcpp.h>
#include "get_simulated_means.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_estimated_prob(double omega, double omega0, const IntegerVector& k_indices, const IntegerVector& class_table, std::string k, const CharacterVector& assignment,  const NumericMatrix& D, int B, const NumericVector& observed_means, std::string prior) {
  
  // bookkeeping
  int Nk = k_indices.size();
  int K = class_table.size();
  // int N = assignment.size();
  
  // vector for output
  NumericVector prob(Nk);
  
  
  // // temporary variables for deciding whether or not to sample from prior
  // LogicalVector ft = LogicalVector::create(false, true);
  // LogicalVector draw_from_prior(B);
  // LogicalVector::iterator prior_itr = draw_from_prior.begin();
  // NumericVector prior_prob = NumericVector::create(K, omega) / (K + omega);
  
  // create probability table from class assignment table and make sure we don't draw the current class
  NumericVector class_prob = as<NumericVector>(class_table);
  class_prob(class_table.findName(k)) = 0;
  CharacterVector class_names = class_table.names();
  
  // temporary variables for holding means
  NumericVector simulated_means(B);
  NumericVector::iterator sm_itr = simulated_means.begin();
  CharacterVector simulated_k(B);
  
  for(int i = 0; i < Nk; i++) {
    sm_itr = simulated_means.begin();
    // determine whether each observation comes from Fik a prior distribution
    // draw_from_prior = sample(ft, B, true, prior_prob);
    simulated_k = sample(class_names, B, true, class_prob);
    
    simulated_means = get_simulated_means(omega, omega0, K, simulated_k, k, i, assignment, k_indices, D, prior);
    
    prob(i) = mean(simulated_means < observed_means(i));
  }
  
  return prob;
}


/***R

assignment <- rep(1:4, each = 5)

D <- 1- diag(20)
D[lower.tri(D)] <- runif(190)
D <- D + t(D)
k <- '1'
omega <- 1
mu.junk <- colSums(D) / 19
B <- 5

tmp <- get_estimated_prob(1, 1, 0:4, table(assignment), k, assignment, D, B, mu.junk[1:4], "Fk")

*/