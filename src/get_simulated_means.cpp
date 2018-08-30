#include <Rcpp.h>
#include "misc.h"
#include "samplers.h"

using namespace Rcpp;

double generate_mean(double omega, double omega0, int N, int K, int i, const IntegerVector& s_indices, const NumericMatrix& D, const IntegerVector k_indices, std::string prior) {

  int Nk = k_indices.size();
  
  IntegerVector indices = sample_Fik(N, K, i, omega, omega0, k_indices, s_indices, prior);

  // vector for distances
  NumericVector d(Nk);
  for(int j = 0; j < Nk; j++) d(j) = D[indices(j)];
  
  return mean(d);

}

// [[Rcpp::export]]
NumericVector get_simulated_means(double omega, double omega0, int K, const Rcpp::CharacterVector& simulated_ks,  std::string k, int i, const Rcpp::CharacterVector& assignment, const Rcpp::IntegerVector& k_indices, const Rcpp::NumericMatrix& D, std::string prior){

  // shortcuts
  int N = assignment.size();
  
  // setup output vector of simulated means
  NumericVector simulated_means(simulated_ks.size());
  NumericVector::iterator means_itr = simulated_means.begin();
 
  CharacterVector::const_iterator sim_k_itr = simulated_ks.begin();

  for(; sim_k_itr != simulated_ks.end(); ++sim_k_itr, ++means_itr) {

    IntegerVector sim_indices = find_string_locs(assignment, as<std::string>(*sim_k_itr));
    if(sim_indices.size() == 0) return 0;
    
    *means_itr = generate_mean(omega, omega0, N, K, i, sim_indices, D, k_indices, prior);
  }
  
  return simulated_means;
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

get_simulated_means(1, 1, 5, c('3', '2'), k, 0,  assignment, 0:4,  D, "Fk")


*/