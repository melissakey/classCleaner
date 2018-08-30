#include <Rcpp.h>
#include "misc.h"
#include "get_observed_means.h"
#include "get_estimated_prob.h"

using namespace Rcpp;

//' Bayesian Filtering Probablity
//' Calculate the probability of P(i in Group k | data) for all entities in Group k using simulated draws from a hierarchical Dirichlet sampling scheme.
//' @param assignment A character vector containing the group to which each entity has been classified
//' @param D A numeric vector containing the distance or quasi-distance between the entities.
//' @param k A string identifying the group on which to apply the filtering.
//' @param B The number of simulated draws to use in determining the distribution of the means.
//'
//' @details
//' For each entity in group k, this function estimates the probability the mean distance between it and another entity in group k has the same distribution as a mean distance from it and another group.
//' Under the null hypothesis, the mean distance is no different, and we remove it from the group.
//' Under the alterantive hypothesis, the mean distance is closer, and we conclude that the entity truly belongs to the group.
//'
//' @export
//' @examples
//' # we use a matrix of indices to illustrate how the draws are made:
//' set.seed(23)
//' N <- 20
//' K <- 4
//' assignment <- sample(LETTERS[1:K], N, replace = TRUE)
//'
//' m <- matrix(runif(N^2), N, N)
//' diag(m) = 0
//' m <- m %*% m
//' k <- 'A'
//' i <- 1
//' probs <- estimate_probability(assignment, m, "A", 1000)
//' hist(means)


// [[Rcpp::export]]
NumericVector bayesian_prob(double omega, double omega0, const CharacterVector& assignment, const NumericMatrix& D, std::string k, int B, std::string prior) {


  // **********************************************
  //
  // Step 1: calculate empirical estiamte of mean
  // for each entity in the selected class
  //
  // **********************************************
  
  // get vector sizes
  int N = assignment.size();
  
  // identify indices within selected class
  IntegerVector k_indices = find_string_locs(assignment, k);
  // int N1 = k_indices.size();
  
  NumericVector mu_k = get_observed_means(N, omega, omega0, k_indices, D);

  IntegerVector class_table = table(assignment);
  NumericVector probs = get_estimated_prob(omega, omega0, k_indices, class_table, k, assignment, D, B, mu_k, prior);
  
  return probs; // DataFrame::create(_["k_ind"] = k_indices, _["i"] = seq_len(N1), _["mu"] = muk, _["probs"] = probs);
}
