#include <Rcpp.h>
#include "misc.h"
#include "samplers.h"

using namespace Rcpp;

double generate_mean(int N, int N1, int i, const IntegerVector& sim_vector, const NumericMatrix& D, const IntegerVector& k_indices, bool bypass) {

  IntegerVector sampling_vec(sim_vector.size() + N1);
  const IntegerVector::iterator draws_start = sampling_vec.begin() + sim_vector.size();
  IntegerVector::iterator current_draw = sampling_vec.begin();
  IntegerVector::const_iterator sim_itr = sim_vector.begin();

  // Initialize sampling vec with current data
  for(; current_draw != draws_start; ++current_draw, ++sim_itr){
    *current_draw = *sim_itr;
  }

  // Now draw the samples
  for(; current_draw != sampling_vec.end(); ++current_draw){
    // Rcout << sampling_vec << "\n";
    if(bypass) {
      *current_draw = sample_Fi(N, i, k_indices);
    }
    else {
      *current_draw = -1;
      *current_draw = sample_Fik(N, i, sampling_vec.begin(), current_draw, k_indices);
    }
  }
  // Rcout << sampling_vec << "\n";

  IntegerVector indices(2);
  current_draw = draws_start;
  for(int j = 0; j < N1; j++, ++current_draw){
    indices(0) = *current_draw % N;
    indices(1) = *current_draw / N;
    *current_draw = k_indices[indices(1)] * N + indices(0);
  }
  // Rcout << sampling_vec << "\n";

  NumericVector x(N1);
  NumericVector::iterator x_itr = x.begin();
  for(current_draw = draws_start; current_draw !=  sampling_vec.end(); ++current_draw, ++x_itr){
    *x_itr = D[*current_draw];
  }
  // Rcout << x << "\n";
  return mean(x);

  // Rcout << sampling_vec << "\n";
}

NumericVector sample_means(const CharacterVector& P_ks, std::string k, int i, const CharacterVector& assignment, const IntegerVector& k_indices, const NumericMatrix& D) {

  NumericVector output(P_ks.size());
  bool bypass;

  int N = assignment.size();
  int N1 = k_indices.size();

  CharacterVector::const_iterator Pk_itr = P_ks.begin();
  IntegerVector sim_indices;
  // IntegerVector::iterator pep_itr = k_indices.begin();

  NumericVector::iterator means_itr = output.begin();

  for(; Pk_itr != P_ks.end(); ++Pk_itr, ++means_itr) {
    sim_indices = find_string_locs(assignment, as<std::string>(*Pk_itr));
    bypass = as<std::string>(*Pk_itr) == k;
    *means_itr = generate_mean(N, N1, i, sim_indices, D, k_indices, bypass);
  }

  return output;
}

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
NumericVector bayesian_prob(const CharacterVector& assignment, const NumericMatrix& D, std::string k, int B) {

  IntegerVector k_indices = find_string_locs(assignment, k);
  int N1 = k_indices.size();
  NumericVector muk(N1);

  // count number of entities per group
  IntegerVector groups_table = table(assignment);
  CharacterVector all_groups = groups_table.names();

  int index = groups_table.findName(k);
  groups_table(index) = 1;
  // all_groups(index) = "-1";
  groups_table.names() = all_groups;

  NumericVector groups_num = as<NumericVector>(groups_table);

  IntegerVector::iterator mu_itr = k_indices.begin();
  NumericVector mean_vec(B);
  NumericVector probs(N1);

  for(int i = 0; i < N1; ++i) {
    // calculate observed mean
    muk(i) = 0;
    mu_itr = k_indices.begin();
    for(; mu_itr != k_indices.end(); ++mu_itr){
      muk(i) += D[assignment.size()* k_indices(i) + *mu_itr];
    }
    muk(i) = muk(i)/ (N1 - 1);

    // draw k's to sample from:
    CharacterVector k_rep = sample(all_groups, B, true, groups_num);
    mean_vec = sample_means(k_rep, k, i, assignment, k_indices, D);

    for(NumericVector::iterator bi = mean_vec.begin(); bi < mean_vec.end(); ++bi){
      probs(i) += *bi < muk(i);
    }
    probs(i) = probs(i) / B;
  }


  return probs; // DataFrame::create(_["k_ind"] = k_indices, _["i"] = seq_len(N1), _["mu"] = muk, _["probs"] = probs);
}
