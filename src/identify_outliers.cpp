#include <progress.hpp>
#include <progress_bar.hpp>
#include <Rcpp.h>
#include "bayesian_probability.h"
#include "empirical_prob.h"
#include "misc.h"
// [[Rcpp::depends(RcppProgress)]]

using namespace Rcpp;

//' @title Identify Outliers
//' 
//' @description Calculate the probability that each entity in a group with at least 2 entities actually belongs
//' @param assignment A character vector containing the group to which each entity has been classified
//' @param D A numeric vector containing the distance or quasi-distance between the entities.
//' @param B The number of simulated draws to use in determining the distribution of the means.
//' @param min_group_size The smallest number of entities required to perform filtering.  (Must be at least 2)
//' @param labels entity-specific labels.
//'
//' @details
//' For each entity in group k, this function estimates the probability the mean distance between it and another entity in group k has the same distribution as a mean distance from it and another group.
//' Under the null hypothesis, the mean distance is no different, and we remove it from the group.
//' Under the alterantive hypothesis, the mean distance is closer, and we conclude that the entity truly belongs to the group.
//'
//' @export
// [[Rcpp::export]]
DataFrame identify_outliers(const CharacterVector& assignment, const NumericMatrix& D, int B, int min_group_size, std::string prior, const CharacterVector& labels, bool display_progress = false, double omega = 1, double omega0 = 1) {
  // int N = assignment.size();

  IntegerVector protein_table = table(assignment);
  protein_table = protein_table[protein_table >= min_group_size];
  CharacterVector proteins = protein_table.names();
  int Nr = sum(protein_table);
  int K = protein_table.size();

  CharacterVector::iterator prot = proteins.begin();

  // gather up the results into these vectors:
  CharacterVector protein_id(Nr); CharacterVector::iterator protein_i = protein_id.begin();
  CharacterVector label(Nr); CharacterVector::iterator label_i = label.begin();
  NumericVector prob(Nr); NumericVector::iterator prob_i = prob.begin();
  NumericVector prob2(Nr); NumericVector::iterator prob2_i = prob2.begin();
  IntegerVector index(Nr); IntegerVector::iterator index_i = index.begin();

  // Defined to help iteration
  IntegerVector prot_breaks = cumsum(protein_table); IntegerVector::iterator break_i = prot_breaks.begin();
  int i = 0;
  Progress p(K, display_progress);

  for(; prot != proteins.end(); ++prot, ++break_i) {
    // progress bar
    if(Progress::check_abort()){
      stop("Process interupted by user.\n");
    }
    p.increment();
    
    IntegerVector peptides = find_string_locs(assignment, as<std::string>(*prot)); IntegerVector::iterator peptide_i = peptides.begin();
    NumericVector pep_probs = bayesian_prob(omega, omega0, assignment, D, as<std::string>(*prot), B, prior);
    NumericVector pep_probs2 = empirical_prob(assignment, D, as<std::string>(*prot));
    NumericVector::iterator vec_prob_i = pep_probs.begin();
    NumericVector::iterator vec_prob2_i = pep_probs2.begin();

    for(; i < *break_i; ++i, ++protein_i, ++label_i, ++prob_i, ++prob2_i, ++index_i, ++vec_prob2_i, ++vec_prob_i, ++peptide_i) {
      *protein_i = *prot;
      *label_i = labels(*peptide_i);
      *index_i = *peptide_i;
      *prob_i = *vec_prob_i;
      *prob2_i = *vec_prob2_i;
    }
  }

  return DataFrame::create(
    _["group"] = protein_id,
    _["entity"] = label,
    _["index"] = index,
    _["bayesian_prob"] = prob,
    _["empirical_prob"] = prob2);
}
