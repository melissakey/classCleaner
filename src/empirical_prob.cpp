#include <Rcpp.h>
#include "misc.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector empirical_prob(const CharacterVector& assignment, const NumericMatrix& D, std::string k) {

  IntegerVector k_indices = find_string_locs(assignment, k);
  int N1 = k_indices.size();
  int N = assignment.size();
  // NumericVector muk(N1);

  // count number of entities per group
  IntegerVector groups_table = table(assignment);
  // CharacterVector all_groups = groups_table.names();

  int index = groups_table.findName(k);
  int K = groups_table.size();
  // groups_table(index) = 1;
  // all_groups(index) = "-1";
  // groups_table.names() = all_groups;

  NumericMatrix group_sums(N1, K);
  // NumericMatrix::ConstColumn D_column = D(_, k_indices(0));
  // CharacterVector::const_iterator a_itr = assignment.begin();
  IntegerVector::iterator k_itr = k_indices.begin();
  int kg_index;


  // calculate group-sums
  for(int i = 0; i < N; ++i){
    kg_index = groups_table.findName(as<std::string>(assignment(i)));

    // figure out which column of summation matrix to add stuff to
    NumericMatrix::Column sums_column = group_sums(_, kg_index);
    NumericVector::iterator entry_itr = sums_column.begin();

    for(entry_itr = sums_column.begin(),  k_itr = k_indices.begin();
      entry_itr != sums_column.end();
      ++entry_itr, ++k_itr) {
      *entry_itr += D(i, *k_itr);
    }
  }

  NumericVector muk = group_sums(_, index);
  NumericVector sum_Fi = rowSums(group_sums) - group_sums(_, index);
  double sum_F0 = sum(sum_Fi);
  NumericVector mu_F0 = (sum_F0 - sum_Fi) / ((N - N1) * (N1 - 1));
  NumericVector mu_Fi = (sum_Fi + mu_F0) / (1 + (N - N1));
  NumericMatrix mu_Fik(N1, K - 1);
  LogicalVector tmp(N1);
  IntegerVector probs(N1);
  // Rcout << "sum_Fi = "  << sum_Fi << "\n";
  // Rcout << "mu_F0 = " << mu_F0 << "\n";
  // Rcout << "mu_Fi = " << mu_Fi << "\n";

  muk = (muk + mu_Fi)/N1;


  for(int i = 0; i < K; i++)
  {
    if(i < index){
      mu_Fik(_, i) = (group_sums(_, i) + mu_Fi) / (1 + groups_table(i));
      tmp = mu_Fik(_, i) < muk;
    }
    if(i > index){
      mu_Fik(_, i - 1) = (group_sums(_, i) + mu_Fi) / (1 + groups_table(i));
      tmp = mu_Fik(_, i - 1) < muk;
    }
    probs += as<IntegerVector>(tmp);
  }

  return as<NumericVector>(probs) / K;
}

