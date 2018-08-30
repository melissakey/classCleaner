#include <Rcpp.h>
#include <string>
using namespace Rcpp;


IntegerMatrix sample_F0(int n, int Nk, const IntegerVector& o_indices) {

  
  IntegerMatrix indices(n, 2);

  IntegerVector col_indices = seq(0, Nk - 1);

  IntegerMatrix::Column row_ind = indices(_, 0);
  IntegerMatrix::Column col_ind = indices(_, 1);
  
  col_ind = sample(col_indices, n, true);
  row_ind = sample(o_indices, n, true);
  
  return indices;
}

IntegerVector sample_Fk(int N, int n, int i, double omega0, const IntegerVector& k_indices, const IntegerVector& s_indices) {

  int Nk = k_indices.size();
  int Ns = s_indices.size();

  // define output vector
  IntegerVector sampled_indices(n);

  // for this sampler, we first identify which rows and columns from (the k-indexed columnns of) D we are interested in
  IntegerVector col_indices = seq_len(Nk - 1) - 1;
  for(int j = i; j < Nk - 1; j++) col_indices(j) += 1;

  IntegerVector row_indices(Ns);//  / N;
  IntegerVector::iterator r = row_indices.begin();
  IntegerVector::const_iterator s = s_indices.begin();
  for(; r != row_indices.end(); ++r, ++s) {
    *r = *s % N;
  }

  // standard procedure at this point:
  // temporary variables for deciding whether or not to sample from prior
  IntegerVector sign = IntegerVector::create(-1, 1);
  NumericVector prior_prob = NumericVector::create(omega0, Nk - 1);

  // decide where each sample is generated from
  IntegerVector keep_index = sample(sign, n, true, prior_prob);

  int n0 = sum(keep_index < 0);

  if(n0 < n) {
    IntegerVector row_ind = sample(row_indices, n - n0, true);
    IntegerVector s_ind = sample(col_indices, n - n0, true);
    
    // Rcout << "row_ind = " << row_ind << "\n";
    // Rcout << "s_ind = " << s_ind << "\n";
    s_ind = s_ind * N + row_ind;
    // Rcout << "s_ind = " << s_ind << "\n";
    
    sampled_indices[keep_index > 0] = s_ind;
  }
  if(n0 > 0) {
    IntegerVector all_indices = seq(0, N - 1);
    IntegerVector o_indices = setdiff(setdiff(all_indices, k_indices), s_indices);
    // Rcout << "o_indices: " << o_indices << "\n";

    IntegerMatrix sample_matrix = sample_F0(n0, Nk, o_indices);
    IntegerVector tmp = sample_matrix(_, 0) + sample_matrix(_, 1) * N;
    sampled_indices[keep_index < 0] = tmp;
  }

  return sampled_indices;
}

IntegerVector sample_Fi(int N, int K, int n, int i, double omega0, const IntegerVector& k_indices) {

  int Nk = k_indices.size();

  IntegerVector all_indices = seq(0, N - 1);
  IntegerVector o_indices = setdiff(all_indices, k_indices);
  IntegerVector sampled_indices(n);

  // temporary variables for deciding whether or not to sample from prior
  IntegerVector sign = IntegerVector::create(-1, 1);
  NumericVector prior_prob = NumericVector::create(omega0, K);

  // decide where each sample is generated from
  IntegerVector keep_index = sample(sign, n, true, prior_prob);


 // IntegerVector::iterator sample_index = sampled_indices.begin();
  int n0 = sum(keep_index < 0);

  if(n0 > 0) {
    IntegerMatrix sample_matrix = sample_F0(n0, Nk - 1, o_indices);
    IntegerVector tmp = sample_matrix(_, 1);
    for(int j = 0; j < n0; j++) if(tmp(j) >= i) tmp(j) += 1; 
    tmp = tmp * N + sample_matrix(_, 0);
    
    sampled_indices[keep_index < 0] = tmp;
  }
  if(n0 < n) sampled_indices[keep_index > 0] = sample(o_indices, n - n0, true);

  return sampled_indices;
}

IntegerVector sample_Fik(int N, int K, int i, double omega, double omega0, const IntegerVector& k_indices, const IntegerVector& s_indices, std::string prior) {

  int Ns = s_indices.size();
  int Nk = k_indices.size();

  IntegerVector sampled_indices(Nk);

  // temporary variables for deciding whether or not to sample from prior
  IntegerVector sign = IntegerVector::create(-1, 1);
  NumericVector prior_prob = NumericVector::create(omega, Ns);

  // decide where each sample is generated from
  IntegerVector keep_index = sample(sign, Nk, true, prior_prob);

  // count how many indices need to be drawn from the prior
  int n = sum(keep_index < 0);

  if(n > 0){
    if(prior == "Fi") {
      sampled_indices[keep_index < 0] = sample_Fi(N, K, n, i, omega0, k_indices);
    }
    else if(prior == "Fk") {
      sampled_indices[keep_index < 0] = sample_Fk(N, n, i, omega0, k_indices, s_indices);
    }
    else stop("invalid prior distribution for Fik");
  }
  if(n < Nk) {
    sampled_indices[keep_index > 0] = sample(s_indices, Nk - n, true);
  }
  return sampled_indices;

}

/*** R
# This is R code
N <- 10
Nk <- 3
i <- 0
k_ind <- c(0, 5, 6)
s_ind <- 1:3

# sample_F0(N, 10, 3, i, c(4, 7:9))
tmp0 <- replicate(10, sample_Fk(N, 10, i, 1, k_ind, s_ind))
# tmp1 <- replicate(10, sample_Fi(N, 3, 2, i, 1, k_ind))

# sample_Fik_wrapper(N, i, sim_vector, k_ind)
tmp <- replicate(10000, sample_Fik(N, 3, i, 1, 1, k_ind, s_ind, "Fk"))
# table(tmp)
# 
tmp2 <- replicate(10000, sample_Fik(N, 3, i, 1, 1, k_ind, s_ind, "Fi"))
# table(tmp2)
# 
hist(tmp2, breaks = (-1:max(tmp2)) + 0.5)
hist(tmp, breaks = (-1:max(tmp)) + 0.5)
# 
# tmp %>% table() %>% prop.table()
*/
