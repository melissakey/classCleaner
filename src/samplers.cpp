#include <Rcpp.h>
using namespace Rcpp;


int sample_F0(int N, int i, const IntegerVector::const_iterator& row_indices_start, const IntegerVector::const_iterator& row_indices_end) {
  // IntegerVector row_indices(row_indices_start, row_indices_end);
  
  int Nm1 = row_indices_end - row_indices_start;
  int N1 = N - Nm1;
  
  IntegerVector col_indices = seq(0, N1 - 2);
  for(int j = i; j < N1 - 1; j++) col_indices(j) += 1;
  
  IntegerVector indices((N1 - 1) * Nm1);
  IntegerVector::const_iterator ri = row_indices_start;
  IntegerVector::iterator ci = col_indices.begin();
  IntegerVector::iterator mi = indices.begin();
  
  for(; ci != col_indices.end(); ++ci){
    ri = row_indices_start;
    for(; ri != row_indices_end; ++ri, ++mi){
      *mi = *ci * N + *ri;
    }
  }
  // Rcout << "F0, sampling from: " << indices << "\n";
  return sample(indices, 1, true)(0);
}
int  sample_Fk(int N, int i, const IntegerVector& k_indices, const IntegerVector::iterator& sim_vector_start, const IntegerVector::iterator& sim_vector_end) {
  // Rcout << "Drawing from Fk\n";
  IntegerVector sim_vector(sim_vector_start, sim_vector_end);
  
  IntegerVector col_indices = seq_len(k_indices.size());
  IntegerVector row_entries(sim_vector.size());//  / N;
  IntegerVector::iterator r = row_entries.begin();
  IntegerVector::iterator s = sim_vector_start;
  for(; r != row_entries.end(); ++r, ++s) {
    *r = *s % N;
  }
  
  int out;
  int col_ind = sample(col_indices, 1)(0) - 1;
  // Rcout << "Fk, sampling column: " << col_indices << "\n";
  
  if(col_ind == i) {
    IntegerVector all_indices = seq(-1, N - 1);
    IntegerVector row_indices = setdiff(all_indices, k_indices);
    row_indices.sort();
    
    out = sample_F0(N, i, row_indices.cbegin() + 1, row_indices.cend());
  }
  else {
    int row_ind = sample(row_entries, 1)(0);
    // Rcout << "Fk, sampling row: " << row_entries << "\n";
    
    out = col_ind * N + row_ind;
  }
  // int row_ind = sample()
  // Rcout << "Fk, output index: " << out << "\n";
  
  return out;
  
}
int sample_Fi(int N, int i, const IntegerVector& k_indices) {
  
  IntegerVector all_indices = seq(-1, N - 1);
  IntegerVector row_indices = setdiff(all_indices, k_indices);
  row_indices.sort();
  
  // Rcout << "Fi, sampling from: " << row_indices << "\n";
  int out = sample(row_indices, 1, true)(0);
  
  
  if(out < 0) out = sample_F0(N, i, row_indices.begin() + 1, row_indices.end());  // these are in inverse order (e.g. 5 4 3 2 1 0 -1) - subtracting 1 from the end removes the "-1"
  else out = out + N * i;
  return out;
}
int sample_Fik(int N, int i, const IntegerVector& k_indices, const IntegerVector::iterator& sim_vector_start, const IntegerVector::iterator& sim_vector_end) {
  IntegerVector sim_vector(sim_vector_start, sim_vector_end);
  // Rcout << "Fik, sampling from: " << sim_vector << "\n";
  int out = sample(sim_vector, 1, true)(0);
  if(out < 0) out = sample_Fk(N, i, k_indices, sim_vector_start, sim_vector_end - 1);
  // else out = out + N * i;
  return out;
}

// // [[Rcpp::export]]
// int sample_Fik_wrapper(int N, int i, IntegerVector& k_indices, IntegerVector sim_vector) {
//   int k = sample_Fik(N, i, k_indices, sim_vector.begin(), sim_vector.end());
//   return k;
// }

/*** R
# This is R code
N <- 10
N1 <- 3
i <- 1
k_ind <- c(0, 5, 6)
sim_vector <- c(c(7:9, 2) + N * i, -1)

# sample_Fik_wrapper(N, i, sim_vector, k_ind)
tmp <- replicate(10, sample_Fik_wrapper(N, i, k_ind, sim_vector))
table(tmp)
hist(tmp, breaks = (-1:max(tmp)) + 0.5)
# tmp %>% table() %>% prop.table()
*/
