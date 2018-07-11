#include <Rcpp.h>
using namespace Rcpp;

int sample_F0(int N, int i, const IntegerVector::iterator& row_indices_start, const IntegerVector::iterator& row_indices_end) {
  // IntegerVector row_indices(row_indices_start, row_indices_end);

  int Nm1 = row_indices_end - row_indices_start;
  int N1 = N - Nm1;

  IntegerVector col_indices = seq(0, N1 - 2);
  for(int j = i; j < N1 - 1; j++) col_indices(j) += 1;

  IntegerVector indices((N1 - 1) * Nm1);
  IntegerVector::iterator ri = row_indices_start;
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
int sample_Fi(int N, int i, const IntegerVector& k_indices) {

  IntegerVector all_indices = seq(-1, N - 1);
  IntegerVector row_indices = setdiff(all_indices, k_indices);
  row_indices = row_indices;

  // Rcout << "Fi, sampling from: " << row_indices << "\n";
  int out = sample(row_indices, 1, true)(0);

  if(out < 0) out = sample_F0(N, i, row_indices.begin(), row_indices.end() - 1);
  // else out = out + N * i;
  return out;
}
int sample_Fik(int N, int i, const IntegerVector::iterator& sim_vector_start, const IntegerVector::iterator& sim_vector_end, const IntegerVector& k_indices) {
  IntegerVector sim_vector(sim_vector_start, sim_vector_end);
  // Rcout << "Fik, sampling from: " << sim_vector << "\n";
  int out = sample(sim_vector, 1, true)(0);
  if(out < 0) out = sample_Fi(N, i, k_indices);
  // else out = out + N * i;
  return out;
}
