#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "mvrnorm.h"

// [[Rcpp::export]]
arma::mat sim_by_instance(arma::uword n, const arma::uvec& Nk, const arma::vec& s, const arma::mat& rho){

  arma::uword K = Nk.n_elem;
  arma::uword N = arma::sum(Nk);
  double val;

  arma::uvec K_range = arma::cumsum(Nk);
  arma::mat V(N, N);

  // generate variance-covariance matrix:
  arma::uword k_start = 0, j_start = 0;
  for(arma::uword k = 0; k < K; k ++) {
    if(k == 0) k_start = 0;
    else k_start = K_range(k - 1);

    for(arma::uword j = 0; j <= k; j++) {
      if(j == k){
        val = rho(k,k) * s(k) * s(k);
        V.submat(k_start, k_start, arma::size(Nk(k), Nk(k))).fill(val);
      }
      else {
        if(j == 0) j_start = 0;
        else j_start = K_range(j - 1);

        val = s(k) * s(j) * rho(k, j);

        V.submat(k_start, j_start, arma::size(Nk(k), Nk(j))).fill(val);
        V.submat(j_start, k_start, arma::size(Nk(j), Nk(k))).fill(val);
      }

    }
  }
  V.diag().ones();

  arma::vec mu(N, arma::fill::zeros);

  arma::mat X = mvrnorm(n, mu, V);

return X;

}

