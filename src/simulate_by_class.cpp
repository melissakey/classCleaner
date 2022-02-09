#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "mvrnorm.h"

using namespace Rcpp;

//' @title tmp function
//' @export
// [[Rcpp::export]]
arma::mat sim_by_class( arma::uword n, const arma::uvec&  Nk, const arma::colvec& s, double tau, const arma::mat& rho){

  arma::uword K = Nk.n_elem;
  arma::uword N = sum(Nk);
  arma::mat V(K, K); // = s * s.t();
  // S.diag() = s;

  if(rho.n_elem == 1){
    V.fill(arma::as_scalar(rho));
    V.diag().ones();
    
  }
  else {
    V = rho;
  }

  V = V % (s * s.t());

  arma::vec mu0(K, arma::fill::zeros);
  arma::mat W = mvrnorm(n,mu0, V);


  arma::mat X = arma::randn<arma::mat>(n, N) * tau;

  arma::uword j = 0;
  arma::uword k = 0;
  for(arma::uword i = 0; i < N; ++i, ++j){
    if(j == Nk[k]) {
      ++k;
      j = 0;
    }
    X.col(i) += W.col(k);

  }

  return X;

}

/***R
s <- rep(1, 10)
# Nk <- c(40, 200)
Nk <- rep(5, 10)
n <- 10
rho <- matrix(c(.6, .1, .1, .25), 2, 2)
rho <- matrix(.5, ncol = 10, nrow = 10)
diag(rho) <- 1

# (tmp <- sim_by_instance(n, Nk, s, rho))
# tmp2 <- classCleaner:::sim_by_class(n, Nk, s, 1, rho = rho)
tmp3 <- sim_by_class(n, Nk, s, 1, rho = rho)
*/
