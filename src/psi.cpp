#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector psi(const NumericVector& x, const NumericVector& y) {
  int Nx = x.size();
  int Ny = y.size();
  
  double inc_x = 1. / Nx;
  double inc_y = 1. / Ny;
  double fin = inc_x * inc_y;
  
  arma::vec z = join_cols(as<arma::vec>(x), as<arma::vec>(y));
  arma::vec cat = join_cols(arma::zeros(Nx), arma::ones(Ny));
  arma::uvec indices = sort_index(z);
  arma::uvec::iterator ind_ptr = indices.end();
  
  double pF = 0;
  double pG = 1;

  do {  
    ind_ptr--;
    if(cat(*ind_ptr) < .5) {
      pG = pG - inc_x;
    } else {
      pF = pF + inc_y;
    }
// 
//     Rcout << "pG = " << pG << "\n";
//     Rcout << "pF = " << pF << "\n";
//     Rcout << "index = " << *ind_ptr << " out of " << z.size() << "\n";
//     Rcout << "t = " << z(*ind_ptr) << "\n\n";
  }
  while(pG - pF > fin);
  
  NumericVector result = NumericVector::create(_["t"] = (z(*ind_ptr) + z(*(ind_ptr + 1))) / 2, _["tau"] = pG);

  return result;
}


/*** R

psi <-function(D1, D2){
  F1<-stats::ecdf(D1)
  F2<-stats::ecdf(D2)
  N<-4000
  tt<-seq(min(D1, na.rm = TRUE),max(D2, na.rm = TRUE), length=N)
  
  delta <- abs(1 - F2(tt) - F1(tt))
  tc <-  stats::median(tt[which(delta <= delta[which.min(delta)] + 10e-8)])
  c(
    t = tc,
    tau = F1(tc)
  )
}
x <- c(1:2, 4, 8)
y <- c(4:6)
  
Gx <- ecdf(x)
Fy <- ecdf(y)

nx <- length(x)
ny <- length(y)
  psi_cpp(x, y)
  psi(x, y)

  plot(x, (1:nx)/nx, 
    type = 's',
    ylim = c(0, 1),
    col = 'red'
  )
  lines(y, (ny:1)/ny, 
    type = 'S',
    ylim = c(0, 1),
    col = 'blue'
  )
  
    
  */
