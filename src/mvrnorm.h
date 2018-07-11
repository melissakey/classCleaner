#ifndef MVRNORM_H
#define MVRNORM_H

#include <RcppArmadillo.h>

arma::mat mvrnorm(int n, arma::vec mu, arma::mat sigma);

#endif
