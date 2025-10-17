#include "RcppArmadillo.h"
#include "hotspotRM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec alpha_update(arma::mat x,
                       int sum_n,
                       arma::vec kappa,
                       arma::vec omega,
                       arma::vec beta,
                       arma::vec log_sum_exp,
                       arma::mat alpha_prior_cov_inv){
  
arma::mat cov_alpha = inv_sympd(arma::diagmat(omega) +
                                alpha_prior_cov_inv);
  
arma::vec mean_alpha = cov_alpha*(kappa + omega%(log_sum_exp - x*beta));
  
arma::mat ind_norms = arma::randn(1, sum_n);
arma::vec alpha = mean_alpha + 
                  trans(ind_norms*arma::chol(cov_alpha));
  
return(alpha);
  
}


