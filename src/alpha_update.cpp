#include "RcppArmadillo.h"
#include "MLDM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec alpha_update(arma::mat X_star,
                       int sum_n,
                       arma::vec kappa,
                       arma::vec omega,
                       arma::vec theta,
                       arma::vec log_sum_exp,
                       arma::mat alpha_prior_cov_inv){
  
arma::vec diag_vals = alpha_prior_cov_inv.diag();
arma::vec cov_alpha_diag = 1.00/(omega + diag_vals);

arma::vec mean_alpha = cov_alpha_diag%(kappa + 
                                       omega%(log_sum_exp - X_star*theta));

arma::vec alpha = mean_alpha + 
                  sqrt(cov_alpha_diag)%arma::randn(sum_n);
  
return(alpha);
  
}


