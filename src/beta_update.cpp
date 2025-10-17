#include "RcppArmadillo.h"
#include "hotspotRM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec beta_update(arma::mat x,
                      arma::vec R,
                      int c,
                      int sum_n,
                      int sum_R,
                      arma::vec kappa,
                      arma::vec omega,
                      arma::mat beta,
                      arma::vec alpha,
                      arma::mat Sigma_inv_old,
                      double rho_old,
                      arma::vec log_sum_exp,
                      arma::mat beta_prior_mean_mat,
                      arma::vec beta_prior_mean_vec){
  
arma::mat omega_mat(sum_n, c);
for(int j = 0; j < c; ++j){
   omega_mat.col(j) = omega;
   }

arma::mat x_trans = trans(x);
  
arma::mat cov_beta = inv_sympd(x_trans*(omega_mat%x) +
                               (Sigma_inv_old*(rho_old*sum_R + 1.00 - rho_old)));

arma::vec mean_beta = cov_beta*(x_trans*(kappa + omega%(log_sum_exp - alpha)) +
                                rho_old*(Sigma_inv_old*((beta - beta_prior_mean_mat)*R)) + 
                                (rho_old*sum_R + 1.00 - rho_old)*(Sigma_inv_old*beta_prior_mean_vec));
  
arma::mat ind_norms = arma::randn(1, c);
arma::vec beta_new = mean_beta + 
                     trans(ind_norms*arma::chol(cov_beta));
  
return(beta_new);

}



