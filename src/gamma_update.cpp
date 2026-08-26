#include "RcppArmadillo.h"
#include "MLDM.h"
using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List gamma_update(Rcpp::List v_design_list,
                        int c,
                        int d,
                        int p_v,
                        int p_g,
                        arma::mat theta,
                        arma::mat Sigma_inv_old,
                        arma::mat Q,
                        arma::mat gamma_prior_cov_inv){

int c_star = c*(1 + p_g);
int p_star = p_v*c_star;

arma::mat cov_gamma_inv(p_star, p_star, arma::fill::zeros);
arma::vec mean_gamma_unscaled(p_star, arma::fill::zeros);

for(int j = 0; j < (d-1); ++j){

   arma::mat Vj = Rcpp::as<arma::mat>(v_design_list[j]);
   arma::mat VjT_SigInv = Vj.t()*Sigma_inv_old;
   for(int k = 0; k < (d-1); ++k){
      arma::mat Vk = Rcpp::as<arma::mat>(v_design_list[k]);
      cov_gamma_inv += Q(j,k) * (VjT_SigInv * Vk);
      }

   arma::vec Qrow_theta(c_star, arma::fill::zeros);
   for(int k = 0; k < (d-1); ++k){
      Qrow_theta += Q(j,k) * theta.col(k);
      }
   mean_gamma_unscaled += VjT_SigInv * Qrow_theta;

   }

arma::mat cov_gamma = inv_sympd(cov_gamma_inv + gamma_prior_cov_inv);
arma::vec mean_gamma = cov_gamma * mean_gamma_unscaled;

arma::mat ind_norms = arma::randn(1, p_star);
arma::vec gamma_long = mean_gamma + 
                       trans(ind_norms*arma::chol(cov_gamma));

arma::mat gamma = arma::reshape(gamma_long, 
                                p_v, 
                                c_star);

return Rcpp::List::create(Rcpp::Named("gamma") = gamma,
                          Rcpp::Named("gamma_long") = gamma_long);

}

