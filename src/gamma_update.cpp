#include "RcppArmadillo.h"
#include "hotspotRM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List gamma_update(arma::mat v_design,
                        int c,
                        int p_v,
                        arma::mat beta,
                        arma::mat Sigma_inv_old,
                        arma::mat Q,
                        arma::mat gamma_prior_cov_inv){
  
arma::mat Q_Sigma_inv = kron(Q, Sigma_inv_old);
arma::mat v_design_trans = trans(v_design);

arma::mat cov_gamma = inv_sympd(v_design_trans*(Q_Sigma_inv*v_design) +
                                gamma_prior_cov_inv);

arma::vec mean_gamma = cov_gamma*(v_design_trans*(Q_Sigma_inv*arma::vectorise(beta)));

arma::mat ind_norms = arma::randn(1, (c*p_v));
arma::vec gamma_long = mean_gamma + 
                       trans(ind_norms*arma::chol(cov_gamma));

arma::mat gamma = arma::reshape(gamma_long, 
                                p_v, 
                                c);

return Rcpp::List::create(Rcpp::Named("gamma") = gamma,
                          Rcpp::Named("gamma_long") = gamma_long);
  
}



