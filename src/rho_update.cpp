#include "RcppArmadillo.h"
#include "hotspotRM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List rho_update(int c,
                      int d,
                      arma::mat R,
                      arma::mat beta,
                      arma::mat Sigma_inv,
                      double rho_old,
                      arma::mat Q,
                      double log_deter,
                      arma::mat beta_prior_mean_mat,
                      double metrop_sd_rho_trans,
                      int acctot_rho_trans){

arma::mat Q_old = Q;
double log_deter_old = log_deter;
double rho_trans_old = log(rho_old/(1.00 - rho_old));

double denom = c*0.50*log_deter_old +
               -0.50*dot(arma::vectorise(beta - beta_prior_mean_mat), kron(Q_old, Sigma_inv)*arma::vectorise(beta - beta_prior_mean_mat)) +
               -rho_trans_old +
               -2.00*log(1.00 + exp(-rho_trans_old));
      
double rho_trans = R::rnorm(rho_trans_old, 
                            metrop_sd_rho_trans);
double rho = 1.00/(1.00 + exp(-rho_trans));
arma::mat Q_piece = arma::diagmat(arma::sum(R, 1)) +
                    -R;
Q = rho*Q_piece +
    (1.00 - rho)*eye((d-1), (d-1));
double sign = 0.00;     
log_det(log_deter, 
        sign, 
        Q);
   
double numer = c*0.50*log_deter +
               -0.50*dot(arma::vectorise(beta - beta_prior_mean_mat), kron(Q, Sigma_inv)*arma::vectorise(beta - beta_prior_mean_mat)) +
               -rho_trans +
               -2.00*log(1.00 + exp(-rho_trans));  
      
double ratio = exp(numer - denom);   
double acc = 1;
if(ratio < R::runif(0.00, 1.00)){
        
  Q = Q_old;
  log_deter = log_deter_old;
  rho = rho_old;
  acc = 0;
      
  }
    
acctot_rho_trans = acctot_rho_trans + 
                   acc;

return Rcpp::List::create(Rcpp::Named("rho") = rho,
                          Rcpp::Named("acctot_rho_trans") = acctot_rho_trans,
                          Rcpp::Named("Q") = Q,
                          Rcpp::Named("log_deter") = log_deter);

}



