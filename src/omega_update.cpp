#include "RcppArmadillo.h"
#include "hotspotRM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec omega_update(arma::mat X_star,
                       arma::vec m,
                       arma::vec theta,
                       arma::vec alpha,
                       arma::vec log_sum_exp){

arma::vec psi = X_star*theta + 
                alpha;

arma::vec omega = rcpp_pgdraw(m,
                              (psi - log_sum_exp));

return(omega);

}
































































