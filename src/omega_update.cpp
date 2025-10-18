#include "RcppArmadillo.h"
#include "hotspotRM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec omega_update(arma::mat x,
                       arma::vec m,
                       arma::vec beta,
                       arma::vec alpha,
                       arma::vec log_sum_exp){

arma::vec psi = x*beta + 
                alpha;

arma::vec omega = rcpp_pgdraw(m,
                              (psi - log_sum_exp));

return(omega);

}
































































