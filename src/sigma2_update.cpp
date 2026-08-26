#include "RcppArmadillo.h"
#include "MLDM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_update(int k,
                     int counter,
                     int n,
                     int d,
                     arma::mat alpha,
                     double shape_update,
                     double b_sigma2){

arma::mat alpha_k = alpha.rows(counter, (counter + n - 1));
double rate_update = arma::accu(alpha_k % alpha_k)/2.00 +
                     b_sigma2;
double sigma2 = 1.00/R::rgamma(shape_update,
                               (1.00/rate_update));

return(sigma2);

}





