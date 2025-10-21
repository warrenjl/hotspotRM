#include "RcppArmadillo.h"
#include "hotspotRM.h"
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

double rate_update = 0.00;
for(int l = 0; l < (d-1); ++l){
   rate_update = rate_update +
                 arma::as_scalar(sum(pow(alpha(arma::span(counter, (n + counter - 1)), l), 2)));
   }
rate_update = rate_update/2.00 +
              b_sigma2;
double sigma2 = 1.00/R::rgamma(shape_update,
                               (1.00/rate_update));

return(sigma2);

}





