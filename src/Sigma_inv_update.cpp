#include "RcppArmadillo.h"
#include "hotspotRM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat Sigma_inv_update(int c,
                           int d,
                           int p_g,
                           Rcpp::List v_design_list,
                           arma::mat theta,
                           arma::vec gamma_long,
                           arma::mat Q,
                           arma::mat Sigma_inv_scale_inv,
                           double Sigma_inv_df){

arma::mat B(c*(1 + p_g), (d-1));
for(int j = 0; j < (d-1); ++j){
   B.col(j) = theta.col(j) + 
              -Rcpp::as<arma::mat>(v_design_list[j])*gamma_long;
   }

arma::mat temp = inv_sympd(B*Q*B.t() + Sigma_inv_scale_inv);

double df = d + 
            -1 +
            Sigma_inv_df;

//Bartlett Decomposition
arma::mat L = arma::chol(temp,
                         "lower");
arma::mat A((c*(1 + p_g)), (c*(1 + p_g))); A.fill(0.00);
for(int j = 0; j < (c*(1 + p_g)); ++j){
  
   A(j,j) = sqrt(R::rchisq(df - j));
   for(int k = 0; k < j; ++k){
      A(j,k) = R::rnorm(0.00,
                        sqrt(1.00));
      }
  
   }

arma::mat Sigma_inv = (L*A)*trans(L*A);

return(Sigma_inv);

}


