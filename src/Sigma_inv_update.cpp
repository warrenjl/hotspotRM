#include "RcppArmadillo.h"
#include "hotspotRM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat Sigma_inv_update(int c,
                           int d,
                           Rcpp::List v_design_list,
                           arma::mat beta,
                           arma::vec gamma_long,
                           arma::mat Q,
                           arma::mat Sigma_inv_scale_inv,
                           double Sigma_inv_df){

arma::mat temp(c,c); temp.fill(0.00);
for(int j = 0; j < (d-1); ++ j){
   for(int k = 0; k < (d-1); ++k){
      temp = temp +
             ((beta.col(j) - (Rcpp::as<arma::mat>(v_design_list[j])*gamma_long))*trans(beta.col(k) - (Rcpp::as<arma::mat>(v_design_list[k])*gamma_long)))*Q(k,j);
      }
   } 
temp = temp +
       Sigma_inv_scale_inv;
temp = inv_sympd(temp);

double df = d + 
            -1 +
            Sigma_inv_df;

//Bartlett Decomposition
arma::mat L = arma::chol(temp,
                         "lower");
arma::mat A(c,c); A.fill(0.00);
for(int j = 0; j < c; ++j){
  
   A(j,j) = sqrt(R::rchisq(df - j));
   for(int k = 0; k < j; ++k){
      A(j,k) = R::rnorm(0.00,
                        sqrt(1.00));
      }
  
   }

arma::mat Sigma_inv = (L*A)*trans(L*A);

return(Sigma_inv);

}

  
  

  



