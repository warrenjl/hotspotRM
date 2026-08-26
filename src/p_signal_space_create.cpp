#include "RcppArmadillo.h"
#include "MLDM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::cube p_signal_space_create(arma::cube beta){

//p_signal_space
int c = beta.n_rows;
int d = beta.n_cols + 1;
int mcmc_samples = beta.n_slices;
arma::cube p_signal_space(c, d, mcmc_samples); p_signal_space.fill(0.00);
for(int j = 0; j < mcmc_samples; ++j){
  
   arma::vec denom = exp(beta.slice(j).col(0));
   for(int k = 1; k < (d-1); ++k){
      denom = denom +
              exp(beta.slice(j).col(k));
      }
   p_signal_space.slice(j).col(0) = 1.00/(1.00 + denom);
   for(int k = 0; k < (d-1); ++k){
      p_signal_space.slice(j).col(k+1) = exp(beta.slice(j).col(k))/(1.00 + denom);
      }

   }
  
return(p_signal_space);

}





