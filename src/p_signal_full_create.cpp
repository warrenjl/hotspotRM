#include "RcppArmadillo.h"
#include "MLDM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::cube p_signal_full_create(arma::cube beta,
                                arma::cube gamma,
                                arma::mat v){
//p_signal_full
int c = beta.n_rows;
int d = beta.n_cols + 1;
int mcmc_samples = beta.n_slices;

arma::cube p_signal_full(c, d, mcmc_samples); p_signal_full.fill(0.00);
for(int j = 0; j < mcmc_samples; ++j){

   arma::mat gamma_beta = gamma.slice(j).cols(0, c - 1);   
   arma::mat correction  = v * gamma_beta;                  
   arma::mat beta_full  = beta.slice(j) - correction.t();  

   arma::vec denom = exp(beta_full.col(0));
   for(int k = 1; k < (d-1); ++k){
      denom = denom +
              exp(beta_full.col(k));
      }
   p_signal_full.slice(j).col(0) = 1.00/(1.00 + denom);
   for(int k = 0; k < (d-1); ++k){
      p_signal_full.slice(j).col(k+1) = exp(beta_full.col(k))/(1.00 + denom);
      }
   }

return(p_signal_full);

}



