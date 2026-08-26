#include "RcppArmadillo.h"
#include "MLDM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::cube psi1_create(arma::cube p_signal){

//p_signal_space
int c = p_signal.n_rows;
int d = p_signal.n_cols;
int mcmc_samples = p_signal.n_slices;

//psi1 
arma::cube psi1(c, d, mcmc_samples); psi1.fill(0.00);
for(int j = 0; j < mcmc_samples; ++j){
  
   arma::rowvec group_avg = arma::mean(p_signal.slice(j), 0); 
   psi1.slice(j) = p_signal.slice(j) - arma::repmat(group_avg, c, 1);
  
   
   }

return(psi1);

}





