#include "RcppArmadillo.h"
#include "MLDM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::cube psi2_create(arma::cube p_signal){

//p_signal_space
int c = p_signal.n_rows;
int d = p_signal.n_cols;
int mcmc_samples = p_signal.n_slices;

//psi2
arma::cube psi2(c, d, mcmc_samples); psi2.fill(0.00);
for(int j = 0; j < mcmc_samples; ++j){
  
   arma::vec region_avg = arma::mean(p_signal.slice(j), 1);    
   psi2.slice(j) = p_signal.slice(j) - arma::repmat(region_avg, 1, d);
   
   }

return(psi2);

}





