#include "RcppArmadillo.h"
#include "hotspotRM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List hotspotRM1(int mcmc_samples,
                      arma::mat z,
                      arma::vec n,
                      arma::vec m,
                      arma::mat v,
                      arma::mat R,
                      double metrop_sd_rho_trans,
                      Rcpp::Nullable<double> sigma2_gamma_prior = R_NilValue,
                      Rcpp::Nullable<arma::mat> Sigma_inv_scale_prior = R_NilValue,
                      Rcpp::Nullable<double> Sigma_inv_df_prior = R_NilValue,
                      Rcpp::Nullable<arma::mat> beta_init = R_NilValue,
                      Rcpp::Nullable<arma::mat> gamma_init = R_NilValue,
                      Rcpp::Nullable<arma::mat> Sigma_inv_init = R_NilValue,
                      Rcpp::Nullable<double> rho_init = R_NilValue){
 
//Defining Parameters and Quantities of Interest
int d = z.n_cols;
int c = n.size();
int p_v = v.n_cols;
int sum_n = sum(n);

//kappa
arma::mat M = arma::repmat(m, 1, (d-1));
arma::mat kappa = z.cols(1, (d-1)) +
                  -0.50*M;
       
//x           
arma::mat x(sum_n, c); x.fill(0.00);
int counter = 0;
for(int j = 0; j < c; ++j){
  
   int start = counter;
   int end = counter +
             n[j] +
             -1;
   x(arma::span(start, end), j).ones();
   counter = counter +
             n[j];
  
   }
  
//v_design
Rcpp::List v_design_list(d-1);
for(int j = 0; j < (d-1); ++j){
  
   arma::mat block(c, (p_v*c), arma::fill::zeros);
   arma::rowvec vrow = v.row(j);  
  
   for(int k = 0; k < c; ++k){
      block.submat(k, (k*p_v), k, (k*p_v + p_v - 1)) = vrow;
      }
   v_design_list[j] = block;

   }

int total_rows = (d-1)*c;
arma::mat v_design(total_rows, (p_v*c), arma::fill::zeros);
for(int j = 0; j < (d-1); ++j){
  
   arma::mat block = Rcpp::as<arma::mat>(v_design_list[j]);
   v_design.rows((j*c), ((j + 1)*c - 1)) = block;

   }

//Parameters
arma::cube beta(c, (d-1), mcmc_samples); beta.fill(0.00);
arma::cube gamma(p_v, c, mcmc_samples); gamma.fill(0.00);
arma::vec gamma_long(p_v*c); gamma_long.fill(0.00);
arma::mat alpha(sum_n, (d-1)); alpha.fill(0.00);
arma::cube Sigma_inv(c, c, mcmc_samples); Sigma_inv.fill(0.00);
arma::vec rho(mcmc_samples); rho.fill(0.00);
arma::mat omega(sum_n, (d-1)); omega.fill(0.00);
arma::cube p_signal(c, d, mcmc_samples); p_signal.fill(0.00);

//Prior Information
double sigma2_gamma = 10000.00;
if(sigma2_gamma_prior.isNotNull()){
  sigma2_gamma = Rcpp::as<double>(sigma2_gamma_prior);
  }
arma::mat gamma_prior_cov_inv = arma::eye((c*p_v), (c*p_v))/sigma2_gamma;

arma::mat Sigma_inv_scale(c, c); Sigma_inv_scale.eye();
if(Sigma_inv_scale_prior.isNotNull()){
  Sigma_inv_scale = Rcpp::as<arma::mat>(Sigma_inv_scale_prior);
  }
arma::mat Sigma_inv_scale_inv = inv_sympd(Sigma_inv_scale);

double Sigma_inv_df = c +
                      1.00;
if(Sigma_inv_df_prior.isNotNull()){
  Sigma_inv_df = Rcpp::as<double>(Sigma_inv_df_prior);
  }

//Initial Values
beta.slice(0).zeros();
if(beta_init.isNotNull()){
  beta.slice(0) = Rcpp::as<arma::mat>(beta_init);
  }

gamma.slice(0).fill(0.00);
if(gamma_init.isNotNull()){
  gamma.slice(0) = Rcpp::as<arma::mat>(gamma_init);
  }
gamma_long = arma::vectorise(gamma.slice(0));
arma::vec prod = v_design*gamma_long;
arma::mat beta_prior_mean = arma::reshape(prod, 
                                          c, 
                                          (d-1));

arma::mat temp_mat(sum_n, d); temp_mat.fill(0.00);
for(int j = 0; j < (d-1); ++j){
  temp_mat.col(j + 1) = x*beta.slice(0).col(j) +
                        alpha.col(j);
  }

Sigma_inv.slice(0) = arma::eye(c, c);
if(Sigma_inv_init.isNotNull()){
  Sigma_inv.slice(0) = Rcpp::as<arma::mat>(Sigma_inv_init);
  }

rho(0) = 0.50;
if(rho_init.isNotNull()){
  rho(0) = Rcpp::as<double>(rho_init);
  }
arma::mat Q_piece = arma::diagmat(arma::sum(R, 1)) +
                    -R;
arma::mat Q = rho(0)*Q_piece +
              (1.00 - rho(0))*eye((d-1), (d-1));
double log_deter = 0.00; 
double sign = 0.00;     
log_det(log_deter, sign, Q);

//p_signal
arma::vec denom = exp(beta.slice(0).col(0));
for(int j = 1; j < (d-1); ++j){
   denom = denom +
           exp(beta.slice(0).col(j));
   }
p_signal.slice(0).col(0) = 1.00/(1.00 + denom);
for(int j = 0; j < (d-1); ++j){
   p_signal.slice(0).col(j+1) = exp(beta.slice(0).col(j))/(1.00 + denom);
   }

//Metropolis Settings
int acctot_rho_trans = 1;

//Main Sampling Loop
for(int j = 1; j < mcmc_samples; ++j){
    
   beta.slice(j) = beta.slice(j-1);
   for(int k = 0; k < (d-1); ++k){
     
      //log_sum_exp
      arma::mat temp_mat_reduced = temp_mat;
      temp_mat_reduced.shed_col(k+1);
     
      arma::vec log_sum_exp(temp_mat_reduced.n_rows);
      for(arma::uword l = 0; l < temp_mat_reduced.n_rows; ++l){
        
         arma::rowvec row = temp_mat_reduced.row(l);
         double max_val = row.max();
         log_sum_exp(l) = max_val + 
                          log(sum(exp(row - max_val)));
     
         }
     
      //omega
      omega.col(k) = omega_update(x,
                                  m,
                                  beta.slice(j).col(k),
                                  alpha.col(k),
                                  log_sum_exp);
      
      //beta
      beta.slice(j).col(k) = beta_update(x,
                                         R.col(k),
                                         c,
                                         sum_n,
                                         sum(R.col(k)),
                                         kappa.col(k),
                                         omega.col(k),
                                         beta.slice(j),
                                         alpha.col(k),
                                         Sigma_inv.slice(j-1),
                                         rho(j-1),
                                         log_sum_exp,
                                         beta_prior_mean,
                                         beta_prior_mean.col(k));
      temp_mat.col(k+1) = x*beta.slice(j).col(k) +
                          alpha.col(k);
      
      }
   
   //Centering on the fly
   for(int k = 0; k < c; ++k){
     
      double spatial_mean = mean(beta.slice(j).row(k) - beta_prior_mean.row(k));
      beta.slice(j).row(k) = beta.slice(j).row(k) + 
                             -spatial_mean;
       
      }
   for(int k = 0; k < (d-1); ++k){
      temp_mat.col(k+1) = x*beta.slice(j).col(k);
      }
   
   //gamma
   Rcpp::List gamma_output = gamma_update(v_design,
                                          c,
                                          p_v,
                                          beta.slice(j),
                                          Sigma_inv.slice(j-1),
                                          Q,
                                          gamma_prior_cov_inv);
   gamma.slice(j) = Rcpp::as<arma::mat>(gamma_output[0]);
   gamma_long = Rcpp::as<arma::vec>(gamma_output[1]);
   arma::vec prod = v_design*gamma_long;
   beta_prior_mean = arma::reshape(prod,  
                                   c, 
                                   (d-1));
   
   //Sigma_inv
   Sigma_inv.slice(j) = Sigma_inv_update(c,
                                         d,
                                         v_design_list,
                                         beta.slice(j),
                                         gamma_long,
                                         Q,
                                         Sigma_inv_scale_inv,
                                         Sigma_inv_df);
   
   //rho
   Rcpp::List rho_output = rho_update(c,
                                      d,
                                      R,
                                      beta.slice(j),
                                      Sigma_inv.slice(j),
                                      rho(j-1),
                                      Q,
                                      log_deter,
                                      beta_prior_mean,
                                      metrop_sd_rho_trans,
                                      acctot_rho_trans);
   rho(j) = Rcpp::as<double>(rho_output[0]);
   acctot_rho_trans = Rcpp::as<int>(rho_output[1]);
   Q = Rcpp::as<arma::mat>(rho_output[2]);
   log_deter = Rcpp::as<double>(rho_output[3]);
   
   //p_signal
   arma::vec denom = exp(beta.slice(j).col(0));
   for(int k = 1; k < (d-1); ++k){
      denom = denom +
              exp(beta.slice(j).col(k));
      }
   p_signal.slice(j).col(0) = 1.00/(1.00 + denom);
   for(int k = 0; k < (d-1); ++k){
      p_signal.slice(j).col(k+1) = exp(beta.slice(j).col(k))/(1.00 + denom);
      }
   
   //Progress
   if((j + 1) % 10 == 0){ 
     Rcpp::checkUserInterrupt();
     }
    
   if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
      
      double completion = round(100*((j + 1)/(double)mcmc_samples));
      Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
      
      double accrate_rho_trans = round(100*(acctot_rho_trans/(double)j));
      Rcpp::Rcout << "rho Acceptance: " << accrate_rho_trans << "%" << std::endl;
      
      Rcpp::Rcout << "********************" << std::endl;
      
      }
    
     }
  
   return Rcpp::List::create(Rcpp::Named("beta") = beta,
                             Rcpp::Named("gamma") = gamma,
                             Rcpp::Named("Sigma_inv") = Sigma_inv,
                             Rcpp::Named("rho") = rho,
                             Rcpp::Named("p_signal") = p_signal,
                             Rcpp::Named("acctot_rho_trans") = acctot_rho_trans);
  
   }


