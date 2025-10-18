#ifndef __hotspotRM__
#define __hotspotRM__

arma::vec rcpp_pgdraw(arma::vec b, 
                      arma::vec c);

arma::vec omega_update(arma::mat x,
                       arma::vec m,
                       arma::vec beta,
                       arma::vec alpha,
                       arma::vec log_sum_exp);

arma::vec beta_update(arma::mat x,
                      arma::vec R,
                      int c,
                      int sum_n,
                      int sum_R,
                      arma::vec kappa,
                      arma::vec omega,
                      arma::mat beta,
                      arma::vec alpha,
                      arma::mat Sigma_inv_old,
                      double rho_old,
                      arma::vec log_sum_exp,
                      arma::mat beta_prior_mean_mat,
                      arma::vec beta_prior_mean_vec);

arma::vec alpha_update(arma::mat x,
                       int sum_n,
                       arma::vec kappa,
                       arma::vec omega,
                       arma::vec beta,
                       arma::vec log_sum_exp,
                       arma::mat alpha_prior_cov_inv);

Rcpp::List gamma_update(arma::mat v_design,
                        int c,
                        int p_v,
                        arma::mat beta,
                        arma::mat Sigma_inv_old,
                        arma::mat Q,
                        arma::mat gamma_prior_cov_inv);

double sigma2_update(int k,
                     int counter,
                     int n,
                     int d,
                     arma::mat alpha,
                     double shape_update,
                     double b_sigma2);

arma::mat Sigma_inv_update(int c,
                           int d,
                           Rcpp::List v_design_list,
                           arma::mat beta,
                           arma::vec gamma_long,
                           arma::mat Q,
                           arma::mat Sigma_inv_scale_inv,
                           double Sigma_inv_df);

Rcpp::List rho_update(int c,
                      int d,
                      arma::mat R,
                      arma::mat beta,
                      arma::mat Sigma_inv,
                      double rho_old,
                      arma::mat Q,
                      double log_deter,
                      arma::mat beta_prior_mean_mat,
                      double metrop_sd_rho_trans,
                      int acctot_rho_trans);

Rcpp::List hotspotRM(int mcmc_samples,
                     arma::mat z,
                     arma::vec n,
                     arma::vec m,
                     arma::mat v,
                     arma::mat R,
                     double metrop_sd_rho_trans,
                     Rcpp::Nullable<double> sigma2_gamma_prior,
                     Rcpp::Nullable<double> a_sigma2_prior,
                     Rcpp::Nullable<double> b_sigma2_prior,
                     Rcpp::Nullable<arma::mat> Sigma_inv_scale_prior,
                     Rcpp::Nullable<double> Sigma_inv_df_prior,
                     Rcpp::Nullable<arma::mat> beta_init,
                     Rcpp::Nullable<arma::mat> gamma_init,
                     Rcpp::Nullable<arma::mat> alpha_init,
                     Rcpp::Nullable<arma::vec> sigma2_init,
                     Rcpp::Nullable<arma::mat> Sigma_inv_init,
                     Rcpp::Nullable<double> rho_init);

#endif // __hotspotRM__
