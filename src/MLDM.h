#ifndef __MLDM__
#define __MLDM__

arma::cube p_signal_space_create(arma::cube beta);

arma::cube p_signal_full_create(arma::cube beta,
                                arma::cube gamma,
                                arma::mat v);

arma::cube psi1_create(arma::cube p_signal);

arma::cube psi2_create(arma::cube p_signal);
  
arma::vec rcpp_pgdraw(arma::vec b, 
                      arma::vec c);

arma::vec omega_update(arma::mat X_star,
                       arma::vec m,
                       arma::vec theta,
                       arma::vec alpha,
                       arma::vec log_sum_exp);

arma::vec theta_update(arma::mat X_star,
                       arma::vec R,
                       int c,
                       int p_g,
                       arma::vec n,
                       int sum_R,
                       arma::vec kappa,
                       arma::vec omega,
                       arma::mat theta,
                       arma::vec alpha,
                       arma::mat Sigma_inv_old,
                       double rho_old,
                       arma::vec log_sum_exp,
                       arma::mat theta_prior_mean_mat,
                       arma::vec theta_prior_mean_vec);

arma::vec alpha_update(arma::mat X_star,
                       int sum_n,
                       arma::vec kappa,
                       arma::vec omega,
                       arma::vec theta,
                       arma::vec log_sum_exp,
                       arma::mat alpha_prior_cov_inv);

Rcpp::List gamma_update(Rcpp::List v_design_list,
                        int c,
                        int d,
                        int p_v,
                        int p_g,
                        arma::mat theta,
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
                           int p_g,
                           Rcpp::List v_design_list,
                           arma::mat theta,
                           arma::vec gamma_long,
                           arma::mat Q,
                           arma::mat Sigma_inv_scale_inv,
                           double Sigma_inv_df);

Rcpp::List rho_update(int c,
                      int d,
                      int p_g,
                      arma::mat R,
                      arma::mat theta,
                      arma::mat Sigma_inv,
                      double rho_old,
                      arma::mat Q,
                      double log_deter,
                      arma::mat theta_prior_mean_mat,
                      double metrop_sd_rho_trans,
                      int acctot_rho_trans);

Rcpp::List MLDM1(int mcmc_samples,
                 arma::mat z,
                 arma::vec n,
                 arma::vec m,
                 arma::mat v,
                 arma::mat R,
                 arma::mat x,
                 double metrop_sd_rho_trans,
                 Rcpp::Nullable<double> sigma2_gamma_prior,
                 Rcpp::Nullable<arma::mat> Sigma_inv_scale_prior,
                 Rcpp::Nullable<double> Sigma_inv_df_prior,
                 Rcpp::Nullable<arma::mat> theta_init,
                 Rcpp::Nullable<arma::mat> gamma_init,
                 Rcpp::Nullable<arma::mat> Sigma_inv_init,
                 Rcpp::Nullable<double> rho_init);
  
Rcpp::List MLDM(int mcmc_samples,
                arma::mat z,
                arma::vec n,
                arma::vec m,
                arma::mat v,
                arma::mat R,
                arma::mat x,
                double metrop_sd_rho_trans,
                Rcpp::Nullable<double> sigma2_gamma_prior,
                Rcpp::Nullable<double> a_sigma2_prior,
                Rcpp::Nullable<double> b_sigma2_prior,
                Rcpp::Nullable<arma::mat> Sigma_inv_scale_prior,
                Rcpp::Nullable<double> Sigma_inv_df_prior,
                Rcpp::Nullable<arma::mat> theta_init,
                Rcpp::Nullable<arma::mat> gamma_init,
                Rcpp::Nullable<arma::mat> alpha_init,
                Rcpp::Nullable<arma::vec> sigma2_init,
                Rcpp::Nullable<arma::mat> Sigma_inv_init,
                Rcpp::Nullable<double> rho_init);

#endif // __MLDM__
