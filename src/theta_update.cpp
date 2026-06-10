#include "RcppArmadillo.h"
#include "hotspotRM.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

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
                       arma::vec theta_prior_mean_vec){

arma::vec kappa_adj = kappa + 
                      omega%(log_sum_exp - alpha);

arma::mat XtOmegaXstar(c, (c*(1 + p_g)), arma::fill::zeros);
arma::vec x_trans_kappa(c, arma::fill::zeros);
int cnt = 0;
for(int j = 0; j < c; ++j){

   arma::uvec idx = arma::regspace<arma::uvec>(cnt, cnt + n(j) - 1);
   arma::vec omega_j = omega.elem(idx);
   XtOmegaXstar(j, j) = arma::sum(omega_j);
   if(p_g > 0){
      arma::mat X_star_j = arma::mat(X_star.rows(idx));
      XtOmegaXstar.row(j).cols(c, c*(1+p_g)-1) = (omega_j.t() * X_star_j.cols(c, c*(1+p_g)-1));
      }
   x_trans_kappa(j) = arma::sum(kappa_adj.elem(idx));
   cnt += n(j);

   }

arma::mat XtOmegaXstar_full;
arma::vec x_trans_kappa_full;

if(p_g > 0){
   arma::mat G = X_star.cols(c, c*(1+p_g)-1);
   arma::mat G_trans = G.t();
   arma::mat GtOmegaXstar = G_trans * (X_star.each_col() % omega);
   XtOmegaXstar_full = arma::join_cols(XtOmegaXstar, GtOmegaXstar);
   x_trans_kappa_full = arma::join_cols(x_trans_kappa, G_trans*(kappa_adj));
   }
else{
   XtOmegaXstar_full = XtOmegaXstar;
   x_trans_kappa_full = x_trans_kappa;
   }

arma::mat cov_theta = inv_sympd(XtOmegaXstar_full +
                                (Sigma_inv_old*(rho_old*sum_R + 1.00 - rho_old)));

arma::vec mean_theta = cov_theta*(x_trans_kappa_full +
                                  rho_old*(Sigma_inv_old*((theta - theta_prior_mean_mat)*R)) + 
                                  (rho_old*sum_R + 1.00 - rho_old)*(Sigma_inv_old*theta_prior_mean_vec));

arma::mat ind_norms = arma::randn(1, (c*(1 + p_g)));
arma::vec theta_new = mean_theta + 
                      trans(ind_norms*arma::chol(cov_theta));

return(theta_new);

}