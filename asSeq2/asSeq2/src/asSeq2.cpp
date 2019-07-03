#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <stdio.h>
#include "shared.h"

// [[Rcpp::depends("RcppArmadillo")]]

// template<typename T>
// void printR_obj(const T& obj){
//   Rcpp::Rcout << obj << std::endl;
// }
// 
// const double LOG2PI = log(2*arma::datum::pi);
// 
// // --------------------
// // Intermediate Functions
// // --------------------
// 
// // [[Rcpp::export]]
// double Rcpp_norm(const arma::vec& a){
//   return arma::norm(a);
// }
// 
// // [[Rcpp::export]]
// double Rcpp_logSumExp(const arma::vec& log_x){
//   
//   if( log_x.n_elem == 1 ){
//     return log_x.at(0);
//   } else {
//     double max_val = max(log_x);
//     
//     arma::vec log_x_2 = log_x - max_val;
//     
//     return log(arma::sum(arma::exp(log_x_2))) + max_val;
//   }
// }
// 
// // [[Rcpp::export]]
// double Rcpp_min_diff(const arma::vec& x){
//   return arma::min(arma::diff(arma::sort(x)));
// }
// 
// // [[Rcpp::export]]
// arma::vec Rcpp_lgy_add_1(const arma::vec& y){
//   return arma::lgamma(y + 1);
// }

/* ---------------------------
 * TReC
 ---------------------------*/

using namespace Rcpp;

//[[Rcpp::export]]
void compute_offset(const double& bxj, const arma::vec& z,
                    arma::vec& offsets){
  
  arma::uword ii, pp;
  double tmp1 = std::log((1+std::exp(bxj))*0.5);
  
  for(ii =0; ii<z.n_elem; ii++){
    if(z.at(ii) == 2){
      offsets.at(ii) = bxj;
    }else if(z.at(ii) == 1){
      offsets.at(ii) = tmp1;
    }
  }
}

// [[Rcpp::export]]
double Rcpp_loglikNB_bxj(const double& phi,const arma::vec& mu1,
                         const arma::vec& y, const arma::vec& lgy1){
  //lgy1 = std::lgamma(y + 1)
  arma::uword ii;
  double vphi = 1.0/phi;
  double loglik = 0.0;
  
  for(ii = 0; ii< y.n_elem; ii++){
    if(y.at(ii) > 0){
      loglik +=  y.at(ii) * std::log(mu1.at(ii));
    }
    loglik += - (vphi+y.at(ii)) * std::log(vphi+mu1.at(ii));
  }
  
  return loglik;
}

// [[Rcpp::export]]
double Rcpp_loglik_pois(const arma::vec& mu1, const arma::vec& y,
                        const arma::vec& lgy1){
  //lgy1 = std::lgamma(y + 1)
  double loglik = 0.0;
  arma::uword ii;
  
  for(ii =0; ii<y.n_elem; ii++){
    if(y.at(ii) > 0){
      loglik +=  y.at(ii) * std::log(mu1.at(ii));
    }
    //loglik += y.at(ii)*std::log(mu1.at(ii)) - mu1.at(ii) - lgy1.at(ii);
    loglik += - mu1.at(ii) - lgy1.at(ii);
    //loglik += R::dpois(y.at(ii), mu1.at(ii), true);
  }
  return loglik;
}



// [[Rcpp::export]]
double Rcpp_logLTReC(const double& bxj, const arma::vec& y,
                     const arma::mat& X, const arma::vec& z,
                     const arma::vec& BETA, const double& phi,
                     const bool& fam_nb,
                     const arma::vec& lgy1, arma::vec& mu){
  // z is the genotype vector take value 0,1,2 (same as x in the R code)
  // lgy1 = lgamma(y+1)
  arma::uword ii, pp;
  double tmp1 = std::log((1+std::exp(bxj))*0.5);
  
  for(ii =0; ii<y.n_elem; ii++){
    if(z.at(ii) == 2){
      mu.at(ii) = std::exp( arma::dot(X.row(ii).t(), BETA) + bxj);
    }else if(z.at(ii) == 1){
      mu.at(ii) = std::exp( arma::dot(X.row(ii).t(), BETA) + tmp1);
    }else{
      mu.at(ii) = std::exp( arma::dot(X.row(ii).t(), BETA));
    }
  }
  
  if(fam_nb){
    return Rcpp_loglikNB(y, phi, lgy1, mu);
  }else{
    return Rcpp_loglik_pois(mu, y, lgy1);
  }
}

// [[Rcpp::export]]
double Rcpp_logLTReC_bxj(const double& bxj, const arma::vec& y,
                         const arma::mat& X, const arma::vec& z,
                         const arma::vec& BETA, const double& phi,
                         const bool& fam_nb,
                         const arma::vec& lgy1, arma::vec& mu){
  //ignored some coeffcient terms that is irrelavent to mu
  // z is the genotype vector take value 0,1,2 (same as x in the R code)
  // lgy1 = lgamma(y+1)
  arma::uword ii, pp;
  double tmp1 = std::log((1+std::exp(bxj))*0.5);
  
  for(ii =0; ii<y.n_elem; ii++){
    if(z.at(ii) == 2){
      mu.at(ii) = std::exp( arma::dot(X.row(ii).t(), BETA) + bxj);
    }else if(z.at(ii) == 1){
      mu.at(ii) = std::exp( arma::dot(X.row(ii).t(), BETA) + tmp1);
    }else{
      mu.at(ii) = std::exp( arma::dot(X.row(ii).t(), BETA));
    }
  }
  
  if(fam_nb){
    return Rcpp_loglikNB_bxj(phi, mu, y, lgy1);
  }else{
    return Rcpp_loglik_pois(mu, y, lgy1);
  }
}
// [[Rcpp::export]]
arma::vec Rcpp_grad_hess_bxj_trec(const double& bxj, const arma::vec& y,
                                  const arma::vec& z, const arma::vec& mu,
                                  const double& phi, const bool& fam_nb){
  arma::uword ii;
  arma::vec grad_hess = arma::zeros<arma::vec>(2);
  double df_dmu =0.0, dmu_db = 0.0, df_dmu2 =0.0;
  double tmp1 = std::exp(bxj)/(1.0+std::exp(bxj));
  
  for(ii =0; ii<y.n_elem; ii++){
    if(fam_nb){
      df_dmu  = y.at(ii)/mu.at(ii) - (1.0+phi*y.at(ii))/(1.0+phi*mu.at(ii));
      df_dmu2 = -y.at(ii)/std::pow(mu.at(ii), 2.0) +
        phi*(1.0+phi*y.at(ii))/std::pow(1.0+phi*mu.at(ii), 2.0);
    }else{
      df_dmu = y.at(ii)/mu.at(ii) - 1.0;
      df_dmu2 = -y.at(ii)/std::pow(mu.at(ii), 2.0);
    }
    if(z.at(ii) == 2){
      dmu_db = mu.at(ii);
    }else if(z.at(ii) == 1){
      dmu_db = mu.at(ii)*tmp1;
    }else{
      dmu_db = 0.0;
    }
    
    grad_hess.at(0) += df_dmu*dmu_db;
    grad_hess.at(1) += (df_dmu2*dmu_db*dmu_db + df_dmu*dmu_db);
  }
  return grad_hess;
}

// [[Rcpp::export]]
double Rcpp_trec_grad_bxj(const double& bxj, const arma::vec& y,
                          const arma::vec& z, const arma::vec& mu,
                          const double& phi, const bool& fam_nb){
  
  arma::uword ii;
  double grad = 0.0, df_dmu =0.0, dmu_db = 0.0, df_dmu2 =0.0;
  double tmp1 = std::exp(bxj)/(1.0+std::exp(bxj));
  
  for(ii =0; ii<y.n_elem; ii++){
    if(fam_nb){
      df_dmu  = y.at(ii)/mu.at(ii) - (1.0+phi*y.at(ii))/(1.0+phi*mu.at(ii));
    }else{
      df_dmu = y.at(ii)/mu.at(ii) - 1.0;
    }
    if(z.at(ii) == 2){
      dmu_db = mu.at(ii);
    }else if(z.at(ii) == 1){
      dmu_db = mu.at(ii)*tmp1;
    }else{
      dmu_db = 0.0;
    }
    
    grad += df_dmu*dmu_db;
  }
  return grad;
}

// [[Rcpp::export]]
Rcpp::List Rcpp_trec_bxj_BFGS(const double& bxj0, const arma::vec& y,
                              const arma::mat& X, const arma::vec& z,
                              const arma::vec& BETA, double& phi,
                              const bool& fam_nb, const arma::vec& lgy1,
                              const arma::uword& max_iter = 4e3,
                              const double& eps = 1e-7,const bool& show = true){
  
  arma::uword num_params = 1;
  arma::uword iter = 0;
  arma::uword jj,uu;
  arma::uword converge = 0;
  
  arma::vec xk = arma::zeros<arma::vec>(num_params);
  arma::mat inv_Bk = arma::eye<arma::mat>(num_params,num_params);
  arma::vec curr_xk = arma::zeros<arma::vec>(num_params);
  arma::mat I_num_params = arma::eye<arma::mat>(num_params,num_params);
  arma::vec new_xk = arma::zeros<arma::vec>(num_params);
  arma::vec gr_k = arma::zeros<arma::vec>(num_params);
  arma::vec p_k = arma::zeros<arma::vec>(num_params);
  arma::vec s_k = arma::zeros<arma::vec>(num_params);
  arma::vec y_k = arma::zeros<arma::vec>(num_params);
  arma::mat ISYT = arma::zeros<arma::mat>(num_params,num_params);
  arma::vec mu = arma::zeros<arma::vec>(y.n_elem);
  
  double old_LL,new_LL,inv_norm_p_k,tmp_alpha,ys;
  double fnscale = -1.0; // For maximization
  double curr_LL = 0.0;
  xk.at(0) = bxj0;
  
  while(iter < max_iter){
    //calculate direction p_k
    uu = 0;
    old_LL = fnscale * Rcpp_logLTReC_bxj(xk.at(0), y, X, z, BETA,
                                         phi, fam_nb, lgy1, mu);
    gr_k   = fnscale * Rcpp_trec_grad_bxj(xk.at(0), y, z, mu, phi, fam_nb);
    p_k    = -1.0 * inv_Bk * gr_k;
    inv_norm_p_k =  1.0 / std::max(1.0, Rcpp_norm(p_k));
    
    //line search for new xk
    for(jj=0; jj<15; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4, jj);
      new_xk    = xk + tmp_alpha * p_k;
      new_LL    = fnscale * Rcpp_logLTReC_bxj(new_xk.at(0), y, X, z, BETA,
                                              phi, fam_nb, lgy1, mu);
      
      if(new_LL < old_LL){ //minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * Rcpp_trec_grad_bxj(new_xk.at(0), y, z, mu,
                                           phi, fam_nb) - gr_k;
        ys  = arma::dot(y_k, s_k);
        
        if(ys > 0.0){
          // if(show) printR_obj("Update xk and inv_Bk");
          ISYT   = I_num_params - (s_k * y_k.t()) /ys;
          inv_Bk = ISYT * inv_Bk * ISYT.t() + s_k * s_k.t() / ys;
        }else{
          // if(show) printR_obj("Update xk only");
        }
        xk = new_xk;
        old_LL = new_LL;
        uu = 1;
        break;
      }
    }
    
    if(uu==0){
      if(Rcpp_norm(gr_k) > 1.0){
        // if(show) printR_obj("Reset inv_Bk");
        inv_Bk = I_num_params;
      }else{
        // if(show) printR_obj("Failed in search");
        break;
      }
    }
    
    //check convergence
    if(iter > 0){
      if(std::abs(curr_LL - old_LL) < eps &&
         Rcpp_norm(curr_xk - xk) < eps){
        gr_k = Rcpp_trec_grad_bxj(xk.at(0), y, z, mu, phi, fam_nb);
        if(Rcpp_norm(gr_k) < 0.01){
          converge = 1;
          break;
        }
      }
    }
    curr_xk = xk;
    curr_LL = old_LL;
    iter++;
  }
  
  old_LL = Rcpp_logLTReC(xk.at(0), y, X, z, BETA, phi, fam_nb, lgy1, mu);
  //gr_k = Rcpp_NB_reg_grad(y, X, mu, xk);
  return Rcpp::List::create(
    Rcpp::Named("converge", converge),
    Rcpp::Named("LL", old_LL),
    Rcpp::Named("iter", iter),
    Rcpp::Named("norm_GRAD", Rcpp_norm(gr_k)),
    Rcpp::Named("PAR", Rcpp::NumericVector(xk.begin(), xk.end()))
  );
}

// [[Rcpp::export]]
Rcpp::List Rcpp_trec_bxj(const arma::vec& y, const arma::mat& X,
                         const double& bxj, const arma::vec& z,
                         const arma::vec& BETA, double& phi,const bool& fam_nb,
                         const arma::vec& lgy1, const arma::uword& max_iter = 4e3,
                         const double& eps = 1e-7,const bool& show = true){
  // newton raphson's method update bxj 
  arma::uword iter = 0;
  arma::uword jj,uu;
  arma::uword converge = 0;
  arma::vec mu = arma::zeros<arma::vec>(y.n_elem);
  arma::vec g_h = arma::zeros<arma::vec>(2);
  double curr_LL = 0.0;
  double old_LL, new_LL, old_grad, old_hess_grad;
  double old_bxj = bxj, new_bxj = bxj, curr_bxj = bxj;
  
  while(iter < max_iter){
    old_LL = Rcpp_logLTReC(old_bxj, y, X, z, BETA, phi, fam_nb, lgy1, mu);
    g_h    = Rcpp_grad_hess_bxj_trec(old_bxj, y, z, mu, phi, fam_nb);
    old_grad = g_h.at(0);
    old_hess_grad = -1.0 * old_grad/g_h.at(1);
    uu = 0;
    
    for(jj =0; jj <= 15; jj++){
      new_bxj = old_bxj + old_hess_grad / std::pow(4.0,jj);
      new_LL  = Rcpp_logLTReC(new_bxj, y, X, z, BETA, phi, fam_nb, lgy1, mu);
      
      if(new_LL > old_LL){
        old_bxj = new_bxj;
        old_LL  = new_LL;
        uu = 1;
        break;
      }else{
        new_bxj = old_bxj + old_grad / std::pow(4.0,jj);
        new_LL  = Rcpp_logLTReC(new_bxj, y, X, z, BETA, phi, fam_nb, lgy1, mu);
        if(new_LL > old_LL){
          old_bxj = new_bxj;
          old_LL = new_LL;
          uu = 2;
          break;
        }
      }
    }
    // if(show){
    //   if(uu == 0){
    //     printR_obj("Failed update");
    //   } else if(uu == 1){
    //     printR_obj("Newton-Raphson update");
    //   } else {
    //     printR_obj("Gradient-Descent update");
    //   }
    // }
    if( uu == 0 ) break;
    
    if(iter>0){
      if( std::abs(curr_LL - old_LL) < eps && std::abs(curr_bxj - old_bxj) < eps ){
        g_h = Rcpp_grad_hess_bxj_trec(old_bxj, y, z, mu, phi, fam_nb);
        old_grad = g_h.at(0);
        old_hess_grad = -1.0 * old_grad/g_h.at(1);
        if( std::abs(old_grad) < eps && std::abs(old_hess_grad) < eps ){
          converge = 1;
          break;
        }
      }
    }
    curr_bxj = old_bxj;
    curr_LL = old_LL;
    iter++;
  }
  return Rcpp::List::create(
    Rcpp::Named("converge",converge),
    Rcpp::Named("LL",old_LL),
    Rcpp::Named("iter",iter),
    Rcpp::Named("norm_HG",std::abs(old_hess_grad)),
    Rcpp::Named("norm_G",std::abs(old_grad)),
    Rcpp::Named("PAR", old_bxj)
  );
}

/* ---------------------------
 * NB regression
 */
// [[Rcpp::export]]
double Rcpp_NB_reg_LL(const arma::vec& y, const arma::mat& X,
                      const arma::vec& offsets, const arma::vec& PARAMS,
                      const arma::vec& lgy1, arma::vec& mu){
  arma::uword ii;
  //double mu;
  double LL = 0.0;
  arma::uword pp = X.n_cols;
  arma::vec BETA = PARAMS.subvec(0,pp-1);
  double phi = std::exp(PARAMS.at(pp));
  double vphi = 1.0/phi;
  double lgvphi = lgamma(vphi);
  double vphi_lgvphi = vphi*std::log(vphi);
  
  for(ii = 0; ii< y.n_elem; ii++){
    mu.at(ii) = std::exp( arma::dot(X.row(ii).t(),BETA) + offsets.at(ii) );
    if(y.at(ii) > 0){
      LL += lgamma(y.at(ii) + vphi) - lgvphi - lgy1.at(ii) +
        y.at(ii) * std::log(mu.at(ii));
    }
    LL += vphi_lgvphi - (vphi+y.at(ii)) * std::log(vphi+mu.at(ii));
  }
  
  return LL;
  
}

// [[Rcpp::export]]
arma::vec Rcpp_NB_reg_grad(const arma::vec& y, const arma::mat& X,
                           const arma::vec& mu, const arma::vec& PARAMS){
  arma::uword ii;
  arma::uword pp = X.n_cols;
  arma::vec BETA = PARAMS.subvec(0,pp-1);
  double phi = std::exp(PARAMS.at(pp));
  double phi_mu1;
  double vphi = 1.0 / phi;
  arma::vec GRAD = arma::zeros<arma::vec>(pp + 1);
  double digam_vphi = R::digamma(vphi);
  double log_vphi = std::log(vphi);
  double vphi2 = -1.0* pow(vphi, 2.0);
  
  for(ii = 0; ii < y.n_elem; ii++){
    //mu = std::exp( arma::dot(X.row(ii).t(),BETA) + offsets.at(ii) );
    phi_mu1 = (1.0 + phi*mu.at(ii));
    // Part BETA
    GRAD.subvec(0,pp - 1) += (y.at(ii) - mu.at(ii))/phi_mu1 * X.row(ii).t();
    // Part PHI*
    GRAD.at(pp) += vphi2*(R::digamma(y.at(ii) + vphi) - digam_vphi -
      std::log(phi_mu1)) - (vphi+y.at(ii))*mu.at(ii)/phi_mu1 + y.at(ii)*vphi;
  }
  
  return GRAD;
}

// [[Rcpp::export]]
arma::mat Rcpp_NB_reg_Hess(const arma::vec& y, const arma::mat& X,
                           const arma::vec& mu, const arma::vec& PARAMS){
  arma::uword ii;
  arma::uword pp = X.n_cols;
  arma::vec BETA = PARAMS.subvec(0,pp-1);
  double phi = std::exp(PARAMS.at(pp));
  double phi_mu;
  double vphi = 1.0 / phi;
  arma::mat HESS = arma::zeros<arma::mat>(pp + 1,pp + 1);
  arma::vec hess_beta_phi = arma::zeros<arma::vec>(pp);
  double trigam_vphi = R::trigamma(vphi);
  double digam_vphi = R::digamma(vphi);
  double log_vphi = std::log(vphi);
  
  for(ii = 0; ii < y.n_elem; ii++){
    //mu.at(ii) = std::exp( arma::dot(X.row(ii).t(),BETA) + offsets.at(ii) );
    phi_mu =  phi*mu.at(ii);
    // Part2 BETA
    HESS.submat(0,0,pp-1,pp-1) += -1.0 * (vphi + y.at(ii)) /
      std::pow(vphi + mu.at(ii),2.0) *
        vphi * mu.at(ii) * X.row(ii).t() * X.row(ii);
    
    // Part2 BETA*phi*
    hess_beta_phi = -1.0 * (y.at(ii) - mu.at(ii)) / std::pow((1.0+phi_mu),2.0) *
      mu.at(ii) * X.row(ii).t();
    HESS(arma::span(0,pp - 1),pp) += hess_beta_phi;
    HESS(pp,arma::span(0,pp - 1)) += hess_beta_phi.t();
    
    // Part2 phi*
    HESS.at(pp,pp) += 2.0*pow(vphi, 3.0)*( R::digamma(y.at(ii)+vphi) - digam_vphi -
      log(1.0 + phi_mu) ) +
      pow(vphi, 4.0)*( R::trigamma(y.at(ii)+vphi) - trigam_vphi ) +
      2.0 * pow(vphi, 2.0)*mu.at(ii)/(1.0+phi_mu) - y.at(ii)*pow(vphi,2.0)
      + (vphi + y.at(ii))*pow(mu.at(ii), 2.0)/pow(1+phi_mu, 2.0);
  }
  
  return HESS;
}


/* ---------------------------
 * Poisson regression
 */

// [[Rcpp::export]]
double Rcpp_pois_reg_LL(const arma::vec& y, const arma::mat& X,
                        const arma::vec& offsets, const arma::vec& PARAMS,
                        const arma::vec& lgy1, arma::vec& mu){
  //lgy1 = std::lgamma(y + 1)
  double loglik = 0.0;
  arma::uword ii;
  arma::vec BETA = PARAMS.subvec(0, X.n_cols-1);
  
  for(ii =0; ii<y.n_elem; ii++){
    mu.at(ii) = std::exp( arma::dot(X.row(ii).t(),BETA) + offsets.at(ii) );
    loglik += y.at(ii)*std::log(mu.at(ii)) - mu.at(ii) - lgy1.at(ii);
    //loglik += R::dpois(y.at(ii), mu1.at(ii), true);
  }
  return loglik;
}

// [[Rcpp::export]]
arma::vec Rcpp_pois_reg_grad(const arma::vec& y, const arma::mat& X,
                             const arma::vec& mu, const arma::vec& PARAMS){
  arma::vec GRAD = arma::zeros<arma::vec>(X.n_cols) ;
  arma::uword ii;
  arma::vec BETA = PARAMS.subvec(0, X.n_cols-1);
  
  for(ii=0; ii<y.n_elem; ii++){
    GRAD += ((y.at(ii) - mu.at(ii)))* X.row(ii).t();
  }
  return GRAD;
}

// [[Rcpp::export]]
arma::mat Rcpp_pois_reg_Hess(const arma::vec& y, const arma::mat& X,
                             const arma::vec& mu, const arma::vec& PARAMS){
  arma::uword ii;
  arma::mat HESS = arma::zeros<arma::mat>(X.n_cols, X.n_cols);
  
  for(ii = 0; ii<y.n_elem; ii++){
    HESS += mu.at(ii) * X.row(ii).t() *X.row(ii);
  }
  
  return HESS;
}

/* ---------------------------
 * Trec Regression
 */
// [[Rcpp::export]]
double Rcpp_reg_LL(const arma::vec& y, const arma::mat& X,
                   const arma::vec& offsets, const arma::vec& PARAMS,
                   const bool& fam_nb,
                   const arma::vec& lgy1, arma::vec& mu){
  if(fam_nb){
    return Rcpp_NB_reg_LL(y, X, offsets, PARAMS, lgy1, mu);
  }else{
    return Rcpp_pois_reg_LL(y, X, offsets, PARAMS, lgy1, mu);
  }
}

// [[Rcpp::export]]
arma::vec Rcpp_reg_grad(const arma::vec& y, const arma::mat& X,
                        const arma::vec& mu, const arma::vec& PARAMS,
                        const bool& fam_nb){

  if(fam_nb){
    return Rcpp_NB_reg_grad(y, X, mu, PARAMS);
  }else{
    return Rcpp_pois_reg_grad(y, X, mu, PARAMS);
  }

}

// [[Rcpp::export]]
arma::mat Rcpp_reg_Hess(const arma::vec& y, const arma::mat& X,
                        const arma::vec& mu, const arma::vec& PARAMS,
                        const bool& fam_nb){
  
  if(fam_nb){
    return Rcpp_NB_reg_Hess(y, X, mu, PARAMS);
  }else{
    return Rcpp_pois_reg_Hess(y, X, mu, PARAMS);
  }
  
}


// [[Rcpp::export]]
Rcpp::List Rcpp_reg(const arma::vec& y, const arma::mat& X,
                    const arma::vec& offsets, const arma::vec& params0,
                    const bool& fam_nb, const arma::vec& lgy1,
                    const arma::uword& max_iter = 4e3,
                    const double& eps = 1e-7,const bool& show = true){
  // arma::uword N = X.n_rows;
  arma::uword iter = 0;
  arma::uword jj,uu;
  arma::uword converge = 0;
  arma::uword num_params = X.n_cols + 1;
  double curr_LL = 0.0;
  double old_LL,new_LL;
  arma::vec old_PARAMS = params0;
  arma::vec curr_PARAMS = params0;
  arma::vec new_PARAMS = params0;
  arma::vec old_grad = arma::zeros<arma::vec>(num_params);
  arma::vec mu = arma::zeros<arma::vec>(y.n_elem);
  arma::vec old_hess_grad = old_grad;
  
  while(iter < max_iter){
    old_LL = Rcpp_reg_LL(y,X,offsets,old_PARAMS,fam_nb,lgy1, mu);
    old_grad = Rcpp_reg_grad(y,X,mu,old_PARAMS,fam_nb);
    old_hess_grad = -1.0 * arma::inv(Rcpp_reg_Hess(y,X,mu,old_PARAMS,fam_nb))
      * old_grad;
    old_grad = old_grad / std::max(1.0,Rcpp_norm(old_grad));
    old_hess_grad = old_hess_grad / std::max(1.0,Rcpp_norm(old_hess_grad));
    uu = 0;
    for(jj = 0; jj <= 15; jj++){
      new_PARAMS = old_PARAMS + old_hess_grad / std::pow(4.0,jj);
      new_LL = Rcpp_reg_LL(y,X,offsets,new_PARAMS,fam_nb,lgy1,mu);
      if( new_LL > old_LL ){
        old_PARAMS = new_PARAMS;
        old_LL = new_LL;
        uu = 1;
        break;
      } else {
        new_PARAMS = old_PARAMS + old_grad / std::pow(4.0,jj);
        new_LL = Rcpp_reg_LL(y,X,offsets,new_PARAMS,fam_nb,lgy1,mu);
        if(new_LL > old_LL){
          old_PARAMS = new_PARAMS;
          old_LL = new_LL;
          uu = 2;
          break;
        }
      }
    }
    
    // if(show){
    //   if(uu == 0){
    //     printR_obj("Failed update");
    //   } else if(uu == 1){
    //     printR_obj("Newton-Raphson update");
    //   } else {
    //     printR_obj("Gradient-Descent update");
    //   }
    // }
    
    if( uu == 0 ) break;
    
    if(iter > 0){
      if( std::abs(curr_LL - old_LL) < eps &&
          Rcpp_norm(curr_PARAMS - old_PARAMS) < eps ){
        old_grad = Rcpp_reg_grad(y,X,mu,old_PARAMS,fam_nb);
        old_hess_grad = -1.0 * arma::inv(Rcpp_reg_Hess(y,X,mu,old_PARAMS,fam_nb))
          * old_grad;
        if( Rcpp_norm(old_grad) < eps && Rcpp_norm(old_hess_grad) < eps ){
          converge = 1;
          break;
        }
      }
    }
    
    curr_PARAMS = old_PARAMS;
    curr_LL = old_LL;
    iter++;
  }
  
  // return R_NilValue;
  return Rcpp::List::create(
    Rcpp::Named("converge",converge),
    Rcpp::Named("LL",old_LL),
    Rcpp::Named("iter",iter),
    Rcpp::Named("norm_HG",Rcpp_norm(old_hess_grad)),
    Rcpp::Named("norm_G",Rcpp_norm(old_grad)),
    Rcpp::Named("PARAMS",Rcpp::NumericVector(old_PARAMS.begin(),old_PARAMS.end()))
  );
}

// [[Rcpp::export]]
Rcpp::List Rcpp_reg_BFGS(const arma::vec& y, const arma::mat& X,
                         const arma::vec& offsets, const arma::vec& params0,
                         const bool& fam_nb, const arma::vec& lgy1,
                         const arma::uword& max_iter = 4e3,
                         const double& eps = 1e-7, const bool& show = true){

  arma::uword num_params = params0.n_elem;
  arma::uword iter = 0;
  arma::uword jj,uu;
  arma::uword converge = 0;

  arma::vec xk = params0;
  arma::mat inv_Bk = arma::eye<arma::mat>(num_params,num_params);
  arma::vec curr_xk = arma::zeros<arma::vec>(num_params);
  arma::mat I_num_params = arma::eye<arma::mat>(num_params,num_params);
  arma::vec new_xk = arma::zeros<arma::vec>(num_params);
  arma::vec gr_k = arma::zeros<arma::vec>(num_params);
  arma::vec p_k = arma::zeros<arma::vec>(num_params);
  arma::vec s_k = arma::zeros<arma::vec>(num_params);
  arma::vec y_k = arma::zeros<arma::vec>(num_params);
  arma::mat ISYT = arma::zeros<arma::mat>(num_params,num_params);
  arma::vec mu = arma::zeros<arma::vec>(y.n_elem);

  double old_LL,new_LL,inv_norm_p_k,tmp_alpha,ys;
  double fnscale = -1.0; // For maximization
  double curr_LL = 0.0;

  while(iter < max_iter){
    //calculate direction p_k
    uu = 0;
    old_LL = fnscale * Rcpp_reg_LL(y, X, offsets, xk,fam_nb, lgy1, mu);
    gr_k   = fnscale * Rcpp_reg_grad(y, X, mu, xk, fam_nb);
    p_k    = -1.0 * inv_Bk * gr_k;
    inv_norm_p_k =  1.0 / std::max(1.0, Rcpp_norm(p_k));

    //line search for new xk
    for(jj=0; jj<15; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4, jj);
      new_xk    = xk + tmp_alpha * p_k;
      new_LL    = fnscale * Rcpp_reg_LL(y, X, offsets, new_xk, fam_nb, lgy1, mu);

      if(new_LL < old_LL){ //minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * Rcpp_reg_grad(y, X, mu, new_xk,fam_nb) - gr_k;
        ys  = arma::dot(y_k, s_k);

        if(ys > 0.0){
          // if(show) printR_obj("Update xk and inv_Bk");
          ISYT   = I_num_params - (s_k * y_k.t()) /ys;
          inv_Bk = ISYT * inv_Bk * ISYT.t() + s_k * s_k.t() / ys;
        }else{
          // if(show) printR_obj("Update xk only");
        }
        xk = new_xk;
        old_LL = new_LL;
        uu = 1;
        break;
      }
    }

    if(uu==0){
      if(Rcpp_norm(gr_k) > 1.0){
        // if(show) printR_obj("Reset inv_Bk");
        inv_Bk = I_num_params;
      }else{
        // if(show) printR_obj("Failed in search");
        break;
      }
    }

    //check convergence
    if(iter > 0){
      if(std::abs(curr_LL - old_LL) < eps &&
         Rcpp_norm(curr_xk - xk) < eps){
        gr_k = Rcpp_reg_grad(y, X, mu, xk, fam_nb);
        if(Rcpp_norm(gr_k) < 0.01){
          converge = 1;
          break;
        }
      }
    }
    curr_xk = xk;
    curr_LL = old_LL;
    iter++;
  }

  old_LL = Rcpp_reg_LL(y, X, offsets, xk, fam_nb, lgy1, mu);
  //gr_k = Rcpp_NB_reg_grad(y, X, mu, xk);
  return Rcpp::List::create(
    Rcpp::Named("converge", converge),
    Rcpp::Named("LL", old_LL),
    Rcpp::Named("iter", iter),
    Rcpp::Named("norm_GRAD", Rcpp_norm(gr_k)),
    Rcpp::Named("PAR", Rcpp::NumericVector(xk.begin(), xk.end()))
  );
}


// can give initial value and null model regression result

// [[Rcpp::export]]
Rcpp::List Rcpp_trec(const arma::vec& y, const arma::mat& X,
                     const arma::vec& z, const bool& fam_nb,
                     const arma::vec& lgy1, const double& ini_bxj =0.0,
                     const double& LL_null = 0.0,
                     Rcpp::Nullable<Rcpp::NumericVector> ini_reg_par = R_NilValue,
                     const arma::uword& max_iter = 4e3,
                     const double& eps = 1e-7, const bool& show = false){
  arma::uword iter = 0;
  arma::uword converge = 0;
  arma::uword pp = X.n_cols+fam_nb;
  arma::vec offsets = arma::zeros<arma::vec>(z.n_elem);
  Rcpp::List new_reg, new_bxj_fit;
  double new_LL, new_bxj, phi, curr_LL, LL0;
  double curr_bxj = ini_bxj;
  arma::vec curr_reg_par = arma::zeros<arma::vec>(pp);
  arma::vec new_reg_par = arma::zeros<arma::vec>(pp);
  arma::vec BETA = arma::zeros<arma::vec>(pp-1);
  
  // initial regression fit
  
  if(ini_reg_par.isNull()){
    // if(show){
    //   printR_obj("begin initial regression fit");
    // }
    new_reg = Rcpp_reg_BFGS(y, X, offsets, curr_reg_par, fam_nb, lgy1,
                            max_iter, eps, false);
    curr_reg_par = as<arma::vec>(new_reg["PAR"]);
    LL0     = as<double>(new_reg["LL"]);
    curr_LL = LL0;
    
    // if(as<int>(new_reg["converge"]) != 1 ){
    //   printR_obj("initial regression fit dose not converge");
    // }
    
  }else{
    LL0 = LL_null;
    curr_LL = LL_null;
    curr_reg_par = Rcpp::as<arma::vec>(ini_reg_par);
  }
  BETA    = curr_reg_par.subvec(0, X.n_cols-1);
  phi     = std::exp(curr_reg_par.at(pp-1));
  //printR_obj(curr_reg_par);
  
  while(iter < max_iter){
    
    //update bxj
    new_bxj_fit = Rcpp_trec_bxj_BFGS(curr_bxj, y, X,  z, BETA, phi, fam_nb, lgy1,
                                     max_iter, eps, false);
    // newton's method
    // new_bxj_fit = Rcpp_trec_bxj(y, X, curr_bxj, z, BETA, phi, fam_nb, lgy1, 
    //                             max_iter, eps, false); // newton's method
    new_bxj     = as<double>(new_bxj_fit["PAR"]);
    new_LL      = as<double>(new_bxj_fit["LL"]);
    
    // if(as<int>(new_bxj_fit["converge"]) ==0){
    //   Rprintf("failed update bxj at %d iter \n bxj = %.2e\n PHI =%.2e \n",
    //           as<int>(new_bxj_fit["iter"]), new_bxj, phi);
    //   printR_obj(BETA);
    //   converge = 0;
    //   break;
    // }
    // if(show){
    //   Rprintf("TReC: bxj updated after %d iter \n",
    //           as<int>(new_bxj_fit["iter"]));
    //   //printR_obj(new_LL);
    //   //printR_obj(curr_LL);
    // }
    
    if(new_LL < curr_LL - eps){ 
      // log-like decrease and magnitude of it larger than eps
      //printR_obj("likelihood decreased for bxj");
      converge = 0;
      break;
    }
    
    //update BETA, phi
    compute_offset(new_bxj, z, offsets);
    new_reg     = Rcpp_reg_BFGS(y, X, offsets, curr_reg_par, fam_nb, lgy1,
                                max_iter, eps, false);
    
    new_LL      = as<double>(new_reg["LL"]);
    new_reg_par = as<arma::vec>(new_reg["PAR"]);
    
    BETA = new_reg_par.subvec(0, X.n_cols-1);
    phi  = std::exp(new_reg_par.at(pp-1));
    
    // if(show){
    //   Rprintf("TReC: BETA, PHI updated after %d iter \n", 
    //           as<int>(new_reg["iter"]));
    // }
    
    if(new_LL < curr_LL - eps){
      //printR_obj("likelihood decreased for betas");
      converge = 0;
      break;
    }
    //printR_obj(curr_LL);
    
    if(iter > 0){
      if(abs(curr_LL - new_LL) < eps && 
         Rcpp_norm(curr_reg_par - new_reg_par) < eps &&
         std::abs(curr_bxj - new_bxj) < eps){
        // convergence criteria 
        if((curr_LL - new_LL > 0.0 && as<double>(new_reg["norm_GRAD"]) > 0.01)){
          // magnitude of log-like decrease is samll but gradient is large
          // still consider not converged
          converge = 0;
        }else{
          converge = 1;
        }
        break;
        //
      }
      
      
    }
    curr_reg_par = new_reg_par;
    curr_bxj     = new_bxj;
    curr_LL      = new_LL;
    iter++;
    
  }
  if(show){
    Rprintf("TReC converges after %d iter \n", iter);
  }
  double lrt = (new_LL-LL0)*2.0;
  return Rcpp::List::create(
    Rcpp::Named("bxj", new_bxj),
    Rcpp::Named("reg_par", Rcpp::NumericVector(new_reg_par.begin(), new_reg_par.end())),
    Rcpp::Named("LL0", LL0),
    Rcpp::Named("LL", new_LL),
    Rcpp::Named("lrt", lrt),
    Rcpp::Named("pvalue", R::pchisq(lrt, 1, 0, 0)),
    Rcpp::Named("converge", converge),
    Rcpp::Named("iter", iter)
  );
}

/* ---------------------------
 * ASE
 ---------------------------*/

// [[Rcpp::export]]
double Rcpp_loglikBB(const arma::vec& ni, const arma::vec& ni0,
                     const double& bxj, const double& log_theta,
                     const arma::vec& lbc, const arma::vec& zeta){
  
  arma::uword ii;
  double loglik = 0.0;
  double vtheta = std::exp(-log_theta);
  double Pi1 = std::exp(bxj)/(1.0 + std::exp(bxj));
  double aa = Pi1*vtheta;
  double bb = vtheta - aa;
  double lgvt = lgamma(vtheta);
  double lgvab = - lgamma(aa) - lgamma(bb) + lgvt;
  double lgvabH0 = -2 *lgamma(0.5*vtheta) + lgvt;
  
  for(ii=0;ii<ni.n_elem;ii++){
    if(zeta.at(ii) == 1){   //het
      //printR_obj(loglik);
      loglik += lbc.at(ii) + lgamma(aa + ni0.at(ii)) +
        lgamma(bb + ni.at(ii) - ni0.at(ii)) + lgvab -
        lgamma(vtheta + ni.at(ii));
    }else{
      //printR_obj(loglik);
      loglik += lbc.at(ii) + lgamma(0.5*vtheta + ni0.at(ii)) +
        lgamma(0.5*vtheta + ni.at(ii) - ni0.at(ii)) + lgvabH0 -
        lgamma(vtheta + ni.at(ii)) ;
    }
  }
  
  return loglik;
}


// [[Rcpp::export]]
double Rcpp_loglikBB_Pi(const arma::vec& ni, const arma::vec& ni0,
                     const double& bxj, const double& log_theta,
                     const arma::vec& lbc, const arma::vec& zeta){
  // only used to update pi when theta is small 
  arma::uword ii;
  double loglik = 0.0;
  double vtheta = std::exp(-log_theta);
  double Pi1 = std::exp(bxj)/(1.0 + std::exp(bxj));
  double aa = Pi1*vtheta;
  double bb = vtheta - aa;
  //double lgvt = lgamma(vtheta);
  double lgvab = - lgamma(aa) - lgamma(bb); //+ lgvt
  //double lgvabH0 = -2 *lgamma(0.5*vtheta) + lgvt;
  
  for(ii=0;ii<ni.n_elem;ii++){
    if(zeta.at(ii) == 1){   //het
      //printR_obj(loglik);
      loglik += lbc.at(ii) + lgamma(aa + ni0.at(ii)) +
        lgamma(bb + ni.at(ii) - ni0.at(ii)) + lgvab;
    }
  }
  
  return loglik;
}

// [[Rcpp::export]]
arma::vec Rcpp_ase_grad(const arma::vec& ni, const arma::vec& ni0,
                        const double& bxj, const double& log_theta,
                        const arma::vec& zeta){
  // dASE_dlg_theta. 2019-04-26
  arma::uword ii;
  arma::vec grad = arma::zeros<arma::vec>(2);
  double vtheta = std::exp(-log_theta);
  double Pi1 = std::exp(bxj)/(1.0 + std::exp(bxj));
  double aa     = Pi1*vtheta;
  double bb     = vtheta - aa;
  double diaa   = R::digamma(aa);
  double dibb   = R::digamma(bb);
  double divt   = R::digamma(vtheta);
  double divt2  = R::digamma(0.5*vtheta);
  
  //Pi = 0.5
  double diaa_ni0, dibb_ni1;
  
  for(ii=0;ii<ni.n_elem;ii++){
    if(zeta.at(ii) == 1){
      diaa_ni0 = R::digamma(aa + ni0.at(ii));
      dibb_ni1 = R::digamma(bb + ni.at(ii) - ni0.at(ii));
      
      grad.at(0) += diaa_ni0 - dibb_ni1
        - diaa + dibb;
      
      //dlase_dtheta
      grad.at(1) += - Pi1*(diaa_ni0 - diaa) -
      (1.0 - Pi1) * (dibb_ni1 - dibb) -
      (divt - R::digamma(vtheta + ni.at(ii)));
      
    }else{
      
      //dlase_dpi
      //grad.at(0) += diaa_ni0_H0 - dibb_ni1_H0;
      
      //dlase_theta
      grad.at(1) += - 0.5 * R::digamma(0.5*vtheta + ni0.at(ii))
      - 0.5 * R::digamma(0.5*vtheta + ni.at(ii) - ni0.at(ii)) +
        divt2 - (divt - R::digamma(vtheta + ni.at(ii)));
    }
    //printR_obj(grad.at(1));
  }
  grad.at(0) *= vtheta * Pi1/(1.0 + std::exp(bxj)) ;
  // grad.at(1) *= pow(vtheta, 2.0);
  grad.at(1) *= vtheta;
  return grad;
}

// [[Rcpp::export]]
double Rcpp_ase_grad_Pi(const arma::vec& ni, const arma::vec& ni0,
                        const double& bxj, const double& log_theta,
                        const arma::vec& zeta){
  
  arma::uword ii;
  double grad = 0.0;
  double vtheta = std::exp(-log_theta);
  double Pi1 = std::exp(bxj)/(1.0 + std::exp(bxj));
  double aa     = Pi1*vtheta;
  double bb     = vtheta - aa;
  double diaa   = R::digamma(aa);
  double dibb   = R::digamma(bb);
  double divt   = R::digamma(vtheta);
  
  //Pi = 0.5
  double diaa_ni0, dibb_ni1;
  
  for(ii=0;ii<ni.n_elem;ii++){
    if(zeta.at(ii) == 1){
      diaa_ni0 = R::digamma(aa + ni0.at(ii));
      dibb_ni1 = R::digamma(bb + ni.at(ii) - ni0.at(ii));
      
      grad += diaa_ni0 - dibb_ni1 - diaa + dibb;
      
    }
    // else{
    //
    //   grad += R::digamma(0.5*vtheta + ni0.at(ii)) -
    //     R::digamma(0.5*vtheta + ni.at(ii) - ni0.at(ii));
    //
    // }
    //printR_obj(grad.at(1));
  }
  grad *= vtheta * Pi1/(1.0 + std::exp(bxj)) ;
  
  return grad;
}

// [[Rcpp::export]]
double Rcpp_ase_grad_H0(const arma::vec& ni, const arma::vec& ni0,
                        const double& bxj, const double& log_theta,
                        const arma::vec& zeta){
  //under H0 or only update theta
  // dASE_dlg_theta. 2019-04-26
  arma::uword ii;
  double grad   = 0.0;
  double vtheta = std::exp(-log_theta);
  double Pi1 = std::exp(bxj)/(1.0 + std::exp(bxj));
  double aa     = Pi1*vtheta;
  double bb     = vtheta - aa;
  double diaa   = R::digamma(aa);
  double dibb   = R::digamma(bb);
  double divt   = R::digamma(vtheta);
  double divt2  = R::digamma(0.5*vtheta);
  //Pi = 0.5
  double diaa_ni0, dibb_ni1;
  
  for(ii=0;ii<ni.n_elem;ii++){
    if(zeta.at(ii) == 1){
      diaa_ni0 = R::digamma(aa + ni0.at(ii));
      dibb_ni1 = R::digamma(bb + ni.at(ii) - ni0.at(ii));
      
      //dlase_dtheta
      grad += - Pi1*(diaa_ni0 - diaa) -
        (1.0 - Pi1) * (dibb_ni1 - dibb) -
        (divt - R::digamma(vtheta + ni.at(ii)));
      
    }else{
      
      //dlase_theta
      grad += - 0.5 * R::digamma(0.5*vtheta + ni0.at(ii))
      - 0.5 * R::digamma(0.5*vtheta + ni.at(ii) - ni0.at(ii)) +
        divt2 - (divt - R::digamma(vtheta + ni.at(ii)));
    }
    
  }
  // grad *= pow(vtheta, 2.0);
  grad *= vtheta;
  return grad;
}


// // [[Rcpp::export]]
// arma::mat Rcpp_ase_hess(const arma::vec& ni, const arma::vec& ni0,
//                         const double& Pi1, const double& theta,
//                         const arma::vec& zeta){
//   arma::uword ii;
//   arma::mat hess = arma::zeros<arma::mat>(2,2);
//   double vtheta = 1/theta;
//   double aa     = Pi1*vtheta;
//   double bb     = vtheta - aa;
//   double traa   = R::trigamma(aa);
//   double trbb   = R::trigamma(bb);
//   double diaa   = R::digamma(aa);
//   double dibb   = R::digamma(bb);
//   double trvt   = R::trigamma(vtheta);
//   double divt   = R::digamma(vtheta);
//
//   double w = 0.0, dw_dtheta = 0.0;
//
//   double traa_ni0, trbb_ni1, diaa_ni0, dibb_ni1;
//
//   for(ii=0;ii<ni.n_elem;ii++){
//     if(zeta.at(ii) == 1){
//       traa_ni0 = R::trigamma(aa + ni0.at(ii));
//       trbb_ni1 = R::trigamma(bb + ni.at(ii) - ni0.at(ii));
//       diaa_ni0 = R::digamma(aa + ni0.at(ii));
//       dibb_ni1 = R::digamma(bb + ni.at(ii) - ni0.at(ii));
//
//       //pi
//       hess.at(0,0) += traa_ni0 + trbb_ni1
//         - traa - trbb;
//
//       //theta
//       w +=  Pi1*(diaa_ni0 - diaa) +
//       (1.0 - Pi1) * (dibb_ni1 - dibb) +
//       (divt - R::digamma(vtheta + ni.at(ii)));
//
//       dw_dtheta += pow(Pi1, 2.0)*(traa_ni0 - traa) +
//         pow((1 - Pi1), 2.0) * (trbb_ni1 - trbb) +
//         (trvt - R::trigamma(vtheta + ni.at(ii)));
//
//        // pi*theta
//        hess.at(0,1) += - pow(vtheta, 2.0) * ((diaa_ni0 - diaa) -
//          (dibb_ni1 - dibb)) -
//          pow(vtheta, 3.0)*Pi1*(traa_ni0 - traa) +
//          pow(vtheta, 3.0)*(1.0-Pi1)*(trbb_ni1 - trbb) ;
//
//     }else{
//
//       //Pi
//       //hess.at(0,0) += diaa_ni0_H0 - dibb_ni1_H0;
//
//       //theta
//       w +=  0.5 * R::digamma(0.5*vtheta + ni0.at(ii))
//       + 0.5 * R::digamma(0.5*vtheta + ni.at(ii) - ni0.at(ii)) -
//         R::digamma(0.5*vtheta) +
//         (divt - R::digamma(vtheta + ni.at(ii)));
//
//       dw_dtheta += pow(0.5, 2.0)*(R::trigamma(0.5*vtheta + ni0.at(ii))) -
//         0.5*R::trigamma(0.5*vtheta) +
//         pow(0.5, 2.0) * (R::trigamma(0.5*vtheta + ni.at(ii) - ni0.at(ii))) +
//         (trvt - R::trigamma(vtheta + ni.at(ii)));;
//
//       // pi*theta
//       hess.at(0,1) += - pow(vtheta, 2.0) *
//         0.5*(R::digamma(0.5*vtheta + ni0.at(ii)) -
//         R::digamma(0.5*vtheta + ni.at(ii) - ni0.at(ii))) -
//         pow(vtheta, 3.0)*0.5*
//         (R::trigamma(0.5*vtheta + ni0.at(ii)) +
//         R::trigamma(0.5*vtheta + ni.at(ii) - ni0.at(ii)) -
//         2*R::trigamma(0.5*vtheta));
//
//       }
//
//   }
//   hess.at(0,0) *= pow(vtheta, 2.0);
//   //printR_obj(w);
//   //printR_obj(dw_dtheta);
//   hess.at(1,1) = -2.0*pow(vtheta, 3.0)*w + pow(vtheta, 4.0)*dw_dtheta;
//   hess.at(1,0) = hess.at(0,1);
//
//   return hess;
// } // diagonal is not correct
// [[Rcpp::export]]
Rcpp::List Rcpp_ase_Pi_BFGS(const arma::vec& ni, const arma::vec& ni0,
                            const arma::vec& zeta, const double& bxj,
                            const double& lg_theta, const arma::vec& lbc,
                            const arma::uword& max_iter = 4e3,
                            const double& eps = 1e-7, const bool& show = true){
  // used to update Pi when theta is small
  arma::uword num_params = 1;
  arma::uword iter = 0;
  arma::uword jj,uu;
  arma::uword converge = 0;
  
  arma::vec xk = arma::zeros<arma::vec>(num_params);
  arma::mat inv_Bk = arma::eye<arma::mat>(num_params,num_params);
  arma::vec curr_xk = arma::zeros<arma::vec>(num_params);
  arma::mat I_num_params = arma::eye<arma::mat>(num_params,num_params);
  arma::vec new_xk = arma::zeros<arma::vec>(num_params);
  arma::vec gr_k = arma::zeros<arma::vec>(num_params);
  arma::vec p_k = arma::zeros<arma::vec>(num_params);
  arma::vec s_k = arma::zeros<arma::vec>(num_params);
  arma::vec y_k = arma::zeros<arma::vec>(num_params);
  arma::mat ISYT = arma::zeros<arma::mat>(num_params,num_params);
  
  double old_LL,new_LL,inv_norm_p_k,tmp_alpha,ys;
  double fnscale = -1.0; // For maximization
  double curr_LL = 0.0;
  
  xk.at(0) = bxj;
  
  while(iter < max_iter){
    //calculate direction p_k
    uu = 0;
    
    old_LL = fnscale * Rcpp_loglikBB_Pi(ni, ni0,xk.at(0), lg_theta, lbc, zeta);
    
    //if(old_LL < 0) break;
    
    //printR_obj(old_LL);
    gr_k   = fnscale * Rcpp_ase_grad_Pi(ni, ni0, xk.at(0), lg_theta, zeta);
    p_k    = -1.0 * inv_Bk * gr_k;
    inv_norm_p_k =  1.0 / std::max(1.0, Rcpp_norm(p_k));
    
    //line search for new xk
    for(jj=0; jj<15; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4, jj);
      new_xk    = xk + tmp_alpha * p_k;
      new_LL    = fnscale * Rcpp_loglikBB_Pi(ni, ni0, new_xk.at(0),
                                             lg_theta, lbc, zeta);
      
      if(new_LL < old_LL){ //minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * Rcpp_ase_grad_Pi(ni, ni0, new_xk.at(0), 
                                         lg_theta, zeta) - gr_k;
        ys  = arma::dot(y_k, s_k);
        
        if(ys > 0.0){
          // if(show) printR_obj("Update xk and inv_Bk");
          ISYT   = I_num_params - (s_k * y_k.t()) /ys;
          inv_Bk = ISYT * inv_Bk * ISYT.t() + s_k * s_k.t() / ys;
        }else{
          // if(show) printR_obj("Update xk only");
        }
        xk = new_xk;
        old_LL = new_LL;
        uu = 1;
        break;
      }
    }
    //printR_obj(new_LL);
    //printR_obj(xk);
    
    if(uu==0){
      if(Rcpp_norm(gr_k) > 1.0){
        // if(show) printR_obj("Reset inv_Bk");
        inv_Bk = I_num_params;
      }else{
        // if(show) printR_obj("Failed in search");
        break;
      }
    }
    
    //check convergence
    if(iter > 0){
      if(std::abs(curr_LL - old_LL) < eps &&
         Rcpp_norm(curr_xk - xk) < eps){
        gr_k = Rcpp_ase_grad_Pi(ni, ni0, xk.at(0), lg_theta, zeta);
        if(Rcpp_norm(gr_k) < 0.01){
          converge = 1;
          break;
        }
      }
    }
    curr_xk = xk;
    curr_LL = old_LL;
    iter++;
  }
  
  old_LL = Rcpp_loglikBB(ni, ni0, xk.at(0),lg_theta, lbc, zeta);
  
  return Rcpp::List::create(
    Rcpp::Named("converge", converge),
    Rcpp::Named("LL", old_LL),
    Rcpp::Named("iter", iter),
    Rcpp::Named("norm_GRAD", Rcpp_norm(gr_k)),
    Rcpp::Named("PAR", Rcpp::NumericVector(xk.begin(), xk.end()))
  );
}

// [[Rcpp::export]]
Rcpp::List Rcpp_ase_BFGS(const arma::vec& ni, const arma::vec& ni0,
                         const arma::vec& zeta,
                         const arma::vec& params0, const arma::vec& lbc,
                         const arma::uword& max_iter = 4e3,
                         const double& eps = 1e-7, const bool& show = true){
  //optimaize theta and pi
  arma::uword num_params = params0.n_elem;
  arma::uword iter = 0;
  arma::uword jj,uu;
  arma::uword converge = 0;
  
  arma::vec xk = params0;
  arma::mat inv_Bk = arma::eye<arma::mat>(num_params,num_params);
  arma::vec curr_xk = arma::zeros<arma::vec>(num_params);
  arma::mat I_num_params = arma::eye<arma::mat>(num_params,num_params);
  arma::vec new_xk = arma::zeros<arma::vec>(num_params);
  arma::vec gr_k = arma::zeros<arma::vec>(num_params);
  arma::vec p_k = arma::zeros<arma::vec>(num_params);
  arma::vec s_k = arma::zeros<arma::vec>(num_params);
  arma::vec y_k = arma::zeros<arma::vec>(num_params);
  arma::mat ISYT = arma::zeros<arma::mat>(num_params,num_params);
  
  double old_LL,new_LL,inv_norm_p_k,tmp_alpha,ys;
  double fnscale = -1.0; // For maximization
  double curr_LL = 0.0;
  
  while(iter < max_iter){
    
    //calculate direction p_k
    uu = 0;
    
    old_LL = fnscale * Rcpp_loglikBB(ni, ni0, xk.at(0), xk.at(1), lbc, zeta);
    gr_k   = fnscale * Rcpp_ase_grad(ni, ni0, xk.at(0), xk.at(1), zeta);
    
    //if(old_LL < 0 ) break;
    
    p_k    = -1.0 * inv_Bk * gr_k;
    inv_norm_p_k =  1.0 / std::max(1.0, Rcpp_norm(p_k));
    
    //line search for new xk
    for(jj=0; jj<15; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4, jj);
      new_xk    = xk + tmp_alpha * p_k;
      
      new_LL    = fnscale * Rcpp_loglikBB(ni, ni0, new_xk.at(0),
                                          new_xk.at(1), lbc, zeta);
      if(new_LL < old_LL){ //minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * Rcpp_ase_grad(ni, ni0, new_xk.at(0),
                                      new_xk.at(1), zeta) - gr_k;
        ys  = arma::dot(y_k, s_k);
        
        if(ys > 0.0){
          // if(show) printR_obj("Update xk and inv_Bk");
          ISYT   = I_num_params - (s_k * y_k.t()) /ys;
          inv_Bk = ISYT * inv_Bk * ISYT.t() + s_k * s_k.t() / ys;
        }else{
          // if(show) printR_obj("Update xk only");
        }
        xk = new_xk;
        old_LL = new_LL;
        uu = 1;
        break;
      }
    }
    //printR_obj(new_LL);
    //printR_obj(xk);
    
    if(uu==0){
      if(Rcpp_norm(gr_k) > 1.0){
        // if(show) printR_obj("Reset inv_Bk");
        inv_Bk = I_num_params;
      }else{
        // if(show) printR_obj("Failed in search");
        break;
      }
    }
    
    // if theta < 1e-5 stop updating it
    if(xk.at(1) < -11.5){
      Rcpp::List opH1 = Rcpp_ase_Pi_BFGS(ni, ni0, zeta, xk.at(0), xk.at(1),
                                         lbc, max_iter, eps, show);
      converge = opH1["converge"];
      xk.at(0) = as<double>(opH1["PAR"]);
      iter += as<int>(opH1["iter"]);
      gr_k = as<double>(opH1["norm_GRAD"]);
      break;
    }
    
    //check convergence
    if(iter > 0){
      if(std::abs(curr_LL - old_LL) < eps &&
         Rcpp_norm(curr_xk - xk) < eps){
        gr_k = Rcpp_ase_grad(ni, ni0, xk.at(0), xk.at(1), zeta);
        if(Rcpp_norm(gr_k) < 0.01){
          converge = 1;
          break;
        }
      }
    }
    curr_xk = xk;
    curr_LL = old_LL;
    iter++;
    
  }
  
  old_LL = Rcpp_loglikBB(ni, ni0, xk.at(0), xk.at(1), lbc, zeta);
  return Rcpp::List::create(
    Rcpp::Named("converge", converge),
    Rcpp::Named("LL", old_LL),
    Rcpp::Named("iter", iter),
    Rcpp::Named("norm_GRAD", Rcpp_norm(gr_k)),
    Rcpp::Named("PAR", Rcpp::NumericVector(xk.begin(), xk.end()))
  );
}

// [[Rcpp::export]]
Rcpp::List Rcpp_ase_theta_BFGS(const arma::vec& ni, const arma::vec& ni0,
                               const arma::vec& zeta, const double& bxj,
                               const double& lg_theta, const arma::vec& lbc,
                               const arma::uword& max_iter = 4e3,
                               const double& eps = 1e-7, const bool& show = true){
  
  arma::uword num_params = 1;
  arma::uword iter = 0;
  arma::uword jj,uu;
  arma::uword converge = 0;
  
  arma::vec xk = arma::zeros<arma::vec>(num_params);
  arma::mat inv_Bk = arma::eye<arma::mat>(num_params,num_params);
  arma::vec curr_xk = arma::zeros<arma::vec>(num_params);
  arma::mat I_num_params = arma::eye<arma::mat>(num_params,num_params);
  arma::vec new_xk = arma::zeros<arma::vec>(num_params);
  arma::vec gr_k = arma::zeros<arma::vec>(num_params);
  arma::vec p_k = arma::zeros<arma::vec>(num_params);
  arma::vec s_k = arma::zeros<arma::vec>(num_params);
  arma::vec y_k = arma::zeros<arma::vec>(num_params);
  arma::mat ISYT = arma::zeros<arma::mat>(num_params,num_params);
  
  double old_LL,new_LL,inv_norm_p_k,tmp_alpha,ys;
  double fnscale = -1.0; // For maximization
  double curr_LL = 0.0;
  
  xk.at(0) = lg_theta;
  
  while(iter < max_iter){
    //calculate direction p_k
    uu = 0;
    
    old_LL = fnscale * Rcpp_loglikBB(ni, ni0, bxj, xk.at(0), lbc, zeta);
    
    //if(old_LL < 0) break;
    
    //printR_obj(old_LL);
    gr_k   = fnscale * Rcpp_ase_grad_H0(ni, ni0, bxj, xk.at(0), zeta);
    p_k    = -1.0 * inv_Bk * gr_k;
    inv_norm_p_k =  1.0 / std::max(1.0, Rcpp_norm(p_k));
    
    //line search for new xk
    for(jj=0; jj<15; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4, jj);
      new_xk    = xk + tmp_alpha * p_k;
      new_LL    = fnscale * Rcpp_loglikBB(ni, ni0,bxj,
                                          new_xk.at(0), lbc, zeta);
      
      if(new_LL < old_LL){ //minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * Rcpp_ase_grad_H0(ni, ni0, bxj,
                                         new_xk.at(0), zeta) - gr_k;
        ys  = arma::dot(y_k, s_k);
        
        if(ys > 0.0){
          // if(show) printR_obj("Update xk and inv_Bk");
          ISYT   = I_num_params - (s_k * y_k.t()) /ys;
          inv_Bk = ISYT * inv_Bk * ISYT.t() + s_k * s_k.t() / ys;
        }else{
          // if(show) printR_obj("Update xk only");
        }
        xk = new_xk;
        old_LL = new_LL;
        uu = 1;
        break;
      }
    }
    //printR_obj(new_LL);
    //printR_obj(xk);
    
    if(uu==0){
      if(Rcpp_norm(gr_k) > 1.0){
        // if(show) printR_obj("Reset inv_Bk");
        inv_Bk = I_num_params;
      }else{
        // if(show) printR_obj("Failed in search");
        break;
      }
    }
    
    if(xk.at(0) < -11.5){
      converge = 1;
      break;
    }
    
    //check convergence
    if(iter > 0){
      if(std::abs(curr_LL - old_LL) < eps &&
         Rcpp_norm(curr_xk - xk) < eps){
        gr_k = Rcpp_ase_grad_H0(ni, ni0, bxj, xk.at(0), zeta);
        if(Rcpp_norm(gr_k) < 0.01){
          converge = 1;
          break;
        }
      }
    }
    curr_xk = xk;
    curr_LL = old_LL;
    iter++;
  }
  
  old_LL = Rcpp_loglikBB(ni, ni0, bxj, xk.at(0), lbc, zeta);
  
  return Rcpp::List::create(
    Rcpp::Named("converge", converge),
    Rcpp::Named("LL", old_LL),
    Rcpp::Named("iter", iter),
    Rcpp::Named("norm_GRAD", Rcpp_norm(gr_k)),
    Rcpp::Named("PAR", Rcpp::NumericVector(xk.begin(), xk.end()))
  );
}

// [[Rcpp::export]]
Rcpp::List Rcpp_ase(const arma::vec& ni, const arma::vec& ni0,
                    const arma::vec& zeta, const arma::vec& lbc,
                    const double& theta0 = 0.0,
                    const double& LL_null_ase = 0.0,
                    const arma::uword& max_iter = 4e3,
                    const double& eps = 1e-7, const bool& show = true){
  
  Rcpp::List opH0, opH1;
  arma::vec par = arma::zeros<arma::vec>(2);
  double LL0, lg_theta;
  
  if(theta0 < 0.0){
    opH0 = Rcpp_ase_theta_BFGS(ni, ni0, zeta, 0.0, -2.0, lbc, 
                               max_iter, eps, show);
    lg_theta = as<double>(opH0["PAR"]);
    LL0 = as<double>(opH0["LL"]);
  }else{
    lg_theta = std::log(theta0);
    LL0 = LL_null_ase;
  }
  
    par.at(0) = 0.0;
    par.at(1) = lg_theta; 
    opH1 = Rcpp_ase_BFGS(ni, ni0, zeta, par, lbc, max_iter, eps, show);


  
  //printR_obj(par);
  
  
  double lrt = (as<double>(opH1["LL"]) - LL0)*2.0;
  return Rcpp::List::create(
    Rcpp::Named("par0", lg_theta),
    Rcpp::Named("par", opH1["PAR"]),
    Rcpp::Named("LL0", LL0),
    Rcpp::Named("LL", opH1["LL"]),
    Rcpp::Named("lrt", lrt),
    Rcpp::Named("pvalue", R::pchisq(lrt, 1, 0, 0)),
    Rcpp::Named("converge",  opH1["converge"] )
  );
}

/* ---------------------------
 * TRECASE
 ---------------------------*/

// [[Rcpp::export]]
double Rcpp_trecase_LL(const double& bxj, const arma::vec& y,
                       const arma::mat& X, const arma::vec& z,
                       const arma::vec& BETA, const double& phi,
                       const bool& fam_nb,const arma::vec& lgy1,
                       arma::vec& mu,
                       const arma::vec& ni, const arma::vec& ni0,
                       const double& log_theta, const arma::vec& lbc,
                       const arma::vec& zeta ){
  // double Pi1 = std::exp(bxj)/(1.0 + exp(bxj));
  
  return(Rcpp_logLTReC(bxj, y, X, z, BETA, phi, fam_nb, lgy1, mu)
           + Rcpp_loglikBB(ni, ni0, bxj, log_theta, lbc, zeta));
  
}

// [[Rcpp::export]]
double Rcpp_trecase_grad_bxj(const double& bxj, const arma::vec& y,
                             const arma::mat& X, const arma::vec& z,
                             const arma::vec& BETA, const double& phi,
                             const bool& fam_nb,const arma::vec& lgy1,
                             const arma::vec& mu,
                             const arma::vec& ni, const arma::vec& ni0,
                             const double& log_theta,
                             const arma::vec& lbc, const arma::vec& zeta){
  // double Pi1 = std::exp(bxj)/(1.0 + exp(bxj));
  arma::vec trec_grad = arma::zeros<arma::vec>(2);
  //trec_grad = Rcpp_grad_hess_bxj_trec(bxj, y, z, mu, phi, fam_nb);
  return(Rcpp_trec_grad_bxj(bxj, y, z, mu, phi, fam_nb)+
         Rcpp_ase_grad_Pi(ni, ni0, bxj, log_theta, zeta));
  
}

// [[Rcpp::export]]
Rcpp::List Rcpp_trecase_BFGS(const double& bxj0, const arma::vec& y,
                             const arma::mat& X, const arma::vec& z,
                             const arma::vec& BETA, const double& phi,
                             const bool& fam_nb, const arma::vec& lgy1,
                             const arma::vec& ni, const arma::vec& ni0,
                             const double& log_theta,
                             const arma::vec& lbc, const arma::vec& zeta,
                             const arma::uword& max_iter = 4e3,
                             const double& eps = 1e-7, const bool& show = true){
  
  arma::uword num_params = 1;
  arma::uword iter = 0;
  arma::uword jj,uu;
  arma::uword converge = 0;
  
  arma::vec xk = arma::zeros<arma::vec>(num_params);
  xk.at(0) = bxj0;
  arma::mat inv_Bk = arma::eye<arma::mat>(num_params,num_params);
  arma::vec curr_xk = arma::zeros<arma::vec>(num_params);
  arma::mat I_num_params = arma::eye<arma::mat>(num_params,num_params);
  arma::vec new_xk = arma::zeros<arma::vec>(num_params);
  arma::vec gr_k = arma::zeros<arma::vec>(num_params);
  arma::vec p_k = arma::zeros<arma::vec>(num_params);
  arma::vec s_k = arma::zeros<arma::vec>(num_params);
  arma::vec y_k = arma::zeros<arma::vec>(num_params);
  arma::mat ISYT = arma::zeros<arma::mat>(num_params,num_params);
  arma::vec mu = arma::zeros<arma::vec>(y.n_elem);
  
  double old_LL,new_LL,inv_norm_p_k,tmp_alpha,ys;
  double fnscale = -1.0; // For maximization
  double curr_LL = 0.0;
  
  while(iter < max_iter){
    //calculate direction p_k
    uu = 0;
    old_LL = fnscale * Rcpp_trecase_LL(xk.at(0), y, X, z, BETA, phi, fam_nb, lgy1,
                                       mu, ni, ni0, log_theta, lbc, zeta);
    gr_k   = fnscale * Rcpp_trecase_grad_bxj(xk.at(0), y, X, z, BETA, phi, fam_nb,
                                             lgy1, mu, ni, ni0,
                                             log_theta, lbc, zeta);
    p_k    = -1.0 * inv_Bk * gr_k;
    inv_norm_p_k =  1.0 / std::max(1.0, Rcpp_norm(p_k));
    
    //line search for new xk
    for(jj=0; jj<15; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4, jj);
      new_xk    = xk + tmp_alpha * p_k;
      new_LL    = fnscale * Rcpp_trecase_LL(new_xk.at(0), y, X, z, BETA, phi, fam_nb,
                                            lgy1, mu, ni, ni0, log_theta, lbc, zeta);
      if(new_LL < old_LL){ //minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * Rcpp_trecase_grad_bxj(new_xk.at(0), y, X, z, BETA, phi, fam_nb,
                                              lgy1, mu, ni, ni0,
                                              log_theta, lbc, zeta) - gr_k;
        ys  = arma::dot(y_k, s_k);
        
        if(ys > 0.0){
          // if(show) printR_obj("Update xk and inv_Bk");
          ISYT   = I_num_params - (s_k * y_k.t()) /ys;
          inv_Bk = ISYT * inv_Bk * ISYT.t() + s_k * s_k.t() / ys;
        }else{
          // if(show) printR_obj("Update xk only");
        }
        xk = new_xk;
        old_LL = new_LL;
        uu = 1;
        break;
      }
    }
    
    if(uu==0){
      if(Rcpp_norm(gr_k) > 1.0){
        // if(show) printR_obj("Reset inv_Bk");
        inv_Bk = I_num_params;
      }else{
        // if(show) printR_obj("Failed in search");
        break;
      }
    }
    
    //check convergence
    if(iter > 0){
      if(std::abs(curr_LL - old_LL) < eps &&
         Rcpp_norm(curr_xk - xk) < eps){
        gr_k = Rcpp_trecase_grad_bxj(xk.at(0), y, X, z, BETA, phi, fam_nb,
                                     lgy1, mu, ni, ni0,
                                     log_theta, lbc, zeta);
        if(Rcpp_norm(gr_k) < 0.01){
          converge = 1;
          break;
        }
      }
    }
    curr_xk = xk;
    curr_LL = old_LL;
    iter++;
  }
  
  old_LL = Rcpp_trecase_LL(new_xk.at(0), y, X, z, BETA, phi, fam_nb,
                           lgy1, mu, ni, ni0, log_theta, lbc, zeta);
  
  return Rcpp::List::create(
    Rcpp::Named("converge", converge),
    Rcpp::Named("LL", old_LL),
    Rcpp::Named("iter", iter),
    Rcpp::Named("norm_GRAD", Rcpp_norm(gr_k)),
    Rcpp::Named("PAR", xk.at(0))
  );
}

// [[Rcpp::export]]
Rcpp::List Rcpp_trecase(const arma::vec& y, const arma::mat& X,
                        const arma::vec& z, const bool& fam_nb,
                        const arma::vec& lgy1,
                        const arma::vec& ni, const arma::vec& ni0,
                        const arma::vec& zeta, const arma::vec& lbc,
                        const double& ini_bxj =0.0,
                        const double& LL_null = 0.0,
                        Rcpp::Nullable<Rcpp::NumericVector> ini_reg_par = R_NilValue,
                        const double& theta0 = 0.0,
                        const double& LL_null_ase = 0.0,
                        const arma::uword& max_iter = 4e3,
                        const double& eps = 1e-7, const bool& show = false){
  arma::uword iter1 = 0;
  arma::uword converge = 0;
  arma::uword pp = X.n_cols+fam_nb;
  arma::vec offsets = arma::zeros<arma::vec>(z.n_elem);
  Rcpp::List new_reg, new_bxj_fit, new_theta_fit, curr_theta_fit,
  ase_fit, trec_fit;
  double new_LL, new_bxj, phi, curr_LL, LL0;
  double curr_bxj = 0.0;
  arma::vec curr_reg_par = arma::zeros<arma::vec>(pp);
  arma::vec new_reg_par = arma::zeros<arma::vec>(pp);
  arma::vec BETA = arma::zeros<arma::vec>(pp-1);
  arma::vec parAse = arma::zeros<arma::vec>(2);
  double new_lg_theta  = 0.1;
  double curr_lg_theta = 0.1;
  
  // if(show){
  //   printR_obj("begin ase fit and initial theta fit");
  // }
  //ase fit (initial theta fit)
  ase_fit = Rcpp_ase(ni, ni0, zeta, lbc, theta0, LL_null_ase,
                     max_iter, eps, false);
  double aselrt = as<double>(ase_fit["lrt"]);
  curr_lg_theta = as<double>(ase_fit["par0"]);
  LL0 = as<double>(ase_fit["LL0"]);
  
  // if(show){
  //   printR_obj("begin trec fit and initial NB regression fit");
  // }
  
  trec_fit = Rcpp_trec(y, X, z, fam_nb, lgy1, ini_bxj, LL_null, ini_reg_par,
                       max_iter, eps, false);
  curr_bxj = as<double>(trec_fit["bxj"]);
  curr_reg_par = as<arma::vec>(trec_fit["reg_par"]);
  LL0     += as<double>(trec_fit["LL0"]);
  
  curr_LL = LL0;
  BETA    = curr_reg_par.subvec(0, X.n_cols-1);
  phi     = std::exp(curr_reg_par.at(pp-1));
  
  if(as<int>(ase_fit["converge"]) !=1 | as<int>(trec_fit["converge"]) !=1){
    return Rcpp::List::create(
      Rcpp::Named("ASE_par", ase_fit["par"]),
      Rcpp::Named("ASE_lrt", aselrt),
      Rcpp::Named("ASE_pval", R::pchisq(aselrt, 1, 0, 0)),
      Rcpp::Named("ASE_converge", as<int>(ase_fit["converge"]) ),
      Rcpp::Named("Trec_converge", as<int>(trec_fit["converge"])),
      Rcpp::Named("reg_par", Rcpp::NumericVector(new_reg_par.begin(), new_reg_par.end())),
      Rcpp::Named("converge", 0),
      Rcpp::Named("Trec_lrt", trec_fit["lrt"]),
      Rcpp::Named("Trec_pval", trec_fit["pvalue"]),
      Rcpp::Named("Trec_bxj", curr_bxj),
      Rcpp::Named("Trec_reg_par", trec_fit["reg_par"]));
  }

  //printR_obj(curr_reg_par);
  while(iter1 < max_iter){
    
    //update bxj trecase
    new_bxj_fit = Rcpp_trecase_BFGS(curr_bxj, y, X, z, BETA, phi, fam_nb, lgy1,
                                    ni, ni0, curr_lg_theta, lbc, zeta,
                                    max_iter, eps, false);
    new_bxj     = as<double>(new_bxj_fit["PAR"]);
    new_LL      = new_bxj_fit["LL"];
    
    // if(new_LL < curr_LL - eps) {
    //   printR_obj("likelihood decreased");
    //   printR_obj(new_LL);
    //   printR_obj(curr_LL);
    //   printR_obj(as<int>(new_bxj_fit["converge"]));
    //   converge = 0;
    //   break;
    // }
    
    // if(show){
    //   Rprintf("bxj updated after %d iter \n", as<int>(new_bxj_fit["iter"]));
    // }
    
    // update theta (if theta > 1e-5)
    // double Pi1 = std::exp(new_bxj)/(1.0 + std::exp(new_bxj));
    if(new_lg_theta > -11.5){
      new_theta_fit = Rcpp_ase_theta_BFGS(ni, ni0, zeta, new_bxj, new_lg_theta,
                                        lbc, max_iter, eps, false);
      new_lg_theta =  as<double>(new_theta_fit["PAR"]);
    }
    // if(show){
    //   Rprintf("theta updated after %d iter \n", as<int>(new_theta_fit["iter"]));
    // }
    
    //update BETA, phi in trecase
    compute_offset(new_bxj, z, offsets);
    new_reg     = Rcpp_reg_BFGS(y, X, offsets, curr_reg_par, fam_nb, lgy1,
                                max_iter, eps, false);
    new_reg_par = as<arma::vec>(new_reg["PAR"]);
    BETA = new_reg_par.subvec(0, X.n_cols-1);
    phi  = std::exp(new_reg_par.at(pp-1));
    
    // if(show){
    //   Rprintf("BETA, PHI updated after %d iter \n", as<int>(new_reg["iter"]));
    // }
    
    if(iter1 > 0){
      if( std::abs(curr_LL - new_LL) < eps &&
          Rcpp_norm(curr_reg_par - new_reg_par) < eps &&
          std::abs(curr_bxj - new_bxj) < eps &&
          std::abs(curr_lg_theta - new_lg_theta) < eps){
        if(((curr_LL - new_LL > 0.0 && 
           as<double>(new_bxj_fit["norm_GRAD"]) > 0.01)) ){
          // magnitude of log-like decrease is samll but gradient is large
          // still consider not converged
          converge = 0;
        }else{
          converge = 1;
        }
        break;
        
      }
    }
    curr_reg_par = new_reg_par;
    curr_bxj     = new_bxj;
    curr_LL      = new_LL;
    curr_lg_theta = new_lg_theta;
    iter1++;
    //rintR_obj(iter1);
  }
  double lrt = (new_LL-LL0)*2.0 ;
  double CisTrans_lrt = -2.0*(new_LL - as<double>(ase_fit["LL"]) -
                              as<double>(trec_fit["LL"]));
  if(show){
    Rprintf("TReCASE done after %d iter \n", iter1);
  }  
  return Rcpp::List::create(
    Rcpp::Named("ASE_par", ase_fit["par"]),
    Rcpp::Named("ASE_lrt", aselrt),
    Rcpp::Named("ASE_pval", R::pchisq(aselrt, 1, 0, 0)),
    Rcpp::Named("ASE_converge", as<int>(ase_fit["converge"]) ),
    Rcpp::Named("Trec_converge", as<int>(trec_fit["converge"])),
    Rcpp::Named("bxj", new_bxj),
    Rcpp::Named("lg_theta", new_lg_theta),
    Rcpp::Named("reg_par", Rcpp::NumericVector(new_reg_par.begin(), new_reg_par.end())),
    // Rcpp::Named("LL0", LL0),
    // Rcpp::Named("LL", new_LL),
    Rcpp::Named("lrt", lrt),
    Rcpp::Named("pval", R::pchisq(lrt, 1, 0, 0)),
    Rcpp::Named("converge", converge),
    Rcpp::Named("Trec_lrt", trec_fit["lrt"]),
    Rcpp::Named("Trec_pval", trec_fit["pvalue"]),
    Rcpp::Named("Trec_bxj", trec_fit["bxj"]),
    Rcpp::Named("Trec_reg_par", trec_fit["reg_par"]),
    Rcpp::Named("CisTrans_lrt", CisTrans_lrt),
    Rcpp::Named("CisTrans_pval", R::pchisq(CisTrans_lrt, 1, 0, 0)),
    Rcpp::Named("iter1", iter1)
  );
}


/* ---------------------------
 * SNP-GENE (converge: 1-converged, 0-not converge)
 ---------------------------*/

/*
 * for each gene choose the most significant snp to run trecase 
 */
// // [[Rcpp::export]]
// void Rcpp_trecase_mtest(const arma::mat& Y, const arma::mat& Y1,
//                         const arma::mat& Y2, const arma::mat& Z,
//                         const arma::mat& XX, const arma::vec& SNP_pos,
//                         const arma::uvec& sChr, const bool& fam_nb,
//                         const arma::vec& gene_start, const arma::vec& gene_end,
//                         const arma::uvec& gChr,
//                         const char* file_trec = "trec.txt",
//                         const char* file_trecase = "trecase.txt",
//                         const double& cis_window=1e5,
//                         const bool& useASE = 1, const arma::uword& min_ASE_total=8,
//                         const arma::uword& min_nASE=10, const double& eps=1e-5,
//                         const arma::uword& max_iter=4000L,const bool& show=false){
//   arma::uword gg, ss, ii, ssMin, xi;
//   arma::uword ssBegin = 0;
//   double nSam = Y.n_rows;
//   double pp = XX.n_cols+fam_nb;
//   //arma::vec y  = arma::zeros<arma::vec>(nSam);
//   //arma::vec y1 = arma::zeros<arma::vec>(nSam);
//   //arma::vec y2 = arma::zeros<arma::vec>(nSam);
//   Rcpp::List res_trec, res_trecase;
// 
//   //create files for TReC and TReCASE results
//   FILE * f1, * f2;
// 
//   f1 = fopen(file_trec, "w");
// 
//   fprintf(f1, "GeneRowID\tMarkerRowID\tTReC_b\tTReC_Chisq\tTReC_Pvalue\tTreC_Conv\t");
//   for(xi=0;xi<XX.n_cols;xi++){
//     fprintf(f1, "beta%d\t", xi);
//   }
//   if(fam_nb){
//     fprintf(f1, "phi\n");
//   }
// 
//   if(useASE){
// 
//     f2 = fopen(file_trecase, "w");
// 
//     fprintf(f2, "GeneRowID\tMarkerRowID\tTReC_b\tTReC_Chisq\tTReC_Pvalue\t");
//     for(xi=0;xi<XX.n_cols;xi++){
//       fprintf(f2, "TReC_beta%d\t", xi);
//     }
//     if(fam_nb){
//       fprintf(f2, "TReC_phi\t");
//     }
//     fprintf(f2, "Joint_b\tJoint_Chisq\tJoint_Pvalue\t");
//     for(xi=0;xi<XX.n_cols;xi++){
//       fprintf(f2, "Joint_beta%d\t", xi);
//     }
//     if(fam_nb){
//       fprintf(f2, "Joint_phi\t");
//     }
//     fprintf(f2, "Joint_theta\tConverge\tCisTrans_Chisq\tCisTrans_Pvalue\n");
// 
//   }
// 
//   //gene-snp pair
//   for(gg=0; gg<Y.n_cols; gg++){
// 
//     Rprintf("gene %d \n", gg+1);
//     arma::vec y  = Y.col(gg);
// 
//     double ini_bxj = 0.0;
//     double ptmp = 1.0;
// 
//     arma::vec lgy1 = Rcpp_lgy_add_1(y); //lgamma(y + 1)
//     arma::vec offsets  = arma::zeros<arma::vec>(nSam);
//     arma::vec curr_reg_par  = arma::zeros<arma::vec>(pp);
// 
//     //null trec model for one gene
//     //printR_obj("fit TReC null model");
//     Rcpp::List new_reg = Rcpp_reg_BFGS(y, XX, offsets, curr_reg_par, fam_nb, lgy1,
//                                        max_iter, eps, false);
//     Rcpp::NumericVector ini_reg_par = new_reg["PAR"];
//     double LL_null     = as<double>(new_reg["LL"]);
//     //Rprintf("LL_null %.4f \n", LL_null);
// 
//     if(as<int>(new_reg["converge"]) != 1){
//         Rprintf(" baseline TReC for gene %d model does not converge \n",
//                 gg+1);
//         continue;
//     }
// 
//     for(ss = ssBegin; ss < Z.n_cols ; ss++){
// 
//       if(gChr.at(gg) != sChr.at(ss)){
//         ssBegin = ss;
//         break;
//       }
//       //chr : geneInfo and SNPinfo have to be in order
//       // ssChr ++;
// 
//       if(SNP_pos.at(ss) > gene_start.at(gg) - cis_window &&
//          SNP_pos.at(ss) < gene_end.at(gg)   + cis_window){
// 
//         arma::vec zz2 = Z.col(ss);
//         arma::vec zz  = Z.col(ss);
// 
//         for(ii=0;ii<nSam;ii++){
//           if(zz2.at(ii)==2){
//             zz.at(ii) = 1;
//           }else if(zz2.at(ii)==3){
//             zz.at(ii) = 2;
//           }
//         }
// 
//         res_trec = Rcpp_trec(y, XX, zz, fam_nb, lgy1, ini_bxj, LL_null,
//                              ini_reg_par, max_iter, eps, show);
//         NumericVector reg_pars = res_trec["reg_par"];
//         //ini_bxj = as<double>(res_trec["bxj"]); //initial value of next snp
//         if(as<int>(new_reg["converge"]) != 1){
//           Rprintf("TReC model for snp %d does not converge \n",
//                   ss+1);
//         }
//         //printR_obj(as<double>(res_trec["LL"]));
//         //printR_obj(reg_pars);
// 
//         if(as<double>(res_trec["pvalue"]) < ptmp &&
//            as<int>(res_trec["converge"]) == 1 ){
//           ptmp = as<double>(res_trec["pvalue"]);
//           ssMin = ss;
//         }
// 
//         //write out trec result
//         fprintf(f1, "%d\t%d\t%.2e\t%.2e\t%.4e\t%d\t",
//                 gg+1,ss+1, as<double>(res_trec["bxj"]),
//                 as<double>(res_trec["lrt"]),
//                 as<double>(res_trec["pvalue"]),as<int>(res_trec["converge"]));
//         for(xi=0;xi<(pp-1);xi++){
//           fprintf(f1, "%.2e\t", reg_pars[xi]);
//         }
//         if(fam_nb){
//           fprintf(f1, "%.2e\n", exp(reg_pars[pp-1]));
//         }
// 
//         //printR_obj(as<int>(res_trec["converge"]));
//         //printR_obj("TREC");
//       }
//     }
// 
//     if(useASE){
// 
//       //Trecase
//       arma::vec ni = arma::zeros<arma::vec>(nSam);
//       arma::vec ni0 = arma::zeros<arma::vec>(nSam);
//       arma::vec lbc = arma::zeros<arma::vec>(nSam);
//       arma::vec zeta = arma::zeros<arma::vec>(nSam);
//       arma::uword h1 = 0, h0 = 0;
// 
//       arma::vec zz2 = Z.col(ssMin);
//       arma::vec zz = Z.col(ssMin);
//       arma::vec y1 = Y1.col(gg);
//       arma::vec y2 = Y2.col(gg);
// 
//       for(ii=0;ii<nSam;ii++){
// 
//         double nTi = y1.at(ii) + y2.at(ii);
// 
//         if(zz2.at(ii)==2){
//           zz.at(ii) = 1;
//         }else if(zz2.at(ii)==3){
//           zz.at(ii) = 2;
//         }
// 
//         if(nTi < min_ASE_total){
//           continue;
//         }
// 
//         ni.at(h0) = nTi;
//         if(zz2.at(ii)==0){
//           ni0.at(h0) = y2.at(ii);
//           zeta.at(h0) = 0;
//         }else if(zz2.at(ii)==1){
//           ni0.at(h0) = y2.at(ii);
//           zeta.at(h0) = 1;
//           h1++;
//         }else if(zz2.at(ii)==2){
//           ni0.at(h0) = y1.at(ii);
//           zeta.at(h0) = 1;
//           h1++;
//         }else{
//           ni0.at(h0) = y2.at(ii);
//           zeta.at(h0) = 0;
//         }
// 
//         lbc.at(h0) = R::lchoose(ni.at(h0), ni0.at(h0));
//         h0++;
//       }
//       // printR_obj(h0);
//       // printR_obj(lbc.subvec(0, 10));
//       // printR_obj(ni0.subvec(0, 10));
// 
// 
//       if(h1 < min_nASE){
//         printR_obj("sample size of heterzygous genotype is not enough");
//       }else{
//         res_trecase = Rcpp_trecase(y, XX, zz, fam_nb, lgy1, ni.subvec(0, h0-1),
//                                    ni0.subvec(0, h0-1), zeta.subvec(0, h0-1),
//                                    lbc.subvec(0, h0-1), ini_bxj, LL_null,
//                                    ini_reg_par, max_iter, eps, show) ;
//         NumericVector Trec_reg_pars = res_trecase["Trec_reg_par"];
//         NumericVector Trecase_reg_pars = res_trecase["reg_par"];
// 
//         // write out trec result
//         // printR_obj(as<int>(res_trecase["converge"]));
//         // printR_obj("TRECASE");
// 
//         fprintf(f2, "%d\t%d\t%.2e\t%.2e\t%.4e\t",
//                 gg+1,ssMin+1, as<double>(res_trecase["Trec_bxj"]),
//                 as<double>(res_trecase["Trec_lrt"]),
//                 as<double>(res_trecase["Trec_pval"]));
//         for(xi=0;xi<XX.n_cols;xi++){
//           fprintf(f2, "%.2e\t", Trec_reg_pars[xi]);
//         }
//         if(fam_nb){
//           fprintf(f2, "%.2e\t", exp(Trec_reg_pars[XX.n_cols]));
//         }
//         fprintf(f2, "%.2e\t%.2e\t%.4e\t",
//                 as<double>(res_trecase["bxj"]),
//                 as<double>(res_trecase["lrt"]),
//                 as<double>(res_trecase["pval"]));
//         for(xi=0;xi<XX.n_cols;xi++){
//           fprintf(f2, "%.2e\t", Trecase_reg_pars[xi]);
//         }
//         if(fam_nb){
//           fprintf(f2, "%.2e\t", exp(Trecase_reg_pars[XX.n_cols]));
//         }
//         fprintf(f2, "%.2e\t%d\t%.2e\t%.4e\n",
//                 exp(as<double>(res_trecase["lg_theta"])),
//                 as<int>(res_trecase["converge"]),
//                 as<double>(res_trecase["CisTrans_lrt"]),
//                 as<double>(res_trecase["CisTrans_pval"]) );
// 
//       }
// 
//     }
// 
//   }
//   fclose(f1);
//   if(useASE){
//     fclose(f2);
//   }
// }

// [[Rcpp::export]]
void Rcpp_ase_mtest(const arma::mat& Y1, const arma::mat& Y2,
                    const arma::mat& Z, const arma::vec& SNP_pos,
                    const arma::uvec& sChr,
                    const arma::vec& gene_start, const arma::vec& gene_end,
                    const arma::uvec& gChr,
                    const char* file_ase = "ase.txt", const double& cis_window=1e5,
                    const int& min_ASE_total=8,
                    const int& min_nASE=5, const int& min_nASE_het=5 ,
                    const double& eps=1e-5,
                    const arma::uword& max_iter=4000L,const bool& show=false){
  arma::uword gg, ss, ii;
  arma::uword ssBegin = 0;
  Rcpp::List res_ase;
  FILE * f2;
  int nSam  = Y1.n_rows;
  int nGene = Y1.n_cols;
  int aseNull_conv;
  double theta0, LL_Null;
  
  f2 = fopen(file_ase, "w");
  fprintf(f2,"GeneRowID\tMarkerRowID\tASE_b\tASE_theta\tASE_Chisq\tASE_Pvalue\tASE_Conv\n");
  
  for(gg=0; gg<nGene; gg++){
    Rprintf("gene %d \n", gg+1);
    
    arma::vec y1 = Y1.col(gg);
    arma::vec y2 = Y2.col(gg);
    arma::vec yA = y1 + y2;
    arma::vec zeta = arma::zeros<arma::vec>(nSam);
    arma::vec ni = arma::zeros<arma::vec>(nSam);
    arma::vec lbc = arma::zeros<arma::vec>(nSam);
    arma::vec ni0 = arma::zeros<arma::vec>(nSam);
    arma::uword h0 = 0;
    
    for(ii=0;ii<nSam;ii++){
      if(yA.at(ii) >= min_ASE_total){
        ni.at(h0) = yA.at(ii);
        // lbc.at(h0) = R::lchoose(ni.at(h0), y1.at(ii));
        ni0.at(h0) = y1.at(ii);
        h0++;
      }
    }
    if(h0 < min_nASE){
      aseNull_conv = 0;
    }else{
      Rcpp::List ase_Null = Rcpp_ase_theta_BFGS(ni.subvec(0, h0-1),
                                              ni0.subvec(0, h0-1), 
                                              zeta.subvec(0, h0-1),
                                              0.0, -2.0, lbc.subvec(0, h0-1), 
                                              max_iter, eps, show);
      theta0 = exp(as<double>(ase_Null["PAR"]));
      LL_Null = as<double>(ase_Null["LL"]);
      aseNull_conv = as<int>(ase_Null["converge"]);
    }
    if(aseNull_conv != 1){
      // Rprintf("baseline ASE for gene %d model does not converge \n",
      //         gg+1);
      continue;
    }
    
    
    for(ss = ssBegin; ss < Z.n_cols ; ss++){
      
      if(gChr.at(gg) != sChr.at(ss)){
        ssBegin = ss;
        break;
      }
      //chr : geneInfo and SNPinfo have to be in order
      // ssChr ++;
      
      if(SNP_pos.at(ss) > gene_start.at(gg) - cis_window &&
         SNP_pos.at(ss) < gene_end.at(gg)   + cis_window){
        
        arma::vec zz2 = Z.col(ss);
        arma::vec zeta = arma::zeros<arma::vec>(nSam);
        arma::uword h1 = 0, h0 = 0;
        
        for(ii=0;ii<nSam;ii++){
          if(yA.at(ii) >= min_ASE_total){
            if(zz2.at(ii)==1){
              ni0.at(h0) = y2.at(ii);
              zeta.at(h0) = 1;
              h1++;
            }else if(zz2.at(ii)==2){
              ni0.at(h0) = y1.at(ii);
              zeta.at(h0) = 1;
              h1++;
            }
            h0++;
          }
        }
        
        if(h0 < min_nASE){
          // Rprintf("sample size of heterzygous genotype is not enoug \n",
          //         ss+1);
        }else{
          //begin ase fit
          res_ase = Rcpp_ase(ni.subvec(0, h0-1), ni0.subvec(0, h0-1),
                             zeta.subvec(0, h0-1), lbc.subvec(0, h0-1),
                             theta0, LL_Null,
                             max_iter, eps, show);
          NumericVector ase_pars = res_ase["par"];
          fprintf(f2, "%d\t%d\t%.2e\t%.2e\t%.4e\t%.4e\t%d\n",
                  gg+1,ss+1, ase_pars[0], ase_pars[1],
                                                  as<double>(res_ase["lrt"]),
                                                  as<double>(res_ase["pvalue"]),
                                                  as<int>(res_ase["converge"]));
          
        }
        
      }
    }
  }
  fclose(f2);
  
}

/*
 * run all gene-snp pair for trec / trecase 
 * Fail code: 0: TReCASE does not converge.
 * 2: ASE does not converge, 4: TReC does not converge, 
 * 6: both ASE and TReC do not converge. 
 * 7: not enouth AS sample / heterzegous SNPs 
 * 8: ASE null Fail 
 * 10: both ASE null and TReC not converge
 */
// [[Rcpp::export]]
void Rcpp_trecase_mtest(const arma::mat& Y, const arma::mat& Y1,
                        const arma::mat& Y2, const arma::mat& Z,
                        const arma::mat& XX, const arma::vec& SNP_pos,
                        const arma::uvec& sChr, const bool& fam_nb,
                        const arma::vec& gene_start, const arma::vec& gene_end,
                        const arma::uvec& gChr,
                        const List& GeneSnpList, 
                        const char* file_trec = "trec.txt",
                        const char* file_trecase = "trecase.txt",
                        const double& transTestP = 0.01,
                        const double& cis_window=1e5,
                        const bool& useASE = 1, const int& min_ASE_total=8,
                        const int& min_nASE= 5, const int& min_nASE_het=5,
                        const double& eps=1e-5,
                        const arma::uword& max_iter=4000L,const bool& show=false){
  arma::uword gg, ss, ii, xi, h0, h1;
  arma::uword ssBegin = 0;
  double nSam = Y.n_rows;
  double pp = XX.n_cols+fam_nb;
  double theta0, LL_null, LL_null_ase, ini_bxj;
  int aseNull_conv;
  arma::vec y       = arma::zeros<arma::vec>(nSam);
  arma::vec y1      = arma::zeros<arma::vec>(nSam);
  arma::vec y2      = arma::zeros<arma::vec>(nSam);
  
  arma::vec lgy1    = arma::zeros<arma::vec>(nSam);
  arma::vec offsets = arma::zeros<arma::vec>(nSam);
  arma::vec yA      = arma::zeros<arma::vec>(nSam);
  arma::vec zeta    = arma::zeros<arma::vec>(nSam);
  arma::vec ni      = arma::zeros<arma::vec>(nSam);
  arma::vec lbc     = arma::zeros<arma::vec>(nSam);
  arma::vec ni0     = arma::zeros<arma::vec>(nSam);
  
  arma::vec curr_reg_par  = arma::zeros<arma::vec>(pp);
  Rcpp::List res_trec, res_trecase, ase_Null, new_reg;
  //create files for TReC and TReCASE results
  FILE * f1, * f2;
  
  if(useASE){
    
    f2 = fopen(file_trecase, "w");
    
    fprintf(f2, "GeneRowID\tMarkerRowID\tTReC_b\tTReC_Chisq\tTReC_Pvalue\t");
    for(xi=0;xi<XX.n_cols;xi++){
      fprintf(f2, "TReC_beta%d\t", xi);
    }
    if(fam_nb){
      fprintf(f2, "TReC_phi\t");
    }
    fprintf(f2,"ASE_b\tASE_theta\tASE_Chisq\tASE_Pvalue\t");
    fprintf(f2, "Joint_b\tJoint_Chisq\tJoint_Pvalue\t");
    for(xi=0;xi<XX.n_cols;xi++){
      fprintf(f2, "Joint_beta%d\t", xi);
    }
    if(fam_nb){
      fprintf(f2, "Joint_phi\t");
    }
    fprintf(f2, "Joint_theta\tConverge\tCisTrans_Chisq\tCisTrans_Pvalue\tp.final\n");
    
  }else{
    
    f1 = fopen(file_trec, "w");
    
    fprintf(f1, "GeneRowID\tMarkerRowID\tTReC_b\tTReC_Chisq\tTReC_Pvalue\tTReC_Conv\t");
    for(xi=0;xi<XX.n_cols;xi++){
      fprintf(f1, "beta%d\t", xi);
    }
    if(fam_nb){
      fprintf(f1, "phi\n");
    }
  }
  
  //gene-snp pair
  if(GeneSnpList.length() > 0){
    for(gg=0; gg<GeneSnpList.length(); gg++){
      
      if(GeneSnpList[gg]==R_NilValue){
        continue;
      }
      Rprintf("gene %d \n", gg+1);
      y  = Y.col(gg);
      
      ini_bxj = 0.0;
      lgy1 = Rcpp_lgy_add_1(y); //lgamma(y + 1)
      offsets.zeros();
      curr_reg_par.zeros();
      
      //null trec model for one gene
      //printR_obj("fit TReC null model");
      new_reg = Rcpp_reg_BFGS(y, XX, offsets, curr_reg_par, fam_nb, lgy1,
                              max_iter, eps, false);
      Rcpp::NumericVector ini_reg_par = new_reg["PAR"];
      LL_null     = as<double>(new_reg["LL"]);
      //Rprintf("LL_null %.4f \n", LL_null);
      
      if(as<int>(new_reg["converge"]) != 1){
        Rprintf("baseline TReC for gene %d model does not converge \n",
                gg+1);
        continue;
      }
      
      if(useASE){
        y1 = Y1.col(gg);
        y2 = Y2.col(gg);
        yA = y1 + y2;
        zeta = arma::zeros<arma::vec>(nSam);
        ni = arma::zeros<arma::vec>(nSam);
        lbc = arma::zeros<arma::vec>(nSam);
        ni0 = arma::zeros<arma::vec>(nSam);
        h0 = 0;
        
        for(ii=0;ii<nSam;ii++){
          
          if(yA.at(ii) >= min_ASE_total){
            ni.at(h0) = yA.at(ii);
            // lbc.at(h0) = R::lchoose(ni.at(h0), y1.at(ii));
            ni0.at(h0) = y1.at(ii);
            h0++;
          }
        }
        if(h0 < min_nASE){
          aseNull_conv = 7;
        }else{
          Rcpp::List ase_Null = Rcpp_ase_theta_BFGS(ni.subvec(0, h0-1),
                                                    ni0.subvec(0, h0-1), 
                                                    zeta.subvec(0, h0-1),
                                                    0.0, -2.0, lbc.subvec(0, h0-1), 
                                                    max_iter, eps, show);
          theta0 = exp(as<double>(ase_Null["PAR"]));
          LL_null_ase = as<double>(ase_Null["LL"]);
          aseNull_conv = as<int>(ase_Null["converge"]) ;
        }
        
        // if(aseNull_conv != 1){
        //   Rprintf("baseline ASE for gene %d model does not converge \n",
        //           gg+1);
        //   continue;
        // }
      }
      
      arma::vec ssVec = GeneSnpList[gg]; 
      
      for(ssBegin = 0; ssBegin < ssVec.n_elem; ssBegin++){
        
        ss = ssVec.at(ssBegin)-1;
        arma::vec zz2 = Z.col(ss);
        arma::vec zz  = Z.col(ss);
        
        if(useASE & aseNull_conv == 1){
          
          arma::vec zz2 = Z.col(ss);
          zeta.zeros();
          arma::uword h1 = 0, h0 = 0;
          
          for(ii=0;ii<nSam;ii++){
            if(zz2.at(ii)==2){
              zz.at(ii) = 1;
            }else if(zz2.at(ii)==3){
              zz.at(ii) = 2;
            }
            
            if(yA.at(ii) >= min_ASE_total){
              if(zz2.at(ii)==1){
                ni0.at(h0) = y2.at(ii);
                zeta.at(h0) = 1;
                h1++;
              }else if(zz2.at(ii)==2){
                ni0.at(h0) = y1.at(ii);
                zeta.at(h0) = 1;
                h1++;
              }
              h0++;
            }
          }
          
          if(h1 < min_nASE_het){
            // can still run trec model
            // Rprintf("sample size of heterzygous genotype is not enough for SNP %d\n",
            //         ss+1);
            
            res_trec = Rcpp_trec(y, XX, zz, fam_nb, lgy1, ini_bxj, LL_null,
                                 ini_reg_par, max_iter, eps, show);
            NumericVector reg_pars = res_trec["reg_par"];
            // if(as<int>(new_reg["converge"]) != 1){
            //   Rprintf("TReC model for snp %d does not converge \n",
            //           ss+1);
            // }
            
            fprintf(f2, "%d\t%d\t%.2e\t%.2e\t%.4e\t",
                    gg+1,ss+1, as<double>(res_trec["bxj"]),
                    as<double>(res_trec["lrt"]),
                    as<double>(res_trec["pvalue"]));
            for(xi=0;xi<(pp-1);xi++){
              fprintf(f2, "%.2e\t", reg_pars[xi]);
            }
            if(fam_nb){
              fprintf(f2, "%.2e\t", exp(reg_pars[pp-1]));
            }
            fprintf(f2,"NA\tNA\tNA\tNA\t");
            fprintf(f2, "NA\tNA\tNA\t");
            for(xi=0;xi<XX.n_cols;xi++){
              fprintf(f2, "NA\t");
            }
            if(fam_nb){
              fprintf(f2, "NA\t");
            }
            fprintf(f2, "NA\t%d\tNA\tNA\t%.4e\n", 
                    4+3*as<int>(res_trec["converge"]),
                    as<double>(res_trec["pvalue"]));
            
            
          }else{
            res_trecase = Rcpp_trecase(y, XX, zz, fam_nb, lgy1, ni.subvec(0, h0-1),
                                       ni0.subvec(0, h0-1), zeta.subvec(0, h0-1),
                                       lbc.subvec(0, h0-1), ini_bxj, LL_null,
                                       ini_reg_par,theta0, LL_null_ase,
                                       max_iter, eps, show) ;
            NumericVector Trec_reg_pars = res_trecase["Trec_reg_par"];
            NumericVector Trecase_reg_pars = res_trecase["reg_par"];
            NumericVector ase_par = res_trecase["ASE_par"];
            int trecase_conver = as<int>(res_trecase["converge"]);
            double pfinal = as<double>(res_trecase["Trec_pval"]);
            if(as<int>(res_trecase["ASE_converge"]) != 1){
              // Rprintf("ASE model for snp %d does not converge \n",
              //         ss+1);
              trecase_conver += 2;
            }
            if(as<int>(res_trecase["Trec_converge"]) != 1){
              // Rprintf("TReC model for snp %d does not converge \n",
              //         ss+1);
              trecase_conver += 4;
              
            }
            
            // write out trec result
            // printR_obj(as<int>(res_trecase["converge"]));
            // printR_obj("TRECASE");
            fprintf(f2, "%d\t%d\t%.2e\t%.2e\t%.4e\t",
                    gg+1,ss+1, as<double>(res_trecase["Trec_bxj"]),
                    as<double>(res_trecase["Trec_lrt"]),
                    as<double>(res_trecase["Trec_pval"]));
            for(xi=0;xi<XX.n_cols;xi++){
              fprintf(f2, "%.2e\t", Trec_reg_pars[xi]);
            }
            if(fam_nb){
              fprintf(f2, "%.2e\t", exp(Trec_reg_pars[XX.n_cols]));
            }
            fprintf(f2,"%.2e\t%.2e\t%.2e\t%.4e\t", 
                    ase_par[0], exp(ase_par[1]), 
                    as<double>(res_trecase["ASE_lrt"]), 
                    as<double>(res_trecase["ASE_pval"]) );
            if(trecase_conver != 1){
              fprintf(f2, "NA\tNA\tNA\t");
              for(xi=0;xi<XX.n_cols;xi++){
                fprintf(f2, "NA\t");
              }
              if(fam_nb){
                fprintf(f2, "NA\t");
              }
              fprintf(f2, "NA\t%d\tNA\tNA\t%.4e\n", trecase_conver,
                      as<double>(res_trecase["Trec_pval"]));
            }else{  
              
              if(as<double>(res_trecase["CisTrans_pval"]) >= transTestP){
                pfinal = as<double>(res_trecase["pval"]);
              }
              fprintf(f2, "%.2e\t%.2e\t%.4e\t",
                      as<double>(res_trecase["bxj"]),
                      as<double>(res_trecase["lrt"]),
                      as<double>(res_trecase["pval"]));
              for(xi=0;xi<XX.n_cols;xi++){
                fprintf(f2, "%.2e\t", Trecase_reg_pars[xi]);
              }
              if(fam_nb){
                fprintf(f2, "%.2e\t", exp(Trecase_reg_pars[XX.n_cols]));
              }
              fprintf(f2, "%.2e\t%d\t%.2e\t%.4e\t%.4e\n",
                      exp(as<double>(res_trecase["lg_theta"])),
                      trecase_conver,
                      as<double>(res_trecase["CisTrans_lrt"]),
                      as<double>(res_trecase["CisTrans_pval"]),
                      pfinal);
              
            }}
          
        }else{
          // run trec model (useASE = 0 or didn't pass null model ASE )
          for(ii=0;ii<nSam;ii++){
            if(zz2.at(ii)==2){
              zz.at(ii) = 1;
            }else if(zz2.at(ii)==3){
              zz.at(ii) = 2;
            }
          }
          
          res_trec = Rcpp_trec(y, XX, zz, fam_nb, lgy1, ini_bxj, LL_null,
                               ini_reg_par, max_iter, eps, show);
          NumericVector reg_pars = res_trec["reg_par"];
          int Conv2 = res_trec["converge"];
          //ini_bxj = as<double>(res_trec["bxj"]); //initial value of next snp
          // 
          // if(as<int>(res_trec["converge"]) != 1){
          //    Rprintf("TReC model for snp %d does not converge \n",
          //            ss+1);
          //    }
          
          if(aseNull_conv == 0){
            Conv2 = 10 - 2*Conv2;
          }else{
            Conv2 = 4 + 3*Conv2;
          }
          //printR_obj(as<double>(res_trec["LL"]));
          //printR_obj(reg_pars);
          
          if(useASE){ //write out trecase result 
            // (when ase null model is not coverged)
            fprintf(f2, "%d\t%d\t%.2e\t%.2e\t%.4e\t",
                    gg+1,ss+1, as<double>(res_trec["bxj"]),
                    as<double>(res_trec["lrt"]),
                    as<double>(res_trec["pvalue"]));
            for(xi=0;xi<(pp-1);xi++){
              fprintf(f2, "%.2e\t", reg_pars[xi]);
            }
            if(fam_nb){
              fprintf(f2, "%.2e\t", exp(reg_pars[pp-1]));
            }
            fprintf(f2,"NA\tNA\tNA\tNA\t");
            fprintf(f2, "NA\tNA\tNA\t");
            for(xi=0;xi<XX.n_cols;xi++){
              fprintf(f2, "NA\t");
            }
            if(fam_nb){
              fprintf(f2, "NA\t");
            }
            fprintf(f2, "NA\t%d\tNA\tNA\t%.4e\n", 
                    Conv2, 
                    as<double>(res_trec["pvalue"]));
            
          }else{ //write out trec result
            fprintf(f1, "%d\t%d\t%.2e\t%.2e\t%.4e\t%d\t",
                    gg+1,ss+1, as<double>(res_trec["bxj"]),
                    as<double>(res_trec["lrt"]),
                    as<double>(res_trec["pvalue"]),
                    as<int>(res_trec["converge"]));
            for(xi=0;xi<(pp-1);xi++){
              fprintf(f1, "%.2e\t", reg_pars[xi]);
            }
            if(fam_nb){
              fprintf(f1, "%.2e\n", exp(reg_pars[pp-1]));
            }
          }
          
          
          //printR_obj(as<int>(res_trec["converge"]));
          //printR_obj("TREC");
        }
        
        
      }
      
    }
    
  }else{
    for(gg=0; gg<Y.n_cols; gg++){
    
    Rprintf("gene %d \n", gg+1);
    y  = Y.col(gg);
    
    
    ini_bxj = 0.0;
    lgy1 = Rcpp_lgy_add_1(y); //lgamma(y + 1)
    offsets.zeros();
    curr_reg_par.zeros();
    
    //null trec model for one gene
    //printR_obj("fit TReC null model");
    new_reg = Rcpp_reg_BFGS(y, XX, offsets, curr_reg_par, fam_nb, lgy1,
                            max_iter, eps, false);
    Rcpp::NumericVector ini_reg_par = new_reg["PAR"];
    LL_null     = as<double>(new_reg["LL"]);
    //Rprintf("LL_null %.4f \n", LL_null);
    
    if(as<int>(new_reg["converge"]) != 1){
      Rprintf("baseline TReC for gene %d model does not converge \n",
              gg+1);
      continue;
    }
    
    if(useASE){
      y1 = Y1.col(gg);
      y2 = Y2.col(gg);
      yA = y1 + y2;
      zeta = arma::zeros<arma::vec>(nSam);
      ni = arma::zeros<arma::vec>(nSam);
      lbc = arma::zeros<arma::vec>(nSam);
      ni0 = arma::zeros<arma::vec>(nSam);
      h0 = 0;
      
      for(ii=0;ii<nSam;ii++){
        
        if(yA.at(ii) >= min_ASE_total){
          ni.at(h0) = yA.at(ii);
          // lbc.at(h0) = R::lchoose(ni.at(h0), y1.at(ii));
          ni0.at(h0) = y1.at(ii);
          h0++;
        }
      }
      if(h0 < min_nASE){
        aseNull_conv = 7;
      }else{
        Rcpp::List ase_Null = Rcpp_ase_theta_BFGS(ni.subvec(0, h0-1),
                                                  ni0.subvec(0, h0-1), 
                                                  zeta.subvec(0, h0-1),
                                                  0.0, -2.0, lbc.subvec(0, h0-1), 
                                                  max_iter, eps, show);
        theta0 = exp(as<double>(ase_Null["PAR"]));
        LL_null_ase = as<double>(ase_Null["LL"]);
        aseNull_conv = as<int>(ase_Null["converge"]) ;
      }
      
      // if(aseNull_conv != 1){
      //   Rprintf("baseline ASE for gene %d model does not converge \n",
      //           gg+1);
      //   continue;
      // }
    }
    
    for(ss = ssBegin; ss < Z.n_cols ; ss++){
      
      if(gChr.at(gg) != sChr.at(ss)){
        ssBegin = ss;
        break;
      }
      if(ss % 5000 == 0 & show){
        Rprintf("Begin analysis for SNP %d  \n", ss+1);
      }
      //chr : geneInfo and SNPinfo have to be in order
      // ssChr ++;
      
      if(SNP_pos.at(ss) > gene_start.at(gg) - cis_window &&
         SNP_pos.at(ss) < gene_end.at(gg)   + cis_window){
        
        arma::vec zz2 = Z.col(ss);
        arma::vec zz  = Z.col(ss);
        
        if(useASE & aseNull_conv == 1){
          
          arma::vec zz2 = Z.col(ss);
          zeta.zeros();
          arma::uword h1 = 0, h0 = 0;
          
          for(ii=0;ii<nSam;ii++){
            if(zz2.at(ii)==2){
              zz.at(ii) = 1;
            }else if(zz2.at(ii)==3){
              zz.at(ii) = 2;
            }
            
            if(yA.at(ii) >= min_ASE_total){
              if(zz2.at(ii)==1){
                ni0.at(h0) = y2.at(ii);
                zeta.at(h0) = 1;
                h1++;
              }else if(zz2.at(ii)==2){
                ni0.at(h0) = y1.at(ii);
                zeta.at(h0) = 1;
                h1++;
              }
              h0++;
            }
          }
          
          if(h1 < min_nASE_het){
            // can still run trec model
            // Rprintf("sample size of heterzygous genotype is not enough for SNP %d\n",
            //         ss+1);
            
            res_trec = Rcpp_trec(y, XX, zz, fam_nb, lgy1, ini_bxj, LL_null,
                                 ini_reg_par, max_iter, eps, show);
            NumericVector reg_pars = res_trec["reg_par"];
            // if(as<int>(new_reg["converge"]) != 1){
            //   Rprintf("TReC model for snp %d does not converge \n",
            //           ss+1);
            // }
            
            fprintf(f2, "%d\t%d\t%.2e\t%.2e\t%.4e\t",
                    gg+1,ss+1, as<double>(res_trec["bxj"]),
                    as<double>(res_trec["lrt"]),
                    as<double>(res_trec["pvalue"]));
            for(xi=0;xi<(pp-1);xi++){
              fprintf(f2, "%.2e\t", reg_pars[xi]);
            }
            if(fam_nb){
              fprintf(f2, "%.2e\t", exp(reg_pars[pp-1]));
            }
            fprintf(f2,"NA\tNA\tNA\tNA\t");
            fprintf(f2, "NA\tNA\tNA\t");
            for(xi=0;xi<XX.n_cols;xi++){
              fprintf(f2, "NA\t");
            }
            if(fam_nb){
              fprintf(f2, "NA\t");
            }
            fprintf(f2, "NA\t%d\tNA\tNA\t%.4e\n", 
                    4+3*as<int>(res_trec["converge"]),
                    as<double>(res_trec["pvalue"]));
            
            
          }else{
            res_trecase = Rcpp_trecase(y, XX, zz, fam_nb, lgy1, ni.subvec(0, h0-1),
                                       ni0.subvec(0, h0-1), zeta.subvec(0, h0-1),
                                       lbc.subvec(0, h0-1), ini_bxj, LL_null,
                                       ini_reg_par,theta0, LL_null_ase,
                                       max_iter, eps, show) ;
            NumericVector Trec_reg_pars = res_trecase["Trec_reg_par"];
            NumericVector Trecase_reg_pars = res_trecase["reg_par"];
            NumericVector ase_par = res_trecase["ASE_par"];
            int trecase_conver = as<int>(res_trecase["converge"]);
            double pfinal = as<double>(res_trecase["Trec_pval"]);
            if(as<int>(res_trecase["ASE_converge"]) != 1){
              // Rprintf("ASE model for snp %d does not converge \n",
              //         ss+1);
              trecase_conver += 2;
            }
            if(as<int>(res_trecase["Trec_converge"]) != 1){
              // Rprintf("TReC model for snp %d does not converge \n",
              //         ss+1);
              trecase_conver += 4;
              
            }
            
            // write out trec result
            // printR_obj(as<int>(res_trecase["converge"]));
            // printR_obj("TRECASE");
            fprintf(f2, "%d\t%d\t%.2e\t%.2e\t%.4e\t",
                    gg+1,ss+1, as<double>(res_trecase["Trec_bxj"]),
                    as<double>(res_trecase["Trec_lrt"]),
                    as<double>(res_trecase["Trec_pval"]));
            for(xi=0;xi<XX.n_cols;xi++){
              fprintf(f2, "%.2e\t", Trec_reg_pars[xi]);
            }
            if(fam_nb){
              fprintf(f2, "%.2e\t", exp(Trec_reg_pars[XX.n_cols]));
            }
            fprintf(f2,"%.2e\t%.2e\t%.2e\t%.4e\t", 
                    ase_par[0], exp(ase_par[1]), 
                    as<double>(res_trecase["ASE_lrt"]), 
                    as<double>(res_trecase["ASE_pval"]) );
            if(trecase_conver != 1){
              fprintf(f2, "NA\tNA\tNA\t");
              for(xi=0;xi<XX.n_cols;xi++){
                fprintf(f2, "NA\t");
              }
              if(fam_nb){
                fprintf(f2, "NA\t");
              }
              fprintf(f2, "NA\t%d\tNA\tNA\t%.4e\n", trecase_conver,
                      as<double>(res_trecase["Trec_pval"]));
            }else{  
              
              if(as<double>(res_trecase["CisTrans_pval"]) >= transTestP){
                pfinal = as<double>(res_trecase["pval"]);
              }
              fprintf(f2, "%.2e\t%.2e\t%.4e\t",
                      as<double>(res_trecase["bxj"]),
                      as<double>(res_trecase["lrt"]),
                      as<double>(res_trecase["pval"]));
              for(xi=0;xi<XX.n_cols;xi++){
                fprintf(f2, "%.2e\t", Trecase_reg_pars[xi]);
              }
              if(fam_nb){
                fprintf(f2, "%.2e\t", exp(Trecase_reg_pars[XX.n_cols]));
              }
              fprintf(f2, "%.2e\t%d\t%.2e\t%.4e\t%.4e\n",
                      exp(as<double>(res_trecase["lg_theta"])),
                      trecase_conver,
                      as<double>(res_trecase["CisTrans_lrt"]),
                      as<double>(res_trecase["CisTrans_pval"]),
                      pfinal);
              
            }}
          
        }else{
          // run trec model (useASE = 0 or didn't pass null model ASE )
          for(ii=0;ii<nSam;ii++){
            if(zz2.at(ii)==2){
              zz.at(ii) = 1;
            }else if(zz2.at(ii)==3){
              zz.at(ii) = 2;
            }
          }
          
          res_trec = Rcpp_trec(y, XX, zz, fam_nb, lgy1, ini_bxj, LL_null,
                               ini_reg_par, max_iter, eps, show);
          NumericVector reg_pars = res_trec["reg_par"];
          int Conv2 = res_trec["converge"];
          //ini_bxj = as<double>(res_trec["bxj"]); //initial value of next snp
         // 
         // if(as<int>(res_trec["converge"]) != 1){
         //    Rprintf("TReC model for snp %d does not converge \n",
         //            ss+1);
         //    }
          
          if(aseNull_conv == 0){
            Conv2 = 10 - 2*Conv2;
          }else{
            Conv2 = 4 + 3*Conv2;
          }
          //printR_obj(as<double>(res_trec["LL"]));
          //printR_obj(reg_pars);
          
          if(useASE){ //write out trecase result 
            // (when ase null model is not coverged)
            fprintf(f2, "%d\t%d\t%.2e\t%.2e\t%.4e\t",
                    gg+1,ss+1, as<double>(res_trec["bxj"]),
                    as<double>(res_trec["lrt"]),
                    as<double>(res_trec["pvalue"]));
            for(xi=0;xi<(pp-1);xi++){
              fprintf(f2, "%.2e\t", reg_pars[xi]);
            }
            if(fam_nb){
              fprintf(f2, "%.2e\t", exp(reg_pars[pp-1]));
            }
            fprintf(f2,"NA\tNA\tNA\tNA\t");
            fprintf(f2, "NA\tNA\tNA\t");
            for(xi=0;xi<XX.n_cols;xi++){
              fprintf(f2, "NA\t");
            }
            if(fam_nb){
              fprintf(f2, "NA\t");
            }
            fprintf(f2, "NA\t%d\tNA\tNA\t%.4e\n", 
                    Conv2, 
                    as<double>(res_trec["pvalue"]));
            
          }else{ //write out trec result
            fprintf(f1, "%d\t%d\t%.2e\t%.2e\t%.4e\t%d\t",
                    gg+1,ss+1, as<double>(res_trec["bxj"]),
                    as<double>(res_trec["lrt"]),
                    as<double>(res_trec["pvalue"]),
                    as<int>(res_trec["converge"]));
            for(xi=0;xi<(pp-1);xi++){
              fprintf(f1, "%.2e\t", reg_pars[xi]);
            }
            if(fam_nb){
              fprintf(f1, "%.2e\n", exp(reg_pars[pp-1]));
            }
          }
          
          
          //printR_obj(as<int>(res_trec["converge"]));
          //printR_obj("TREC");
        }
      }
      
    }
    
  }
  }
  if(useASE){
    fclose(f2);
  }else{
    fclose(f1);
  }
}




