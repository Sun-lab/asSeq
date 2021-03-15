#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <stdio.h>
#include "shared.h"

// [[Rcpp::depends("RcppArmadillo")]]
// 
// template<typename T>
// void printR_obj(const T& obj){
//   Rcpp::Rcout << obj << std::endl;
// }

// ----------------------------------------------------------------------
// Negative Binomial (TREC)
// ----------------------------------------------------------------------
using namespace Rcpp;

// [[Rcpp::export]]
void RcppT_compute_offset(const arma::vec& z, const arma::vec& RHO,
                          double& KAPPA, double& ETA, double& GAMMA,
                          const arma::vec& tau1, const arma::vec& tau2,
                          arma::vec& offsets){
  arma::uword ii;
  for(ii = 0; ii<RHO.size(); ii++){
    if(z.at(ii) == 0){
      offsets.at(ii) = std::log(2*(1 - RHO.at(ii)) +
        (tau1.at(ii)+tau2.at(ii))*RHO.at(ii)*KAPPA);
    }else if(z.at(ii) == 1){
      offsets.at(ii) = std::log((1-RHO.at(ii)) + tau1.at(ii)*RHO.at(ii)*KAPPA +
        (1-RHO.at(ii))*ETA + tau2.at(ii)*RHO.at(ii)*KAPPA*GAMMA);
    }else if(z.at(ii) == 2){
      offsets.at(ii) = std::log((1-RHO.at(ii)) + tau2.at(ii)*RHO.at(ii)*KAPPA +
        (1 - RHO.at(ii))*ETA + tau1.at(ii)*RHO.at(ii)*KAPPA*GAMMA);
    }else{
      offsets.at(ii) = std::log(2*(1 - RHO.at(ii))*ETA +
        (tau1.at(ii)+tau2.at(ii))*RHO.at(ii)*KAPPA*GAMMA);
    }
  }
}

// [[Rcpp::export]]
void RcppT_compute_expXbeta(const arma::mat& X, const arma::vec& BETA, 
                            arma::vec& expXbeta){
  arma::uword ii;
  for(ii = 0; ii< X.n_rows; ii++){
    expXbeta.at(ii) = std::exp(arma::dot(X.row(ii).t(), BETA));
  }
}

/* ---------------------------
 * NB regression (duplicated function)
 */
// [[Rcpp::export]]
double RcppT_reg_LL(const arma::vec& y, const arma::mat& X,
                    const arma::vec& offsets, const arma::vec& PARAMS,
                    const arma::vec& lgy1, arma::vec& mu){
  // duplicate function
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
arma::vec RcppT_reg_grad(const arma::vec& y, const arma::mat& X,
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
Rcpp::List RcppT_reg_BFGS(const arma::vec& y, const arma::mat& X,
                          const arma::vec& offsets, const arma::vec& params0,
                          const arma::vec& lgy1,
                          const arma::uword& max_iter,
                          const double& eps, const bool& show){
  
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
    old_LL = fnscale * RcppT_reg_LL(y, X, offsets, xk, lgy1, mu);
    gr_k   = fnscale * RcppT_reg_grad(y, X, mu, xk);
    p_k    = -1.0 * inv_Bk * gr_k;
    inv_norm_p_k =  1.0 / std::max(1.0, Rcpp_norm(p_k));
    
    //line search for new xk
    for(jj=0; jj<15; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4, jj);
      new_xk    = xk + tmp_alpha * p_k;
      new_LL    = fnscale * RcppT_reg_LL(y, X, offsets, new_xk, lgy1, mu);
      
      if(new_LL < old_LL){ //minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * RcppT_reg_grad(y, X, mu, new_xk) - gr_k;
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
        gr_k = RcppT_reg_grad(y, X, mu, xk);
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
  
  old_LL = RcppT_reg_LL(y, X, offsets, xk, lgy1, mu);
  //gr_k = RcppT_NB_reg_grad(y, X, mu, xk);
  return Rcpp::List::create(
    Rcpp::Named("converge", converge),
    Rcpp::Named("LL", old_LL),
    Rcpp::Named("iter", iter),
    Rcpp::Named("norm_GRAD", Rcpp_norm(gr_k)),
    Rcpp::Named("PAR", Rcpp::NumericVector(xk.begin(), xk.end()))
  );
}

/* ---------------------------
 * TReC eQTL
 */

// [[Rcpp::export]]
double RcppT_loglikNB_KEG(const arma::vec& para, const int& H0, 
                          const arma::vec& y, const arma::vec& z, 
                          const double& phi, const arma::vec& RHO,
                          const arma::vec& tau1, const arma::vec& tau2,
                          const arma::vec& lgy1, const arma::vec& expXbeta,
                          arma::vec& offsets, arma::vec& mu){
  
  //lgy1 = std::lgamma(y + 1)
  double loglik = 0.0;
  arma::uword ii;
  double vphi = 1.0/phi;
  double KAPPA, ETA, GAMMA;
  
  if(H0 == 0){
    KAPPA = exp(para.at(0));
    ETA   = exp(para.at(1));
    GAMMA = exp(para.at(2));
  }else if(H0 == 1){ //para = log(c(KAPPA, GAMMA))
    KAPPA = exp(para.at(0));
    ETA   = 1.0;
    GAMMA = exp(para.at(1));
  }else{ //para = log(c(KAPPA, ETA))
    KAPPA = exp(para.at(0));
    ETA   = exp(para.at(1));
    GAMMA = 1.0;
  }
  
  RcppT_compute_offset(z, RHO, KAPPA, ETA, GAMMA, tau1, tau2, offsets);
  // double lgvphi = lgamma(vphi); //do not neet to compute when update KEG 
  
  for(ii = 0; ii< y.n_elem; ii++){
    mu.at(ii) = expXbeta.at(ii) * exp(offsets.at(ii));
    if(y.at(ii) > 0){
      loglik +=  y.at(ii) * std::log(mu.at(ii));
    }
    loglik += - (vphi+y.at(ii)) * std::log(vphi+mu.at(ii));
  }
  return loglik;
}


// [[Rcpp::export]]
arma::vec RcppT_grad_NB(const arma::vec& para, const int& H0,
                        const arma::vec& y, const arma::vec& z,
                        const arma::vec& RHO, const double& phi,
                        const arma::vec& tau1, const arma::vec& tau2,
                        const arma::vec& expXbeta, 
                        const arma::vec& offsets, const arma::vec& mu){
  
  double dltrec_dmu = 0, dmu_dkappa = 0, dmu_deta = 0, dmu_dgamma = 0;
  arma::uword ii;
  arma::vec grad = arma::zeros<arma::vec>(para.size());
  double KAPPA, ETA, GAMMA;
  
  if(H0 == 0){
    KAPPA = exp(para.at(0));
    ETA   = exp(para.at(1));
    GAMMA = exp(para.at(2));
    
    for(ii =0; ii<z.size(); ii++){
      dltrec_dmu = y.at(ii)/mu.at(ii) - (1 + phi * y.at(ii))/(1 + phi * mu.at(ii));
      
      if(z.at(ii) == 0){
        
        dmu_dkappa = (tau1.at(ii) + tau2.at(ii))*expXbeta.at(ii)*RHO.at(ii);
        dmu_deta = 0;
        dmu_dgamma = 0;
      }else if(z.at(ii) == 1){
        
        dmu_dkappa = (tau1.at(ii) + tau2.at(ii)*GAMMA)*expXbeta.at(ii)*RHO.at(ii);
        dmu_deta   = (1-RHO.at(ii))*expXbeta.at(ii);
        dmu_dgamma = tau2.at(ii)*RHO.at(ii)*KAPPA*expXbeta.at(ii);
        
      }else if(z.at(ii) == 2){
        
        dmu_dkappa = (tau1.at(ii)*GAMMA + tau2.at(ii))*expXbeta.at(ii)*RHO.at(ii);
        dmu_deta   = (1-RHO.at(ii))*expXbeta.at(ii);
        dmu_dgamma = tau1.at(ii)*RHO.at(ii)*KAPPA*expXbeta.at(ii);
        
      }else{
        dmu_dkappa = (tau1.at(ii)+tau2.at(ii))*GAMMA*expXbeta.at(ii)*RHO.at(ii);
        dmu_deta   = 2*(1-RHO.at(ii))*expXbeta.at(ii);
        dmu_dgamma = (tau1.at(ii)+tau2.at(ii))*KAPPA*expXbeta.at(ii)*RHO.at(ii);
      }
      
      grad.at(0) += dltrec_dmu*dmu_dkappa;
      grad.at(1) += dltrec_dmu*dmu_deta;
      grad.at(2) += dltrec_dmu*dmu_dgamma;
      
    }
    grad.at(0) *= KAPPA;
    grad.at(1) *= ETA;
    grad.at(2) *= GAMMA;
    
  }else if(H0 == 1){
    KAPPA = exp(para.at(0));
    ETA   = 1.0;
    GAMMA = exp(para.at(1));
    
    for(ii =0; ii<z.size(); ii++){
      dltrec_dmu = y.at(ii)/mu.at(ii) - (1 + phi * y.at(ii))/(1 + phi * mu.at(ii));
      
      if(z.at(ii) == 0){
        
        dmu_dkappa = (tau1.at(ii) + tau2.at(ii))*expXbeta.at(ii)*RHO.at(ii);
        dmu_dgamma = 0;
      }else if(z.at(ii) == 1){
        
        dmu_dkappa = (tau1.at(ii) + tau2.at(ii)*GAMMA)*expXbeta.at(ii)*RHO.at(ii);
        dmu_dgamma = tau2.at(ii)*RHO.at(ii)*KAPPA*expXbeta.at(ii);
        
      }else if(z.at(ii) == 2){
        
        dmu_dkappa = (tau1.at(ii)*GAMMA + tau2.at(ii))*expXbeta.at(ii)*RHO.at(ii);
        dmu_dgamma = tau1.at(ii)*RHO.at(ii)*KAPPA*expXbeta.at(ii);
        
      }else{
        dmu_dkappa = (tau1.at(ii)+tau2.at(ii))*GAMMA*expXbeta.at(ii)*RHO.at(ii);
        dmu_dgamma = (tau1.at(ii)+tau2.at(ii))*KAPPA*expXbeta.at(ii)*RHO.at(ii);
      }
      
      grad.at(0) += dltrec_dmu*dmu_dkappa;
      grad.at(1) += dltrec_dmu*dmu_dgamma;
      
    }
    grad.at(0) *= KAPPA;
    grad.at(1) *= GAMMA;
    
  }else{
    KAPPA = exp(para.at(0));
    ETA   = exp(para.at(1));
    GAMMA = 1.0;
    
    for(ii =0; ii<z.size(); ii++){
      dltrec_dmu = y.at(ii)/mu.at(ii) - (1 + phi * y.at(ii))/(1 + phi * mu.at(ii));
      
      if(z.at(ii) == 0){
        
        dmu_dkappa = (tau1.at(ii) + tau2.at(ii))*expXbeta.at(ii)*RHO.at(ii);
        dmu_deta = 0;
        
      }else if(z.at(ii) == 1){
        
        dmu_dkappa = (tau1.at(ii) + tau2.at(ii)*GAMMA)*expXbeta.at(ii)*RHO.at(ii);
        dmu_deta   = (1-RHO.at(ii))*expXbeta.at(ii);
        
      }else if(z.at(ii) == 2){
        
        dmu_dkappa = (tau1.at(ii)*GAMMA + tau2.at(ii))*expXbeta.at(ii)*RHO.at(ii);
        dmu_deta   = (1-RHO.at(ii))*expXbeta.at(ii);
        
      }else{
        dmu_dkappa = (tau1.at(ii)+tau2.at(ii))*GAMMA*expXbeta.at(ii)*RHO.at(ii);
        dmu_deta   = 2*(1-RHO.at(ii))*expXbeta.at(ii);
      }
      
      grad.at(0) += dltrec_dmu*dmu_dkappa;
      grad.at(1) += dltrec_dmu*dmu_deta;
      
    }
    grad.at(0) *= KAPPA;
    grad.at(1) *= ETA;
  }
  
  return grad;
}

// [[Rcpp::export]]
Rcpp::List RcppT_trec_KEG_BFGS(const arma::vec& para0, const int& H0, 
                               const arma::vec& y, const arma::vec& z,
                               const arma::vec& RHO, const arma::mat& X,
                               const arma::vec& BETA, const double& phi,
                               const arma::vec& tau1, const arma::vec& tau2,
                               const arma::vec& lgy1,
                               const arma::uword& max_iter = 4e3,
                               const double& eps = 1e-7,const bool& show = true){
  
  arma::uword num_params = para0.size();
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
  arma::vec offsets = arma::zeros<arma::vec>(y.n_elem);
  arma::vec expXbeta = arma::zeros<arma::vec>(y.n_elem);
  
  double old_LL,new_LL,inv_norm_p_k,tmp_alpha,ys;
  double fnscale = -1.0; // For maximization
  double curr_LL = 0.0;
  xk = para0;
  
  RcppT_compute_expXbeta(X, BETA, expXbeta); 
  
  while(iter < max_iter){
    //calculate direction p_k
    uu = 0;
    old_LL = fnscale * RcppT_loglikNB_KEG(xk, H0, y, z, phi, RHO, 
                                          tau1, tau2, lgy1, expXbeta,
                                          offsets, mu);
    gr_k   = fnscale * RcppT_grad_NB(xk, H0, y, z, RHO, phi,
                                     tau1, tau2, expXbeta, offsets, mu) ;
    p_k    = -1.0 * inv_Bk * gr_k;
    inv_norm_p_k =  1.0 / std::max(1.0, Rcpp_norm(p_k));
    
    //line search for new xk
    for(jj=0; jj<15; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4, jj);
      new_xk    = xk + tmp_alpha * p_k;
      new_LL    = fnscale * RcppT_loglikNB_KEG(new_xk, H0, y, z, phi, RHO, 
                                               tau1, tau2, lgy1, expXbeta,
                                               offsets, mu);
      if(new_LL < old_LL){ //minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * RcppT_grad_NB(new_xk, H0, y, z, RHO, phi,
                                      tau1, tau2, expXbeta, offsets, mu) - gr_k;
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
        gr_k = RcppT_grad_NB(xk, H0, y, z, RHO, phi,
                             tau1, tau2, expXbeta, offsets, mu);
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
  
  old_LL = Rcpp_loglikNB(y, phi, lgy1, mu);
  return Rcpp::List::create(
    Rcpp::Named("converge", converge),
    Rcpp::Named("LL", old_LL),
    Rcpp::Named("iter", iter),
    Rcpp::Named("norm_GRAD", Rcpp_norm(gr_k)),
    Rcpp::Named("PAR", Rcpp::NumericVector(xk.begin(), xk.end()))
  );
}

// [[Rcpp::export]]
Rcpp::List RcppT_trec_sfit(const int& H0, const arma::vec& para0,
                           const arma::vec& y, 
                           const arma::vec& z, const arma::vec& RHO, 
                           const arma::mat& X,const arma::vec& tau1, 
                           const arma::vec& tau2, const arma::vec& lgy1, 
                           const arma::uword& max_iter = 4e3,
                           const double& eps = 1e-7, const bool& show = false){
  arma::uword iter = 0;
  arma::uword converge = 0;
  arma::uword pp = X.n_cols + 1;
  arma::vec offsets = arma::zeros<arma::vec>(z.n_elem);
  Rcpp::List new_reg, new_keg_fit;
  double new_LL, curr_LL;
  arma::vec curr_para = para0;
  arma::vec new_para = arma::zeros<arma::vec>(para0.size());
  arma::vec curr_reg_par = arma::zeros<arma::vec>(pp);
  arma::vec new_reg_par = arma::zeros<arma::vec>(pp);
  arma::vec BETA = arma::zeros<arma::vec>(pp-1);
  double KAPPA, ETA, GAMMA;
  double phi = 0.0;
  // initial regression fit
  
  while(iter < max_iter){
    //update BETA, phi
    
    if(H0 == 0){
      KAPPA = exp(new_para.at(0));
      ETA   = exp(new_para.at(1));
      GAMMA = exp(new_para.at(2));
    }else if(H0 == 1){ //para = log(c(KAPPA, GAMMA))
      KAPPA = exp(new_para.at(0));
      ETA   = 1.0;
      GAMMA = exp(new_para.at(1));
    }else{ //para = log(c(KAPPA, ETA))
      KAPPA = exp(new_para.at(0));
      ETA   = exp(new_para.at(1));
      GAMMA = 1.0;
    }
    
    RcppT_compute_offset(z, RHO, KAPPA, ETA, GAMMA, tau1, tau2, offsets);
    new_reg     = RcppT_reg_BFGS(y, X, offsets, curr_reg_par, lgy1,
                                 max_iter, eps, false);
    new_LL      = as<double>(new_reg["LL"]);
    new_reg_par = as<arma::vec>(new_reg["PAR"]);
    
    BETA = new_reg_par.subvec(0, X.n_cols-1);
    phi  = std::exp(new_reg_par.at(pp-1));
    
    if(show){
      Rprintf("TReC: BETA, PHI updated after %d iter \n",
              as<int>(new_reg["iter"]));
    }
    
    //update KEG
    new_keg_fit = RcppT_trec_KEG_BFGS(curr_para, H0, y, z, RHO, X, BETA, phi,
                                      tau1, tau2, lgy1, max_iter, eps, false);
    new_para    = Rcpp::as<arma::vec>(new_keg_fit["PAR"]);
    new_LL      = as<double>(new_keg_fit["LL"]);
    
    if(show){
      Rprintf("TReC: keg updated after %d iter \n",
              as<int>(new_keg_fit["iter"]));
      // printR_obj(new_LL);
      // printR_obj(curr_LL);
    }
    
    if(iter > 0){
      if(abs(curr_LL - new_LL) < eps && 
         Rcpp_norm(curr_reg_par - new_reg_par) < eps &&
         Rcpp_norm(curr_para - new_para) < eps){
        // convergence criteria 
        if((curr_LL - new_LL > 0.0 && as<double>(new_keg_fit["norm_GRAD"]) > 0.01)){
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
    curr_para     = new_para;
    curr_LL      = new_LL;
    iter++;
    
  }
  if(show){
    Rprintf("TReC converges after %d iter \n", iter);
  }
  return Rcpp::List::create(
    Rcpp::Named("PAR", Rcpp::NumericVector(new_para.begin(), new_para.end())),
    Rcpp::Named("reg_par", Rcpp::NumericVector(new_reg_par.begin(), new_reg_par.end())),
    Rcpp::Named("LL", new_LL),
    Rcpp::Named("converge", converge),
    Rcpp::Named("iter", iter)
  );
}

// [[Rcpp::export]]
Rcpp::List RcppT_trec(const arma::vec& y, const arma::vec& z,
                      const arma::vec& RHO,const arma::mat& X, 
                      const arma::vec& tau1, const arma::vec& tau2,
                      const arma::vec& lgy1, 
                      const arma::uword& max_iter = 4e3,
                      const double& eps = 1e-7, const bool& show = false){
  Rcpp::List sfit0, sfit1, sfit2, sfit3;
  double p_eta, p_gamma;
  arma::vec para0 = arma::zeros<arma::vec>(3);
  int converge = 0;
  
  sfit0 = RcppT_trec_sfit(0, para0, y, z, RHO, X, tau1, tau2, lgy1, 
                          max_iter, eps, show); 
  sfit1 = RcppT_trec_sfit(1, para0.subvec(0,1), y, z, RHO, X, tau1, tau2, lgy1, 
                          max_iter, eps, show); 
  sfit2 = RcppT_trec_sfit(2, para0.subvec(0,1), y, z, RHO, X, tau1, tau2, lgy1, 
                          max_iter, eps, show); 
  
  para0.subvec(0,1) = as<arma::vec>(sfit1["PAR"]);
  para0.at(2) = as<arma::vec>(sfit1["PAR"]).at(2); 
  sfit3 = RcppT_trec_sfit(0, para0, y, z, RHO, X, tau1, tau2, lgy1, 
                          max_iter, eps, show); 
  
  if(as<double>(sfit0["LL"]) < as<double>(sfit3["LL"]) && sfit3["converge"] ){
    sfit0 = sfit3;  
  }
  
  converge = sfit0["converge"] && sfit1["converge"] && sfit2["converge"];
  
  Rcpp::NumericVector PAR = as<Rcpp::NumericVector>(sfit0["PAR"]);
  Rcpp::NumericVector reg_par = as<Rcpp::NumericVector>(sfit0["reg_par"]);
  
  p_eta = R::pchisq(2*(as<double>(sfit0["LL"])-as<double>(sfit1["LL"])),
                    1, 0, 0);
  p_gamma = R::pchisq(2*(as<double>(sfit0["LL"])-as<double>(sfit2["LL"])),
                      1, 0, 0);
  
  return Rcpp::List::create(
    Rcpp::Named("p_eta", p_eta),
    Rcpp::Named("p_gamma", p_gamma),
    Rcpp::Named("PAR", Rcpp::NumericVector(PAR.begin(), PAR.end())),
    Rcpp::Named("reg_par", Rcpp::NumericVector(reg_par.begin(), reg_par.end())),
    Rcpp::Named("LL", sfit0["LL"]),
    Rcpp::Named("LL_eta", sfit1["LL"]),
    Rcpp::Named("LL_gamma", sfit2["LL"]),
    Rcpp::Named("converge", converge)
  );
  
}


// ----------------------------------------------------------------------
// ASE
// ----------------------------------------------------------------------

// [[Rcpp::export]]
void RcppT_compite_pi(const arma::vec& z_AS, const arma::vec& RHO_AS, 
                      const double& KAPPA, const double& ETA, const double& GAMMA,
                      const arma::vec& tauB, const arma::vec& tau, 
                      arma::vec& pis){
  arma::uword ii;
  double tmp1, tmp2; 
  
  for(ii=0;ii<z_AS.n_elem;ii++){
    if(z_AS.at(ii) == 0 | z_AS.at(ii) == 3){
      tmp1 = RHO_AS.at(ii)*tauB.at(ii)*KAPPA + 1 - RHO_AS.at(ii);
      tmp2 = RHO_AS.at(ii)*tau.at(ii)*KAPPA + 2*(1 - RHO_AS.at(ii));
      pis.at(ii) = tmp1/tmp2;
      
    }else{
      tmp1 = RHO_AS.at(ii)*tauB.at(ii)*KAPPA*GAMMA + (1 - RHO_AS.at(ii))*ETA;
      tmp2 = RHO_AS.at(ii)*(tau.at(ii) - tauB.at(ii))*KAPPA + 
        (1 - RHO_AS.at(ii)) + tmp1;
      pis.at(ii) = tmp1/tmp2;
    }
  }
}

// [[Rcpp::export]]
double RcppT_loglikBB_KEG(const arma::vec& para, const int& H0,
                          const arma::vec& z_AS, const arma::vec& RHO_AS, 
                          const arma::vec& ni, const arma::vec& ni0,
                          const double& log_theta, const arma::vec& lbc, 
                          const arma::vec& tauB, const arma::vec& tau, 
                          arma::vec& pis){
  //update KEG
  arma::uword ii;
  double loglik = 0.0;
  double vtheta = std::exp(-log_theta);
  double lgvt = lgamma(vtheta);
  double KAPPA, ETA, GAMMA;
  
  if(H0 == 0){
    KAPPA = exp(para.at(0));
    ETA   = exp(para.at(1));
    GAMMA = exp(para.at(2));
  }else if(H0 == 1){ //para = log(c(KAPPA, GAMMA))
    KAPPA = exp(para.at(0));
    ETA   = 1.0;
    GAMMA = exp(para.at(1));
  }else{ //para = log(c(KAPPA, ETA))
    KAPPA = exp(para.at(0));
    ETA   = exp(para.at(1));
    GAMMA = 1.0;
  }
  
  RcppT_compite_pi(z_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau, pis);
  
  for(ii=0;ii<ni.n_elem;ii++){
    double aa = pis.at(ii) * vtheta; 
    double bb = vtheta - aa;
    double lgvab = - lgamma(aa) - lgamma(bb) + lgvt;
    
    //printR_obj(loglik);
    loglik += lbc.at(ii) + lgamma(aa + ni0.at(ii)) +
      lgamma(bb + ni.at(ii) - ni0.at(ii)) + lgvab ;
    //- lgamma(vtheta + ni.at(ii));
    
  }
  
  return loglik;
}

// [[Rcpp::export]]
arma::vec RcppT_grad_BB_KEG(const arma::vec& para, const int& H0,
                            const arma::vec& z_AS, const arma::vec& RHO_AS, 
                            const arma::vec& ni, const arma::vec& ni0,
                            const double& log_theta, const arma::vec& lbc, 
                            const arma::vec& tauB, const arma::vec& tau, 
                            const arma::vec& pis){
  arma::uword ii;
  arma::vec grad = arma::zeros<arma::vec>(para.size());
  double KAPPA, ETA, GAMMA;
  double dlASE_dpi, dpi_dkappa, dpi_deta, dpi_dgamma, tmp1, tmp2, tmp3, tmp4;
  double vtheta = std::exp(-log_theta);
  double lgvt = lgamma(vtheta);
  
  if(H0 == 0){
    KAPPA = exp(para.at(0));
    ETA   = exp(para.at(1));
    GAMMA = exp(para.at(2));
    
    for(ii =0; ii<z_AS.size(); ii++){
      double aa = pis.at(ii) * vtheta; 
      double bb = vtheta - aa;
      double diaa = R::digamma(aa);
      double dibb =  R::digamma(bb);
      double diaa_ni0 = R::digamma(aa + ni0.at(ii));
      double dibb_ni1 = R::digamma(bb + ni.at(ii) - ni0.at(ii));
      
      dlASE_dpi = diaa_ni0 - dibb_ni1 - diaa + dibb;
      
      if(z_AS.at(ii) == 0 | z_AS.at(ii) == 3){
        tmp1 = RHO_AS.at(ii)*tau.at(ii)*KAPPA + 2*(1 - RHO_AS.at(ii));
        tmp2 = RHO_AS.at(ii)*tauB.at(ii)/tmp1;
        tmp3 = RHO_AS.at(ii)*tau.at(ii)*(RHO_AS.at(ii)*tauB.at(ii)*KAPPA + 
          (1 - RHO_AS.at(ii)))/pow(tmp1, 2.0);
        dpi_dkappa = tmp2 - tmp3;
        dpi_deta   = 0;
        dpi_dgamma = 0;
        
      }else{
        tmp1 = RHO_AS.at(ii)*tauB.at(ii)*GAMMA*KAPPA + (1 - RHO_AS.at(ii))*ETA;
        tmp2 = RHO_AS.at(ii)*(tau.at(ii) - tauB.at(ii))*KAPPA +
          (1 - RHO_AS.at(ii)) + tmp1;
        tmp3 = tmp1/pow(tmp2, 2.0);
        tmp4 = 1/tmp2 -tmp3; 
        dpi_dkappa = RHO_AS.at(ii)*tauB.at(ii)*GAMMA*tmp4 -
          RHO_AS.at(ii)*(tau.at(ii) - tauB.at(ii))*tmp3;
        dpi_deta   = (1 - RHO_AS.at(ii))*tmp4;
        dpi_dgamma = RHO_AS.at(ii)*tauB.at(ii)*KAPPA*tmp4;
      }
      
      grad.at(0) += dlASE_dpi*dpi_dkappa;
      grad.at(1) += dlASE_dpi*dpi_deta;
      grad.at(2) += dlASE_dpi*dpi_dgamma;
      
    }
    grad.at(0) *= KAPPA * vtheta;
    grad.at(1) *= ETA * vtheta;
    grad.at(2) *= GAMMA* vtheta;
    
    
  }else if(H0 == 1){ //para = log(c(KAPPA, GAMMA))
    KAPPA = exp(para.at(0));
    ETA   = 1.0;
    GAMMA = exp(para.at(1));
    for(ii =0; ii<z_AS.size(); ii++){
      double aa = pis.at(ii) * vtheta; 
      double bb = vtheta - aa;
      double diaa = R::digamma(aa);
      double dibb =  R::digamma(bb);
      double diaa_ni0 = R::digamma(aa + ni0.at(ii));
      double dibb_ni1 = R::digamma(bb + ni.at(ii) - ni0.at(ii));
      
      dlASE_dpi = diaa_ni0 - dibb_ni1 - diaa + dibb;
      
      if(z_AS.at(ii) == 0 | z_AS.at(ii) == 3){
        tmp1 = RHO_AS.at(ii)*tau.at(ii)*KAPPA + 2*(1 - RHO_AS.at(ii));
        tmp2 = RHO_AS.at(ii)*tauB.at(ii)/tmp1;
        tmp3 = RHO_AS.at(ii)*tau.at(ii)*(RHO_AS.at(ii)*tauB.at(ii)*KAPPA + 
          (1 - RHO_AS.at(ii)))/pow(tmp1, 2.0);
        dpi_dkappa = tmp2 - tmp3;
        dpi_dgamma = 0;
        
      }else{
        tmp1 = RHO_AS.at(ii)*tauB.at(ii)*GAMMA*KAPPA + (1 - RHO_AS.at(ii))*ETA;
        tmp2 = RHO_AS.at(ii)*(tau.at(ii) - tauB.at(ii))*KAPPA +
          (1 - RHO_AS.at(ii)) + tmp1;
        tmp3 = tmp1/pow(tmp2, 2.0);
        tmp4 = 1/tmp2 -tmp3; 
        dpi_dkappa = RHO_AS.at(ii)*tauB.at(ii)*GAMMA*tmp4 -
          RHO_AS.at(ii)*(tau.at(ii) - tauB.at(ii))*tmp3;
        dpi_dgamma = RHO_AS.at(ii)*tauB.at(ii)*KAPPA*tmp4;
      }
      
      grad.at(0) += dlASE_dpi*dpi_dkappa;
      grad.at(1) += dlASE_dpi*dpi_dgamma;
      
    }
    grad.at(0) *= KAPPA * vtheta;
    grad.at(1) *= GAMMA* vtheta;
    
    
  }else{ //para = log(c(KAPPA, ETA))
    KAPPA = exp(para.at(0));
    ETA   = exp(para.at(1));
    GAMMA = 1.0;
    for(ii =0; ii<z_AS.size(); ii++){
      double aa = pis.at(ii) * vtheta; 
      double bb = vtheta - aa;
      double diaa = R::digamma(aa);
      double dibb =  R::digamma(bb);
      double diaa_ni0 = R::digamma(aa + ni0.at(ii));
      double dibb_ni1 = R::digamma(bb + ni.at(ii) - ni0.at(ii));
      
      dlASE_dpi = diaa_ni0 - dibb_ni1 - diaa + dibb;
      
      if(z_AS.at(ii) == 0 | z_AS.at(ii) == 3){
        tmp1 = RHO_AS.at(ii)*tau.at(ii)*KAPPA + 2*(1 - RHO_AS.at(ii));
        tmp2 = RHO_AS.at(ii)*tauB.at(ii)/tmp1;
        tmp3 = RHO_AS.at(ii)*tau.at(ii)*(RHO_AS.at(ii)*tauB.at(ii)*KAPPA + 
          (1 - RHO_AS.at(ii)))/pow(tmp1, 2.0);
        dpi_dkappa = tmp2 - tmp3;
        dpi_deta   = 0;
        
      }else{
        tmp1 = RHO_AS.at(ii)*tauB.at(ii)*GAMMA*KAPPA + (1 - RHO_AS.at(ii))*ETA;
        tmp2 = RHO_AS.at(ii)*(tau.at(ii) - tauB.at(ii))*KAPPA +
          (1 - RHO_AS.at(ii)) + tmp1;
        tmp3 = tmp1/pow(tmp2, 2.0);
        tmp4 = 1/tmp2 -tmp3; 
        dpi_dkappa = RHO_AS.at(ii)*tauB.at(ii)*GAMMA*tmp4 -
          RHO_AS.at(ii)*(tau.at(ii) - tauB.at(ii))*tmp3;
        dpi_deta   = (1 - RHO_AS.at(ii))*tmp4;
      }
      
      grad.at(0) += dlASE_dpi*dpi_dkappa;
      grad.at(1) += dlASE_dpi*dpi_deta;
      
    }
    grad.at(0) *= KAPPA * vtheta;
    grad.at(1) *= ETA * vtheta;
    
  }
  
  return grad;
  
}

// [[Rcpp::export]]
double RcppT_loglikBB_THETA(const arma::vec& ni, const arma::vec& ni0, 
                            const double& log_theta, const arma::vec& pis,
                            const arma::vec& lbc){
  
  arma::uword ii;
  double loglik = 0.0;
  double vtheta = std::exp(-log_theta);
  
  for(ii=0;ii<ni.n_elem;ii++){
    
    double aa = pis.at(ii)*vtheta;
    double bb = vtheta - aa;
    double lgvt = lgamma(vtheta);
    double lgvab = - lgamma(aa) - lgamma(bb) + lgvt;
    
    loglik += lbc.at(ii) + lgamma(aa + ni0.at(ii)) +
      lgamma(bb + ni.at(ii) - ni0.at(ii)) + lgvab -
      lgamma(vtheta + ni.at(ii));
    
  }
  
  return loglik;
}

// [[Rcpp::export]]
double RcppT_grad_BB_THETA(const arma::vec& ni, const arma::vec& ni0,
                           const double& log_theta, const arma::vec& pis){
  // only update theta
  arma::uword ii;
  double grad   = 0.0;
  double vtheta = std::exp(-log_theta);
  double divt   = R::digamma(vtheta);
  double diaa_ni0, dibb_ni1, aa, bb, diaa, dibb;
  
  for(ii=0;ii<ni.n_elem;ii++){
    aa      = pis.at(ii)*vtheta;
    bb      = vtheta - aa;
    diaa    = R::digamma(aa);
    dibb    = R::digamma(bb);
    diaa_ni0 = R::digamma(aa + ni0.at(ii));
    dibb_ni1 = R::digamma(bb + ni.at(ii) - ni0.at(ii));
    
    //dlase_dtheta
    grad += - pis.at(ii)*(diaa_ni0 - diaa) -
      (1.0 - pis.at(ii)) * (dibb_ni1 - dibb) -
      (divt - R::digamma(vtheta + ni.at(ii)));
    
  }
  grad *= vtheta;
  return grad; 
}


// [[Rcpp::export]]
Rcpp::List RcppT_ase_theta_BFGS(const double& lg_theta,
                                const arma::vec& ni, const arma::vec& ni0,
                                const arma::vec& pis, const arma::vec& lbc,
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
    
    old_LL = fnscale * RcppT_loglikBB_THETA(ni, ni0, xk.at(0), pis, lbc);
    
    //if(old_LL < 0) break;
    
    //printR_obj(old_LL);
    gr_k   = fnscale * RcppT_grad_BB_THETA(ni, ni0, xk.at(0), pis);
    p_k    = -1.0 * inv_Bk * gr_k;
    inv_norm_p_k =  1.0 / std::max(1.0, Rcpp_norm(p_k));
    
    //line search for new xk
    for(jj=0; jj<15; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4, jj);
      new_xk    = xk + tmp_alpha * p_k;
      new_LL    = fnscale * RcppT_loglikBB_THETA(ni, ni0,
                                                 new_xk.at(0), pis, lbc);
      
      if(new_LL < old_LL){ //minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * RcppT_grad_BB_THETA(ni, ni0,
                                            new_xk.at(0), pis) - gr_k;
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
        gr_k = RcppT_grad_BB_THETA(ni, ni0, xk.at(0), pis);
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
  
  old_LL = RcppT_loglikBB_THETA(ni, ni0, xk.at(0), pis, lbc);
  
  return Rcpp::List::create(
    Rcpp::Named("converge", converge),
    Rcpp::Named("LL", old_LL),
    Rcpp::Named("iter", iter),
    Rcpp::Named("norm_GRAD", Rcpp_norm(gr_k)),
    Rcpp::Named("PAR", Rcpp::NumericVector(xk.begin(), xk.end()))
  );
}

// [[Rcpp::export]]
Rcpp::List RcppT_ase_KEG_BFGS(const arma::vec& para0, const int& H0, 
                              const arma::vec& z_AS, const arma::vec& RHO_AS,
                              const arma::vec& ni0, const arma::vec& ni, 
                              const double& log_theta, const arma::vec& tauB,
                              const arma::vec& tau, const arma::vec& lbc,
                              const arma::uword& max_iter = 4e3,
                              const double& eps = 1e-7, const bool& show = true){
  arma::uword num_params = para0.size();
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
  
  arma::vec pis = arma::zeros<arma::vec>(z_AS.n_elem);
  
  double old_LL,new_LL,inv_norm_p_k,tmp_alpha,ys;
  double fnscale = -1.0; // For maximization
  double curr_LL = 0.0;
  xk = para0;
  
  while(iter < max_iter){
    //calculate direction p_k
    uu = 0;
    old_LL = fnscale * RcppT_loglikBB_KEG(xk, H0, z_AS, RHO_AS, ni, ni0, 
                                          log_theta, lbc, tauB, 
                                          tau, pis);
    gr_k   = fnscale * RcppT_grad_BB_KEG(xk, H0, z_AS, RHO_AS, ni, ni0,
                                         log_theta, lbc, tauB, 
                                         tau, pis);
    p_k    = -1.0 * inv_Bk * gr_k;
    inv_norm_p_k =  1.0 / std::max(1.0, Rcpp_norm(p_k));
    
    //line search for new xk
    for(jj=0; jj<15; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4, jj);
      new_xk    = xk + tmp_alpha * p_k;
      new_LL    = fnscale * RcppT_loglikBB_KEG(new_xk, H0, z_AS, RHO_AS, ni, ni0, 
                                               log_theta, lbc, tauB, 
                                               tau, pis) ;
      if(new_LL < old_LL){ //minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * RcppT_grad_BB_KEG(new_xk, H0, z_AS, RHO_AS, ni, ni0, 
                                          log_theta, lbc, tauB, 
                                          tau, pis) - gr_k;
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
        gr_k = RcppT_grad_BB_KEG(xk, H0, z_AS, RHO_AS, ni, ni0,
                                 log_theta, lbc, tauB, 
                                 tau, pis);
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
  
  old_LL = RcppT_loglikBB_THETA(ni, ni0, log_theta, pis, lbc);
  return Rcpp::List::create(
    Rcpp::Named("converge", converge),
    Rcpp::Named("LL", old_LL),
    Rcpp::Named("iter", iter),
    Rcpp::Named("norm_GRAD", Rcpp_norm(gr_k)),
    Rcpp::Named("PAR", Rcpp::NumericVector(xk.begin(), xk.end()))
  );
}


// [[Rcpp::export]]
Rcpp::List RcppT_ase_sfit(const int& H0, const arma::vec& para0,
                          const arma::vec& z_AS, const arma::vec& RHO_AS,
                          const arma::vec& ni0, const arma::vec& ni,
                          const arma::vec& tauB,const arma::vec& tau, 
                          const arma::vec& lbc,const arma::uword& max_iter = 4e3,
                          const double& eps = 1e-7, const bool& show = false){
  arma::uword iter = 0;
  arma::uword converge = 0;
  Rcpp::List new_reg, new_keg_fit, new_theta_fit;
  arma::vec curr_para = para0;
  arma::vec new_para = arma::zeros<arma::vec>(para0.size());
  arma::vec pis = arma::zeros<arma::vec>(z_AS.n_elem);
  double curr_LL, new_LL, KAPPA, ETA, GAMMA, new_theta;
  double curr_theta = 0.1;
  // initial regression fit
  
  while(iter < max_iter){
    
    if(H0 == 0){
      KAPPA = exp(new_para.at(0));
      ETA   = exp(new_para.at(1));
      GAMMA = exp(new_para.at(2));
    }else if(H0 == 1){ //para = log(c(KAPPA, GAMMA))
      KAPPA = exp(new_para.at(0));
      ETA   = 1.0;
      GAMMA = exp(new_para.at(1));
    }else{ //para = log(c(KAPPA, ETA))
      KAPPA = exp(new_para.at(0));
      ETA   = exp(new_para.at(1));
      GAMMA = 1.0;
    }
    //update THETA
    RcppT_compite_pi(z_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau, pis);
    
    new_theta_fit = RcppT_ase_theta_BFGS(curr_theta, ni, ni0, pis, lbc);
    new_theta     =  as<double>(new_theta_fit["PAR"]);
    
    if(show){
      Rprintf("ASE: Theta updated after %d iter \n",
              as<int>(new_keg_fit["iter"]));
      // printR_obj(new_LL);
      // printR_obj(curr_LL);
    }
    //update KEG
    new_keg_fit = RcppT_ase_KEG_BFGS(curr_para, H0, z_AS, RHO_AS, 
                                     ni0, ni, new_theta, tauB, 
                                     tau, lbc, max_iter, eps, false);
    new_para    = Rcpp::as<arma::vec>(new_keg_fit["PAR"]);
    new_LL      = as<double>(new_keg_fit["LL"]);
    
    if(show){
      Rprintf("ASE: keg updated after %d iter \n",
              as<int>(new_keg_fit["iter"]));
      // printR_obj(new_LL);
      // printR_obj(curr_LL);
    }
    
    if(iter > 0){
      if(abs(curr_LL - new_LL) < eps && 
         Rcpp_norm(curr_para - new_para) < eps){
        // convergence criteria 
        if((curr_LL - new_LL > 0.0 && as<double>(new_keg_fit["norm_GRAD"]) > 0.01)){
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
    curr_para    = new_para;
    curr_LL      = new_LL;
    curr_theta   = new_theta;
    iter++;
    
  }
  if(show){
    Rprintf("ASE converges after %d iter \n", iter);
  }
  return Rcpp::List::create(
    Rcpp::Named("PAR", Rcpp::NumericVector(new_para.begin(), new_para.end())),
    Rcpp::Named("log_theta", new_theta),
    Rcpp::Named("LL", new_LL),
    Rcpp::Named("converge", converge),
    Rcpp::Named("iter", iter)
  );
}

// [[Rcpp::export]]
Rcpp::List RcppT_ase(const arma::vec& z_AS,const arma::vec& RHO_AS, 
                     const arma::vec& ni0, const arma::vec& ni, 
                     const arma::vec& tauB, const arma::vec& tau, 
                     const arma::vec& lbc, 
                     const arma::uword& max_iter = 4e3,
                     const double& eps = 1e-7, const bool& show = false){
  Rcpp::List sfit0, sfit1, sfit2, sfit3;
  double p_eta, p_gamma;
  arma::vec para0 = arma::zeros<arma::vec>(3);
  int converge = 0;
  
  sfit0 = RcppT_ase_sfit(0, para0, z_AS, RHO_AS,  ni0, ni, tauB, tau, lbc, 
                         max_iter, eps, show); 
  sfit1 = RcppT_ase_sfit(1, para0.subvec(0,1), z_AS, RHO_AS,  ni0, ni, 
                         tauB, tau, lbc, max_iter, eps, show); 
  sfit2 = RcppT_ase_sfit(2, para0.subvec(0,1),  z_AS, RHO_AS,  ni0, ni, 
                         tauB, tau, lbc, max_iter, eps, show); 
  
  para0.subvec(0,1) = as<arma::vec>(sfit1["PAR"]);
  para0.at(2) = as<arma::vec>(sfit1["PAR"]).at(2); 
  sfit3 = RcppT_ase_sfit(0, para0, z_AS, RHO_AS,  ni0, ni, 
                         tauB, tau, lbc, max_iter, eps, show); 
  
  if(as<double>(sfit0["LL"]) < as<double>(sfit3["LL"]) && sfit3["converge"] ){
    sfit0 = sfit3;  
  }
  
  converge = sfit0["converge"] && sfit1["converge"] && sfit2["converge"];
  
  Rcpp::NumericVector PAR = as<Rcpp::NumericVector>(sfit0["PAR"]);
  
  p_eta = R::pchisq(2*(as<double>(sfit0["LL"])-as<double>(sfit1["LL"])),
                    1, 0, 0);
  p_gamma = R::pchisq(2*(as<double>(sfit0["LL"])-as<double>(sfit2["LL"])),
                      1, 0, 0);
  
  return Rcpp::List::create(
    Rcpp::Named("p_eta", p_eta),
    Rcpp::Named("p_gamma", p_gamma),
    Rcpp::Named("PAR", Rcpp::NumericVector(PAR.begin(), PAR.end())),
    Rcpp::Named("log_theta", sfit0["log_theta"]),
    Rcpp::Named("LL", sfit0["LL"]),
    Rcpp::Named("LL_eta", sfit1["LL"]),
    Rcpp::Named("LL_gamma", sfit2["LL"]),
    Rcpp::Named("converge", converge)
  );
  
}


// ----------------------------------------------------------------------
// TRECASE
// ----------------------------------------------------------------------

// [[Rcpp::export]]
double RcppT_TReCASE_LL_KEG(const arma::vec& para, const int& H0,
                            const arma::vec& y, const arma::vec& z,
                            const arma::vec& z_AS, const double& phi,
                            const arma::vec& RHO, const arma::vec& RHO_AS,
                            const arma::vec& tau1,const arma::vec& tau2,
                            const arma::vec& ni0, const arma::vec& ni,
                            const double& log_theta, const arma::vec& tauB,
                            const arma::vec& tau, const arma::vec& lgy1,
                            const arma::vec& expXbeta,const arma::vec& lbc,
                            arma::vec& offsets, 
                            arma::vec& pis, arma::vec& mu){
  double LL_NB, LL_BB;
  LL_NB = RcppT_loglikNB_KEG(para, H0, y, z, phi, RHO, tau1,
                             tau2, lgy1, expXbeta, offsets, mu);
  LL_BB = RcppT_loglikBB_KEG(para, H0, z_AS, RHO_AS, ni, ni0, log_theta,
                             lbc, tauB, tau, pis);
  return LL_NB + LL_BB;
}

// [[Rcpp::export]]
double RcppT_TReCASE_LL(const arma::vec& y, const double& phi, 
                        const arma::vec& lgy1, const arma::vec& mu, 
                        const arma::vec& ni0, const arma::vec& ni,
                        const double& log_theta, const arma::vec& pis,
                        const arma::vec& lbc){
  double LL_NB, LL_BB;
  LL_NB = Rcpp_loglikNB(y, phi, lgy1, mu);
  LL_BB = RcppT_loglikBB_THETA(ni, ni0, log_theta, pis, lbc);
  return LL_NB + LL_BB;
  
}

// [[Rcpp::export]]
arma::vec RcppT_TReCASE_grad_KEG(const arma::vec& para, const int& H0,
                                 const arma::vec& y, const arma::vec& z,
                                 const arma::vec& z_AS, const double& phi,
                                 const arma::vec& RHO, const arma::vec& RHO_AS,
                                 const arma::vec& tau1,const arma::vec& tau2,
                                 const arma::vec& ni0, const arma::vec& ni,
                                 const double& log_theta, const arma::vec& tauB,
                                 const arma::vec& tau, const arma::vec& lgy1,
                                 const arma::vec& expXbeta,const arma::vec& lbc,
                                 const arma::vec& offsets, 
                                 const arma::vec& pis, const arma::vec& mu){
  arma::vec grad_NB = arma::zeros<arma::vec>(para.n_elem);
  arma::vec grad_BB = arma::zeros<arma::vec>(para.n_elem);
  grad_NB = RcppT_grad_NB(para, H0, y, z, RHO, phi, tau1, tau2, expXbeta, 
                          offsets, mu);
  
  grad_BB = RcppT_grad_BB_KEG(para, H0, z_AS, RHO_AS, ni, ni0, log_theta, 
                              lbc, tauB, tau, pis);
  return grad_NB + grad_BB;
  
}

// [[Rcpp::export]]
Rcpp::List RcppT_trecase_KEG_BFGS(const arma::vec& para0, const int& H0, 
                                  const arma::vec& y, const arma::vec& z, 
                                  const arma::vec& z_AS,
                                  const arma::vec& RHO, const arma::vec& RHO_AS,
                                  const arma::mat& X, const arma::vec& BETA, 
                                  const double& phi,
                                  const arma::vec& tau1, const arma::vec& tau2,
                                  const arma::vec& lgy1,  
                                  const arma::vec& ni0, const arma::vec& ni,
                                  const double& log_theta, const arma::vec& tauB,
                                  const arma::vec& tau, const arma::vec& lbc,
                                  const arma::uword& max_iter = 4e3,
                                  const double& eps = 1e-7,const bool& show = true){
  
  arma::uword num_params = para0.size();
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
  arma::vec offsets = arma::zeros<arma::vec>(y.n_elem);
  arma::vec expXbeta = arma::zeros<arma::vec>(y.n_elem);
  arma::vec pis = arma::zeros<arma::vec>(z_AS.n_elem);
  
  double old_LL,new_LL,inv_norm_p_k,tmp_alpha,ys;
  double fnscale = -1.0; // For maximization
  double curr_LL = 0.0;
  xk = para0;
  
  RcppT_compute_expXbeta(X, BETA, expXbeta); 
  
  while(iter < max_iter){
    //calculate direction p_k
    uu = 0;
    old_LL = fnscale * RcppT_TReCASE_LL_KEG(xk, H0, y, z, z_AS, phi, RHO, RHO_AS, 
                                            tau1, tau2, ni0, ni, log_theta, tauB, 
                                            tau, lgy1, expXbeta, lbc, offsets, 
                                            pis, mu);
    gr_k   = fnscale * RcppT_TReCASE_grad_KEG(xk, H0, y, z, z_AS, phi, RHO, 
                                              RHO_AS, tau1, tau2, ni0, ni, 
                                              log_theta, tauB, tau, lgy1, expXbeta, 
                                              lbc, offsets, pis, mu);
    p_k    = -1.0 * inv_Bk * gr_k;
    inv_norm_p_k =  1.0 / std::max(1.0, Rcpp_norm(p_k));
    
    //line search for new xk
    for(jj=0; jj<15; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4, jj);
      new_xk    = xk + tmp_alpha * p_k;
      new_LL    = fnscale * RcppT_TReCASE_LL_KEG(new_xk, H0, y, z, z_AS, phi,
                                                 RHO, RHO_AS, tau1, tau2, 
                                                 ni0, ni, log_theta, tauB, 
                                                 tau, lgy1, expXbeta, lbc,
                                                 offsets, pis, mu);
      if(new_LL < old_LL){ //minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * RcppT_TReCASE_grad_KEG(new_xk, H0, y, z, z_AS, phi, RHO, 
                                               RHO_AS, tau1, tau2, ni0, ni, 
                                               log_theta, tauB, tau, lgy1,
                                               expXbeta, lbc, offsets, pis,
                                               mu) - gr_k;
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
        gr_k = RcppT_TReCASE_grad_KEG(xk, H0, y, z, z_AS, phi, RHO, RHO_AS, 
                                      tau1, tau2, ni0, ni, log_theta, tauB, tau, 
                                      lgy1, expXbeta, lbc, offsets, pis, mu);
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
  
  old_LL = RcppT_TReCASE_LL(y, phi, lgy1, mu, ni0, ni,log_theta, pis, lbc);
  return Rcpp::List::create(
    Rcpp::Named("converge", converge),
    Rcpp::Named("LL", old_LL),
    Rcpp::Named("iter", iter),
    Rcpp::Named("norm_GRAD", Rcpp_norm(gr_k)),
    Rcpp::Named("PAR", Rcpp::NumericVector(xk.begin(), xk.end()))
  );
}


// [[Rcpp::export]]
Rcpp::List RcppT_trecase_sfit(const int& H0, const arma::vec& para0,
                              const arma::vec& y, 
                              const arma::vec& z, const arma::vec& z_AS,
                              const arma::vec& RHO, const arma::vec& RHO_AS,
                              const arma::mat& X, 
                              const arma::vec& tau1, const arma::vec& tau2,
                              const arma::vec& lgy1,  
                              const arma::vec& ni0, const arma::vec& ni,
                              const arma::vec& tauB,const arma::vec& tau, 
                              const arma::vec& lbc,const arma::uword& max_iter = 4e3,
                              const double& eps = 1e-7, const bool& show = false){
  arma::uword iter = 0;
  arma::uword converge = 0;
  arma::uword pp = X.n_cols + 1;
  arma::vec offsets = arma::zeros<arma::vec>(z.n_elem);
  Rcpp::List new_reg, new_keg_fit, new_theta_fit;
  double new_LL, curr_LL;
  arma::vec curr_para = para0;
  arma::vec new_para = arma::zeros<arma::vec>(para0.size());
  arma::vec curr_reg_par = arma::zeros<arma::vec>(pp);
  arma::vec new_reg_par = arma::zeros<arma::vec>(pp);
  arma::vec BETA = arma::zeros<arma::vec>(pp-1);
  arma::vec pis = arma::zeros<arma::vec>(z_AS.n_elem);
  double KAPPA, ETA, GAMMA, new_theta;
  double phi = 0.0;
  double curr_theta = 0.1;
  // initial regression fit
  
  while(iter < max_iter){
    
    if(H0 == 0){
      KAPPA = exp(new_para.at(0));
      ETA   = exp(new_para.at(1));
      GAMMA = exp(new_para.at(2));
    }else if(H0 == 1){ //para = log(c(KAPPA, GAMMA))
      KAPPA = exp(new_para.at(0));
      ETA   = 1.0;
      GAMMA = exp(new_para.at(1));
    }else{ //para = log(c(KAPPA, ETA))
      KAPPA = exp(new_para.at(0));
      ETA   = exp(new_para.at(1));
      GAMMA = 1.0;
    }
    //update THETA
    RcppT_compite_pi(z_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau, pis);
    
    new_theta_fit = RcppT_ase_theta_BFGS(curr_theta, ni, ni0, pis, lbc);
    new_theta     =  as<double>(new_theta_fit["PAR"]);
    
    //update BETA, phi
    RcppT_compute_offset(z, RHO, KAPPA, ETA, GAMMA, tau1, tau2, offsets);
    new_reg     = RcppT_reg_BFGS(y, X, offsets, curr_reg_par, lgy1,
                                 max_iter, eps, false);
    new_LL      = as<double>(new_reg["LL"]);
    new_reg_par = as<arma::vec>(new_reg["PAR"]);
    
    BETA = new_reg_par.subvec(0, X.n_cols-1);
    phi  = std::exp(new_reg_par.at(pp-1));
    
    if(show){
      Rprintf("TReCASE: BETA, PHI updated after %d iter \n",
              as<int>(new_reg["iter"]));
    }
    
    //update KEG
    new_keg_fit = RcppT_trecase_KEG_BFGS(curr_para, H0, y, z, z_AS, RHO, 
                                         RHO_AS, X, BETA, phi, tau1, 
                                         tau2, lgy1, ni0, ni, new_theta, 
                                         tauB, tau, lbc, max_iter, eps, false);
    new_para    = Rcpp::as<arma::vec>(new_keg_fit["PAR"]);
    new_LL      = as<double>(new_keg_fit["LL"]);
    
    if(show){
      Rprintf("TReC: keg updated after %d iter \n",
              as<int>(new_keg_fit["iter"]));
      // printR_obj(new_LL);
      // printR_obj(curr_LL);
    }
    
    if(iter > 0){
      if(abs(curr_LL - new_LL) < eps && 
         Rcpp_norm(curr_reg_par - new_reg_par) < eps &&
         Rcpp_norm(curr_para - new_para) < eps){
        // convergence criteria 
        if((curr_LL - new_LL > 0.0 && as<double>(new_keg_fit["norm_GRAD"]) > 0.01)){
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
    curr_para    = new_para;
    curr_LL      = new_LL;
    curr_theta   = new_theta;
    iter++;
    
  }
  if(show){
    Rprintf("TReCASE converges after %d iter \n", iter);
  }
  return Rcpp::List::create(
    Rcpp::Named("PAR", Rcpp::NumericVector(new_para.begin(), new_para.end())),
    Rcpp::Named("reg_par", Rcpp::NumericVector(new_reg_par.begin(), new_reg_par.end())),
    Rcpp::Named("log_theta", new_theta),
    Rcpp::Named("LL", new_LL),
    Rcpp::Named("converge", converge),
    Rcpp::Named("iter", iter)
  );
}

// [[Rcpp::export]]
Rcpp::List RcppT_trecase(const arma::vec& y, const arma::vec& z,
                         const arma::vec& z_AS, const arma::vec& RHO, 
                         const arma::vec& RHO_AS, const arma::mat& X, 
                         const arma::vec& tau1, const arma::vec& tau2,
                         const arma::vec& lgy1, const arma::vec& ni0, 
                         const arma::vec& ni, const arma::vec& tauB, 
                         const arma::vec& tau, const arma::vec& lbc, 
                         const arma::uword& max_iter = 4e3,
                         const double& eps = 1e-7, const bool& show = false){
  Rcpp::List sfit0, sfit1, sfit2, sfit3;
  double p_eta, p_gamma;
  arma::vec para0 = arma::zeros<arma::vec>(3);
  int converge = 0;
  
  sfit0 = RcppT_trecase_sfit(0, para0, y, z, z_AS, RHO, RHO_AS, X, tau1, tau2, 
                             lgy1, ni0, ni, tauB, tau, lbc, 
                             max_iter, eps, show); 
  sfit1 = RcppT_trecase_sfit(1, para0.subvec(0,1), y, z, z_AS, RHO, RHO_AS, X, tau1, 
                             tau2, lgy1, ni0, ni, tauB, tau, lbc, 
                             max_iter, eps, show); 
  sfit2 = RcppT_trecase_sfit(2, para0.subvec(0,1), y, z, z_AS, RHO, RHO_AS, X, tau1, 
                             tau2, lgy1, ni0, ni, tauB, tau, lbc, 
                             max_iter, eps, show); 
  
  para0.subvec(0,1) = as<arma::vec>(sfit1["PAR"]);
  para0.at(2) = as<arma::vec>(sfit1["PAR"]).at(2); 
  sfit3 = RcppT_trecase_sfit(0, para0, y, z, z_AS, RHO, RHO_AS, X, tau1, tau2, 
                             lgy1, ni0, ni, tauB, tau, lbc, 
                             max_iter, eps, show); 
  
  if(as<double>(sfit0["LL"]) < as<double>(sfit3["LL"]) && sfit3["converge"] ){
    sfit0 = sfit3;  
  }
  
  converge = sfit0["converge"] && sfit1["converge"] && sfit2["converge"];
  
  Rcpp::NumericVector PAR = as<Rcpp::NumericVector>(sfit0["PAR"]);
  Rcpp::NumericVector reg_par = as<Rcpp::NumericVector>(sfit0["reg_par"]);
  
  p_eta = R::pchisq(2*(as<double>(sfit0["LL"])-as<double>(sfit1["LL"])),
                    1, 0, 0);
  p_gamma = R::pchisq(2*(as<double>(sfit0["LL"])-as<double>(sfit2["LL"])),
                      1, 0, 0);
  
  return Rcpp::List::create(
    Rcpp::Named("p_eta", p_eta),
    Rcpp::Named("p_gamma", p_gamma),
    Rcpp::Named("PAR", Rcpp::NumericVector(PAR.begin(), PAR.end())),
    Rcpp::Named("reg_par", Rcpp::NumericVector(reg_par.begin(), reg_par.end())),
    Rcpp::Named("log_theta", sfit0["log_theta"]),
    Rcpp::Named("LL", sfit0["LL"]),
    Rcpp::Named("LL_eta", sfit1["LL"]),
    Rcpp::Named("LL_gamma", sfit2["LL"]),
    Rcpp::Named("converge", converge)
  );
  
}


// ----------------------------------------------------------------------
// TREC & ASE seperatly estimate 
// K, ETA, GAMMA, ETA_ase, GAMMA_ase (KEG_EaseGase / para5 )
// ----------------------------------------------------------------------


// [[Rcpp::export]]
double RcppT_TReC_ASE_LL_KEG(const arma::vec& KEG_EaseGase,
                             const arma::vec& y, const arma::vec& z,
                             const arma::vec& z_AS, const double& phi,
                             const arma::vec& RHO, const arma::vec& RHO_AS,
                             const arma::vec& tau1,const arma::vec& tau2,
                             const arma::vec& ni0, const arma::vec& ni,
                             const double& log_theta, const arma::vec& tauB,
                             const arma::vec& tau, const arma::vec& lgy1,
                             const arma::vec& expXbeta,const arma::vec& lbc,
                             arma::vec& offsets, 
                             arma::vec& pis, arma::vec& mu){
  double LL_NB, LL_BB;
  arma::vec paraBB = arma::zeros<arma::vec>(3);
  paraBB.at(0) = KEG_EaseGase.at(0);
  paraBB.subvec(1,2) = KEG_EaseGase.subvec(3,4);
  LL_NB = RcppT_loglikNB_KEG(KEG_EaseGase.subvec(0,2), 0, y, z, phi, RHO, tau1,
                             tau2, lgy1, expXbeta, offsets, mu);
  LL_BB = RcppT_loglikBB_KEG(paraBB, 0, z_AS, RHO_AS, ni, ni0, log_theta,
                             lbc, tauB, tau, pis);
  return LL_NB + LL_BB;
}



// [[Rcpp::export]]
arma::vec RcppT_TReC_ASE_grad_KEG(const arma::vec& KEG_EaseGase, 
                                  const arma::vec& y, const arma::vec& z,
                                  const arma::vec& z_AS, const double& phi,
                                  const arma::vec& RHO, const arma::vec& RHO_AS,
                                  const arma::vec& tau1,const arma::vec& tau2,
                                  const arma::vec& ni0, const arma::vec& ni,
                                  const double& log_theta, const arma::vec& tauB,
                                  const arma::vec& tau, const arma::vec& lgy1,
                                  const arma::vec& expXbeta,const arma::vec& lbc,
                                  const arma::vec& offsets, 
                                  const arma::vec& pis, const arma::vec& mu){
  arma::vec grad_NB = arma::zeros<arma::vec>(3);
  arma::vec grad_BB = arma::zeros<arma::vec>(3);
  arma::vec paraBB = arma::zeros<arma::vec>(3);
  arma::vec grad = arma::zeros<arma::vec>(5);
  paraBB.at(0) = KEG_EaseGase.at(0);
  paraBB.subvec(1,2) = KEG_EaseGase.subvec(3,4);
  grad_NB = RcppT_grad_NB(KEG_EaseGase.subvec(0,2), 0, y, z, RHO, phi, 
                          tau1, tau2, expXbeta, offsets, mu);
  
  grad_BB = RcppT_grad_BB_KEG(paraBB, 0, z_AS, RHO_AS, ni, ni0, log_theta, 
                              lbc, tauB, tau, pis);
  
  grad.at(0) = grad_NB.at(0) + grad_BB.at(0);
  grad.subvec(1,2) = grad_NB.subvec(1,2);
  grad.subvec(3,4) = grad_BB.subvec(1,2);
  
  return grad;
  
}


// [[Rcpp::export]]
Rcpp::List RcppT_trec_ase_KEG_BFGS(const arma::vec& para0, 
                                   const arma::vec& y, const arma::vec& z, 
                                   const arma::vec& z_AS,
                                   const arma::vec& RHO, const arma::vec& RHO_AS,
                                   const arma::mat& X, const arma::vec& BETA, 
                                   const double& phi,
                                   const arma::vec& tau1, const arma::vec& tau2,
                                   const arma::vec& lgy1,  
                                   const arma::vec& ni0, const arma::vec& ni,
                                   const double& log_theta, const arma::vec& tauB,
                                   const arma::vec& tau, const arma::vec& lbc,
                                   const arma::uword& max_iter = 4e3,
                                   const double& eps = 1e-7,const bool& show = true){
  
  arma::uword num_params = para0.size();
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
  arma::vec offsets = arma::zeros<arma::vec>(y.n_elem);
  arma::vec expXbeta = arma::zeros<arma::vec>(y.n_elem);
  arma::vec pis = arma::zeros<arma::vec>(z_AS.n_elem);
  
  double old_LL,new_LL,inv_norm_p_k,tmp_alpha,ys;
  double fnscale = -1.0; // For maximization
  double curr_LL = 0.0;
  xk = para0;
  
  RcppT_compute_expXbeta(X, BETA, expXbeta); 
  
  while(iter < max_iter){
    //calculate direction p_k
    uu = 0;
    old_LL = fnscale * RcppT_TReC_ASE_LL_KEG(xk,y, z, z_AS, phi, RHO, RHO_AS, 
                                             tau1, tau2, ni0, ni, log_theta, tauB, 
                                             tau, lgy1, expXbeta, lbc, offsets, 
                                             pis, mu);
    gr_k   = fnscale * RcppT_TReC_ASE_grad_KEG(xk, y, z, z_AS, phi, RHO, 
                                               RHO_AS, tau1, tau2, ni0, ni, 
                                               log_theta, tauB, tau, lgy1, expXbeta, 
                                               lbc, offsets, pis, mu);
    p_k    = -1.0 * inv_Bk * gr_k;
    inv_norm_p_k =  1.0 / std::max(1.0, Rcpp_norm(p_k));
    
    //line search for new xk
    for(jj=0; jj<15; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4, jj);
      new_xk    = xk + tmp_alpha * p_k;
      new_LL    = fnscale * RcppT_TReC_ASE_LL_KEG(new_xk, y, z, z_AS, phi,
                                                  RHO, RHO_AS, tau1, tau2, 
                                                  ni0, ni, log_theta, tauB, 
                                                  tau, lgy1, expXbeta, lbc,
                                                  offsets, pis, mu);
      if(new_LL < old_LL){ //minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * RcppT_TReC_ASE_grad_KEG(new_xk, y, z, z_AS, phi, RHO, 
                                                RHO_AS, tau1, tau2, ni0, ni, 
                                                log_theta, tauB, tau, lgy1,
                                                expXbeta, lbc, offsets, pis,
                                                mu) - gr_k;
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
        gr_k = RcppT_TReC_ASE_grad_KEG(xk, y, z, z_AS, phi, RHO, RHO_AS, 
                                       tau1, tau2, ni0, ni, log_theta, tauB, tau, 
                                       lgy1, expXbeta, lbc, offsets, pis, mu);
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
  
  old_LL = RcppT_TReCASE_LL(y, phi, lgy1, mu, ni0, ni,log_theta, pis, lbc);
  return Rcpp::List::create(
    Rcpp::Named("converge", converge),
    Rcpp::Named("LL", old_LL),
    Rcpp::Named("iter", iter),
    Rcpp::Named("norm_GRAD", Rcpp_norm(gr_k)),
    Rcpp::Named("PAR", Rcpp::NumericVector(xk.begin(), xk.end()))
  );
}


// [[Rcpp::export]]
Rcpp::List RcppT_trec_ase(const arma::vec& para0,
                          const arma::vec& y, 
                          const arma::vec& z, const arma::vec& z_AS,
                          const arma::vec& RHO, const arma::vec& RHO_AS,
                          const arma::mat& X, 
                          const arma::vec& tau1, const arma::vec& tau2,
                          const arma::vec& lgy1,  
                          const arma::vec& ni0, const arma::vec& ni,
                          const arma::vec& tauB,const arma::vec& tau, 
                          const arma::vec& lbc,const arma::uword& max_iter = 4e3,
                          const double& eps = 1e-7, const bool& show = false){
  arma::uword iter = 0;
  arma::uword converge = 0;
  arma::uword pp = X.n_cols + 1;
  arma::vec offsets = arma::zeros<arma::vec>(z.n_elem);
  Rcpp::List new_reg, new_keg_fit, new_theta_fit;
  double new_LL, curr_LL;
  arma::vec curr_para = para0;
  arma::vec new_para = arma::zeros<arma::vec>(para0.size());
  arma::vec curr_reg_par = arma::zeros<arma::vec>(pp);
  arma::vec new_reg_par = arma::zeros<arma::vec>(pp);
  arma::vec BETA = arma::zeros<arma::vec>(pp-1);
  arma::vec pis = arma::zeros<arma::vec>(z_AS.n_elem);
  double KAPPA, ETA, GAMMA, ETA_ase, GAMMA_ase, new_theta;
  double phi = 0.0;
  double curr_theta = 0.1;
  // initial regression fit
  
  while(iter < max_iter){
    
    
    KAPPA = exp(new_para.at(0));
    ETA   = exp(new_para.at(1));
    GAMMA = exp(new_para.at(2));
    ETA_ase   = exp(new_para.at(3));
    GAMMA_ase = exp(new_para.at(4));
    
    //update THETA
    RcppT_compite_pi(z_AS, RHO_AS, KAPPA, ETA_ase, GAMMA_ase, tauB, tau, pis);
    
    new_theta_fit = RcppT_ase_theta_BFGS(curr_theta, ni, ni0, pis, lbc);
    new_theta     =  as<double>(new_theta_fit["PAR"]);
    
    //update BETA, phi
    RcppT_compute_offset(z, RHO, KAPPA, ETA, GAMMA, tau1, tau2, offsets);
    new_reg     = RcppT_reg_BFGS(y, X, offsets, curr_reg_par, lgy1,
                                 max_iter, eps, false);
    new_LL      = as<double>(new_reg["LL"]);
    new_reg_par = as<arma::vec>(new_reg["PAR"]);
    
    BETA = new_reg_par.subvec(0, X.n_cols-1);
    phi  = std::exp(new_reg_par.at(pp-1));
    
    if(show){
      Rprintf("TReC_ASE: BETA, PHI updated after %d iter \n",
              as<int>(new_reg["iter"]));
    }
    
    //update KEG
    new_keg_fit = RcppT_trec_ase_KEG_BFGS(curr_para, y, z, z_AS, RHO, 
                                          RHO_AS, X, BETA, phi, tau1, 
                                          tau2, lgy1, ni0, ni, new_theta, 
                                          tauB, tau, lbc, max_iter, eps, false);
    new_para    = Rcpp::as<arma::vec>(new_keg_fit["PAR"]);
    new_LL      = as<double>(new_keg_fit["LL"]);
    
    if(show){
      Rprintf("TReC_ASE: keg updated after %d iter \n",
              as<int>(new_keg_fit["iter"]));
      // printR_obj(new_LL);
      // printR_obj(curr_LL);
    }
    
    if(iter > 0){
      if(abs(curr_LL - new_LL) < eps && 
         Rcpp_norm(curr_reg_par - new_reg_par) < eps &&
         Rcpp_norm(curr_para - new_para) < eps){
        // convergence criteria 
        if((curr_LL - new_LL > 0.0 && as<double>(new_keg_fit["norm_GRAD"]) > 0.01)){
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
    curr_para    = new_para;
    curr_LL      = new_LL;
    curr_theta   = new_theta;
    iter++;
    
  }
  if(show){
    Rprintf("TReC_ASE: converges after %d iter \n", iter);
  }
  return Rcpp::List::create(
    Rcpp::Named("PAR", Rcpp::NumericVector(new_para.begin(), new_para.end())),
    Rcpp::Named("reg_par", Rcpp::NumericVector(new_reg_par.begin(), new_reg_par.end())),
    Rcpp::Named("log_theta", new_theta),
    Rcpp::Named("LL", new_LL),
    Rcpp::Named("converge", converge),
    Rcpp::Named("iter", iter)
  );
}


// ----------------------------------------------------------------------
// Cis-Trans Score test
// ----------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List RcppT_CisTrans_ScoreObs(const arma::vec& para, const arma::vec& y,
                                   const arma::vec& z, const arma::vec& z_AS,
                                   const arma::vec& RHO, const arma::vec& RHO_AS,
                                   const arma::mat& X, const arma::vec& BETA,
                                   const double& phi,
                                   const arma::vec& tau1, const arma::vec& tau2,
                                   const arma::vec& lgy1,
                                   const arma::vec& ni0, const arma::vec& ni,
                                   const double& log_theta, const arma::vec& tauB,
                                   const arma::vec& tau, const arma::vec& lbc){
  
  double KAPPA, ETA, GAMMA;
  arma::uword ii;
  arma::uword n= RHO.n_elem;
  arma::uword n_AS = RHO_AS.n_elem;
  arma::uword beta_npara = BETA.n_elem;
  
  double vphi = 1/phi;
  double divphi = R::digamma(vphi);
  double trvphi = R::trigamma(vphi);
  //TReC
  arma::vec offsets  = arma::zeros<arma::vec>(n);
  arma::vec expXbeta = arma::zeros<arma::vec>(n);
  arma::vec mu       = arma::zeros<arma::vec>(n);
  arma::mat dmu_keg = arma::zeros<arma::mat>(1,3);
  //information matrix
  arma::mat Iee = arma::zeros<arma::mat>(3,3);  //KEG
  arma::mat Ibb = arma::zeros<arma::mat>(beta_npara, beta_npara);  //beta
  arma::mat Ibe = arma::zeros<arma::mat>(beta_npara, 3);  //beta*keg
  arma::mat Iep = arma::zeros<arma::mat>(3, 1);  //beta*phi
  arma::mat Ibp = arma::zeros<arma::mat>(beta_npara, 1);  //beta*keg
  double    Ipp = 0.0; // phi
  
  double dltrec_dmu, d2ltrec_dmu, phi_mu_1, phi_y_1,
  dmu_dkappa, dmu_deta, dmu_dgamma, d2mu_dkappa_dgamma;
  
  //ASE
  double vtheta = 1/exp(log_theta);
  double divtheta = R::digamma(vtheta);
  double trvtheta = R::trigamma(vtheta);
  double Itt = 0.0;
  arma::vec pis  = arma::zeros<arma::vec>(n_AS);
  arma::vec score_alpha  = arma::zeros<arma::vec>(2);
  arma::mat Iaa = arma::zeros<arma::mat>(2, 2);  //alpha
  arma::mat Iae = arma::zeros<arma::mat>(2, 3);  //alpha*KEG
  arma::mat Iat = arma::zeros<arma::mat>(2, 1);  //theta*KEG
  arma::mat Iet = arma::zeros<arma::mat>(3, 1);  //theta*KEG
  
  double dlASE_dpi, d2lASE_dpi, tmp1, tmp2, tmp3, tmp4, tmp_kappa, 
  dpi_dkappa, dpi_deta, dpi_dgamma,
  d2pi_dkappa, d2pi_deta, d2pi_dgamma,
  d2pi_dkappa_eta, d2pi_dkappa_gamma, d2pi_deta_gamma; 
  
  double W, dW_dTHETA, dlASE_dTHETA, d2lASE_dTHETA_dpi, pval;
  KAPPA = exp(para.at(0));
  ETA   = exp(para.at(1));
  GAMMA = exp(para.at(2));
  
  //final information matrix
  arma::mat M1    = arma::zeros<arma::mat>(beta_npara+5, beta_npara+5); 
  arma::mat M2    = arma::zeros<arma::mat>(beta_npara+5, 2); 
  arma::mat OIMat = arma::zeros<arma::mat>(beta_npara+7, beta_npara+7); 
  arma::mat Score = arma::zeros<arma::mat>(1, 1);
  
  /*
   * Information matrix (TReC)
   */
  
  RcppT_compute_offset(z, RHO, KAPPA, ETA, GAMMA, tau1, tau2, offsets);
  RcppT_compute_expXbeta(X, BETA, expXbeta);
  mu = exp(offsets) % expXbeta;
  
  for(ii =0; ii<n; ii++){
    phi_mu_1    = 1.0 + phi*mu.at(ii);
    phi_y_1     = 1.0 + phi*y.at(ii);
    dltrec_dmu  = y.at(ii)/mu.at(ii) - phi_y_1/phi_mu_1;
    d2ltrec_dmu = -y.at(ii)/pow(mu.at(ii), 2.0) + phi*phi_y_1/pow(phi_mu_1,2);
    
    if(z.at(ii) == 0){
      dmu_dkappa = (tau1.at(ii) + tau2.at(ii))*expXbeta.at(ii)*RHO.at(ii);
      dmu_deta   = 0;
      dmu_dgamma = 0;
      
      d2mu_dkappa_dgamma = 0;
      
    }else if(z.at(ii) == 1){
      
      dmu_dkappa = (tau1.at(ii) + tau2.at(ii)*GAMMA)*expXbeta.at(ii)*RHO.at(ii);
      dmu_deta   = (1-RHO.at(ii))*expXbeta.at(ii);
      dmu_dgamma = tau2.at(ii)*RHO.at(ii)*KAPPA*expXbeta.at(ii);
      
      d2mu_dkappa_dgamma = tau2.at(ii)*RHO.at(ii)*expXbeta.at(ii);
      
    }else if(z.at(ii) == 2){
      
      dmu_dkappa = (tau1.at(ii)*GAMMA + tau2.at(ii))*expXbeta.at(ii)*RHO.at(ii);
      dmu_deta   = (1-RHO.at(ii))*expXbeta.at(ii);
      dmu_dgamma = tau1.at(ii)*RHO.at(ii)*KAPPA*expXbeta.at(ii);
      
      d2mu_dkappa_dgamma = tau1.at(ii)*RHO.at(ii)*expXbeta.at(ii);
      
    }else{
      dmu_dkappa = (tau1.at(ii)+tau2.at(ii))*GAMMA*expXbeta.at(ii)*RHO.at(ii);
      dmu_deta   = 2*(1-RHO.at(ii))*expXbeta.at(ii);
      dmu_dgamma = (tau1.at(ii)+tau2.at(ii))*KAPPA*expXbeta.at(ii)*RHO.at(ii);
      
      d2mu_dkappa_dgamma = (tau1.at(ii)+tau2.at(ii))*RHO.at(ii)*expXbeta.at(ii);
      
    }
    dmu_keg.at(0,0) = dmu_dkappa;
    dmu_keg.at(0,1) = dmu_deta;
    dmu_keg.at(0,2) = dmu_dgamma;
    
    Iee.at(0,0) += d2ltrec_dmu*pow(dmu_dkappa, 2.0);  // kappa^2
    Iee.at(1,1) += d2ltrec_dmu*pow(dmu_deta, 2.0);    // eta^2
    Iee.at(2,2) += d2ltrec_dmu*pow(dmu_dgamma, 2.0);  // gamma^2
    
    Iee.at(0,1) += d2ltrec_dmu*dmu_dkappa*dmu_deta; // kappa*eta
    Iee.at(0,2) += d2ltrec_dmu*dmu_dkappa*dmu_dgamma + 
      dltrec_dmu*d2mu_dkappa_dgamma;                //kappa*gamma
    Iee.at(1,2) += d2ltrec_dmu*dmu_deta*dmu_dgamma; //eta*gamma
    
    Ipp += 2.0*pow(vphi, 3.0)*(R::digamma(y.at(ii)+vphi)-divphi-log(phi_mu_1)) +
      pow(vphi, 4.0)*(R::trigamma(y.at(ii)+vphi)-trvphi) +
      pow(vphi, 2.0)*(2.0*mu.at(ii)/phi_mu_1 - y.at(ii)) +
      (vphi+y.at(ii))*pow(mu.at(ii),2.0)/pow(phi_mu_1,2.0);
    
    Ibb -= (mu.at(ii)/phi_mu_1 + 
      phi*mu.at(ii)*(y.at(ii)-mu.at(ii))/pow(phi_mu_1,2.0))*
      X.row(ii).t()*X.row(ii);
    
    Iep -= (y.at(ii)-mu.at(ii))/pow(phi_mu_1,2.0)*dmu_keg.t();
    
    Ibe -= (1/phi_mu_1+(y.at(ii)-mu.at(ii))*phi/pow(phi_mu_1,2.0))*
      X.row(ii).t()*dmu_keg;
    
    Ibp -= ((y.at(ii)-mu.at(ii))*mu.at(ii)/pow(phi_mu_1,2.0))*X.row(ii).t();
  }
  
  // printR_obj(Iee);
  // printR_obj(Ipp);
  // printR_obj(Ibb);
  // printR_obj(Ibp);
  // printR_obj(Ibe);
  // printR_obj(Iep);
  
  /*
   * Information matrix (ASE)
   */
  
  RcppT_compite_pi(z_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau, pis);
  
  for(ii =0; ii<n_AS; ii++){
    double aa = pis.at(ii) * vtheta; 
    double bb = vtheta - aa;
    double diaa = R::digamma(aa);
    double dibb =  R::digamma(bb);
    double diaa_ni0 = R::digamma(aa + ni0.at(ii));
    double dibb_ni1 = R::digamma(bb + ni.at(ii) - ni0.at(ii));
    double triaa = R::trigamma(aa);
    double tribb =  R::trigamma(bb);
    double triaa_ni0 = R::trigamma(aa + ni0.at(ii));
    double tribb_ni1 = R::trigamma(bb + ni.at(ii) - ni0.at(ii));
    
    dlASE_dpi  = vtheta*(diaa_ni0 - dibb_ni1 - diaa + dibb);
    d2lASE_dpi = pow(vtheta, 2.0)*(triaa_ni0 + tribb_ni1 - triaa -
      tribb);
    
    if(z_AS.at(ii) == 0 | z_AS.at(ii) == 3){
      tmp1 = RHO_AS.at(ii)*tau.at(ii)*KAPPA + 2*(1 - RHO_AS.at(ii));
      tmp2 = RHO_AS.at(ii)*tauB.at(ii)/tmp1;
      tmp3 = RHO_AS.at(ii)*tau.at(ii)*(RHO_AS.at(ii)*tauB.at(ii)*KAPPA + 
        (1 - RHO_AS.at(ii)))/pow(tmp1, 2.0);
      dpi_dkappa = tmp2 - tmp3;
      dpi_deta   = 0;
      dpi_dgamma = 0;
      
      d2pi_dkappa = -(2.0*pow(RHO_AS.at(ii),2.0)*(1.0-RHO_AS.at(ii))*
        (2.0*tauB.at(ii) -tau.at(ii))*tau.at(ii))/pow(tmp1, 3.0);
      d2pi_deta = d2pi_dgamma  = d2pi_dkappa_eta = 
        d2pi_dkappa_gamma = d2pi_deta_gamma = 0;
      
    }else{
      tmp1 = RHO_AS.at(ii)*tauB.at(ii)*GAMMA*KAPPA + (1 - RHO_AS.at(ii))*ETA;
      tmp2 = RHO_AS.at(ii)*(tau.at(ii) - tauB.at(ii))*KAPPA +
        (1 - RHO_AS.at(ii)) + tmp1;
      tmp3 = tmp1/pow(tmp2, 2.0);
      tmp4 = 1/tmp2 -tmp3; 
      dpi_dkappa = RHO_AS.at(ii)*tauB.at(ii)*GAMMA*tmp4 -
        RHO_AS.at(ii)*(tau.at(ii) - tauB.at(ii))*tmp3;
      dpi_deta   = (1 - RHO_AS.at(ii))*tmp4;
      dpi_dgamma = RHO_AS.at(ii)*tauB.at(ii)*KAPPA*tmp4;
      
      tmp3 = tmp1/pow(tmp2, 3.0); // changed tmp3
      tmp_kappa = RHO_AS.at(ii)*(tau.at(ii) - tauB.at(ii)) + 
        RHO_AS.at(ii)* tauB.at(ii)*GAMMA;
      d2pi_dkappa = -2.0*GAMMA*RHO_AS.at(ii)*tauB.at(ii)*tmp_kappa/pow(tmp2,2.0) +
        2.0*pow(tmp_kappa,2.0)*tmp3;
      d2pi_deta   = -2.0*pow(1-RHO_AS.at(ii), 2.0)*(1/pow(tmp2, 2.0) - tmp3);
      d2pi_dgamma = -2.0*pow(RHO_AS.at(ii)*tauB.at(ii)*KAPPA,2.0)*
        (1/pow(tmp2, 2.0) - tmp3);
      d2pi_dkappa_eta  = (-RHO_AS.at(ii)*tauB.at(ii)*GAMMA*(1-RHO_AS.at(ii)) -
        tmp_kappa*(1-RHO_AS.at(ii)))/pow(tmp2, 2.0) +
        2.0*(1-RHO_AS.at(ii))*tmp_kappa*tmp3;
      d2pi_dkappa_gamma = RHO_AS.at(ii)*tauB.at(ii)/tmp2 -
        (pow(RHO_AS.at(ii)*tauB.at(ii), 2.0)*KAPPA*GAMMA + 
        tmp_kappa*RHO_AS.at(ii)*tauB.at(ii)*KAPPA)/pow(tmp2, 2.0) -
        RHO_AS.at(ii)*tauB.at(ii)*tmp3*tmp2 +
        2.0*tmp_kappa*RHO_AS.at(ii)*tauB.at(ii)*KAPPA*tmp3;
      d2pi_deta_gamma   = -2.0*KAPPA*(1-RHO_AS.at(ii))*RHO_AS.at(ii)*tauB.at(ii)/pow(tmp2,2.0) +
        2.0*KAPPA*(1.0-RHO_AS.at(ii))*RHO_AS.at(ii)*tauB.at(ii)*tmp3;
    }
    
    Iee.at(0,0) += d2lASE_dpi*pow(dpi_dkappa, 2.0)  + dlASE_dpi*d2pi_dkappa;
    Iae.at(0,1) += d2lASE_dpi*pow(dpi_deta, 2.0)    + dlASE_dpi*d2pi_deta;
    Iae.at(1,2) += d2lASE_dpi*pow(dpi_dgamma, 2.0)  + dlASE_dpi*d2pi_dgamma; 
    Iae.at(0,0) += d2lASE_dpi*dpi_dkappa*dpi_deta   + dlASE_dpi*d2pi_dkappa_eta; 
    Iae.at(1,0) += d2lASE_dpi*dpi_dkappa*dpi_dgamma + dlASE_dpi*d2pi_dkappa_gamma;
    Iae.at(0,2) += d2lASE_dpi*dpi_deta*dpi_dgamma   + dlASE_dpi*d2pi_deta_gamma; 
    
    score_alpha.at(0) += dlASE_dpi*dpi_deta;
    score_alpha.at(1) += dlASE_dpi*dpi_dgamma;
    
    W = -(R::digamma(ni.at(ii)+vtheta)+diaa*pis.at(ii) +dibb*(1.0-pis.at(ii)) -
      diaa_ni0*pis.at(ii)-dibb_ni1*(1.0-pis.at(ii)) - divtheta);
    dW_dTHETA = -pow(vtheta, 2.0)*(pow(pis.at(ii), 2.0)*(triaa_ni0 - triaa) +
      pow(1-pis.at(ii), 2.0)*(tribb_ni1 - tribb) +
      trvtheta - R::trigamma(vtheta+ni.at(ii)));
    Itt += 2.0*pow(vtheta,3.0)*W - pow(vtheta, 2.0)*dW_dTHETA;
    
    
    d2lASE_dTHETA_dpi = - pow(vtheta, 2.0)*(diaa_ni0 - diaa -
      dibb_ni1 + dibb) -
      pow(vtheta,3.0)*pis.at(ii)*(triaa_ni0 - triaa) +
      pow(vtheta,3.0)*(1.0-pis.at(ii))*(tribb_ni1 - tribb);
    // printR_obj(d2lASE_dTHETA_dpi);
    Iet.at(0,0) += d2lASE_dTHETA_dpi*dpi_dkappa;
    Iet.at(1,0) += d2lASE_dTHETA_dpi*dpi_deta;
    Iet.at(2,0) += d2lASE_dTHETA_dpi*dpi_dgamma;
  }
  // printR_obj(Iee);
  
  Iae.at(1,1) = Iae.at(0,2);
  Iaa = Iae.submat(0,1,1,2);
  
  Iee.at(0,1) += Iae.at(0,0);
  Iee.at(1,1) += Iae.at(0,1);
  Iee.at(1,2) += Iae.at(0,2);
  Iee.at(0,2) += Iae.at(1,0);
  Iee.at(2,2) += Iae.at(1,2);
  
  Iee.at(1,0) = Iee.at(0,1);
  Iee.at(2,0) = Iee.at(0,2);
  Iee.at(2,1) = Iee.at(1,2);
  
  Iat = Iet.submat(1,0,2,0);
  // printR_obj(score_alpha);
  // printR_obj(Iee);
  // printR_obj(Iae);
  // printR_obj(Iaa);
  // printR_obj(Itt);
  // printR_obj(Iet);
  // printR_obj(Iat);
  
  /*
   * Final information matrix
   */
  
  arma::uword np = beta_npara -1;
  M1.submat(0,0,np,np)             = Ibb;
  M1.submat(0,np+1,np,np+3)        = Ibe;
  M1.submat(np+1,0,np+3,np)        = Ibe.t();
  M1.submat(0,np+4,np,np+4)        = Ibp;
  M1.submat(np+4,0,np+4,np)        = Ibp.t();
  M1.submat(np+1,np+1,np+3,np+3)   = Iee;
  M1.submat(np+1,np+5,np+3,np+5)   = Iet;
  M1.submat(np+5,np+1,np+5,np+3)   = Iet.t();
  M1.submat(np+1,np+4,np+3,np+4)   = Iep;
  M1.submat(np+4,np+1,np+4,np+3)   = Iep.t();
  M1.at((np+5),(np+5))             = Itt;
  M1.at((np+4),(np+4))             = Ipp;
  M2.rows(np+1,np+3) = Iae.t();
  M2.row(5+np)    = Iat.t();
  // printR_obj(M1);
  // printR_obj(M2);
  
  M1 = -M1;
  M2 = -M2;
  Iaa = -Iaa;
  
  
  OIMat.submat(0,0,(np+5),(np+5))           = M1;
  OIMat.submat(0,(np+6),(np+5),(np+7))      = M2;
  OIMat.submat((np+6),0,(np+7),(np+5))      = M2.t();
  OIMat.submat((np+6),(np+6),(np+7),(np+7)) = Iaa;
  
  
  // printR_obj(OIMat);
  if(!OIMat.is_sympd()){
    Score.at(0,0) = -6.0;
    pval  = 1.0;
  }else{
    Score = score_alpha.t()*((Iaa - M2.t()*M1.i()*M2).i())*score_alpha;
    pval  = R::pchisq(Score.at(0,0),2,0,0);
  }
  
  if(std::isnan(Score.at(0,0))){
    Score.at(0,0) = -6.0;
    pval = 1.0;
  }
  return Rcpp::List::create(
    Rcpp::Named("Score", Score.at(0,0)),
    Rcpp::Named("pval", pval)
  );
}


// [[Rcpp::export]]
void RcppT_ASE_ExpFunc(const double& ni, const double& pis_i, 
                       const double& vtheta, arma::vec& Expvec){
  //   Expvec[0] = Expected value of digamma(vTHETA*pi+A)
  //   Expvec[1] = Expected value of digamma(vTHETA*(1.0-pi)+D-A)
  //   Expvec[2] = Expected value of trigamma(vTHETA*pi+A)
  //   Expvec[3] = Expected value of trigamma(vTHETA*(1.0-pi)+D-A)
  int ni00;
  double lpmf, pmf;
  for(ni00 = 0; ni00<(ni+1); ni00++){
    lpmf = R::lchoose(ni, ni00) + lgamma(ni00+pis_i*vtheta) + 
      lgamma(ni-ni00+(1.0-pis_i)*vtheta) + lgamma(vtheta) -
      lgamma(ni+vtheta)-lgamma(pis_i*vtheta)-
      lgamma((1.0-pis_i)*vtheta);
    pmf = exp(lpmf);
    Expvec.at(0) += R::digamma(vtheta*pis_i + ni00)*pmf;
    Expvec.at(1) += R::digamma(vtheta*(1.0-pis_i) + ni-ni00)*pmf;
    Expvec.at(2) += R::trigamma(vtheta*pis_i + ni00)*pmf;
    Expvec.at(3) += R::trigamma(vtheta*(1.0-pis_i) + ni-ni00)*pmf;
  }
  
}


// [[Rcpp::export]]
Rcpp::List RcppT_CisTrans_Score(const arma::vec& para, const arma::vec& y,
                                const arma::vec& z, const arma::vec& z_AS,
                                const arma::vec& RHO, const arma::vec& RHO_AS,
                                const arma::mat& X, const arma::vec& BETA,
                                const double& phi,
                                const arma::vec& tau1, const arma::vec& tau2,
                                const arma::vec& lgy1,
                                const arma::vec& ni0, const arma::vec& ni,
                                const double& log_theta, const arma::vec& tauB,
                                const arma::vec& tau, const arma::vec& lbc){
  
  double KAPPA, ETA, GAMMA;
  arma::uword ii;
  arma::uword n= RHO.n_elem;
  arma::uword n_AS = RHO_AS.n_elem;
  arma::uword beta_npara = BETA.n_elem;
  
  double vphi = 1/phi;
  double divphi = R::digamma(vphi);
  double trvphi = R::trigamma(vphi);
  //TReC
  arma::vec offsets  = arma::zeros<arma::vec>(n);
  arma::vec expXbeta = arma::zeros<arma::vec>(n);
  arma::vec mu       = arma::zeros<arma::vec>(n);
  arma::mat dmu_keg = arma::zeros<arma::mat>(1,3);
  //information matrix
  arma::mat Iee = arma::zeros<arma::mat>(3,3);  //KEG
  arma::mat Ibb = arma::zeros<arma::mat>(beta_npara, beta_npara);  //beta
  arma::mat Ibe = arma::zeros<arma::mat>(beta_npara, 3);  //beta*keg
  arma::mat Iep = arma::zeros<arma::mat>(3, 1);  //beta*phi
  arma::mat Ibp = arma::zeros<arma::mat>(beta_npara, 1);  //beta*keg
  double    Ipp = 0.0; // phi
  
  double dltrec_dmu, d2ltrec_dmu, phi_mu_1, phi_y_1,
  dmu_dkappa, dmu_deta, dmu_dgamma, d2mu_dkappa_dgamma;
  
  //ASE
  double vtheta = 1/exp(log_theta);
  double divtheta = R::digamma(vtheta);
  double trvtheta = R::trigamma(vtheta);
  double Itt = 0.0;
  arma::vec pis  = arma::zeros<arma::vec>(n_AS);
  arma::vec score_alpha  = arma::zeros<arma::vec>(2);
  arma::mat Iaa = arma::zeros<arma::mat>(2, 2);  //alpha
  arma::mat Iae = arma::zeros<arma::mat>(2, 3);  //alpha*KEG
  arma::mat Iat = arma::zeros<arma::mat>(2, 1);  //theta*KEG
  arma::mat Iet = arma::zeros<arma::mat>(3, 1);  //theta*KEG
  
  //   Expvec[1] = Expected value of digamma(vTHETA*pi+A)
  //   Expvec[2] = Expected value of digamma(vTHETA*(1.0-pi)+D-A)
  //   Expvec[3] = Expected value of trigamma(vTHETA*pi+A)
  //   Expvec[4] = Expected value of trigamma(vTHETA*(1.0-pi)+D-A)
  arma::vec Expvec = arma::zeros<arma::vec>(4);  
  
  double dlASE_dpi, d2lASE_dpi, tmp1, tmp2, tmp3, tmp4, tmp_kappa, 
  dpi_dkappa, dpi_deta, dpi_dgamma,
  d2pi_dkappa, d2pi_deta, d2pi_dgamma,
  d2pi_dkappa_eta, d2pi_dkappa_gamma, d2pi_deta_gamma; 
  
  double W, dW_dTHETA, dlASE_dTHETA, d2lASE_dTHETA_dpi, pval;
  KAPPA = exp(para.at(0));
  ETA   = exp(para.at(1));
  GAMMA = exp(para.at(2));
  
  //final information matrix
  arma::mat M1    = arma::zeros<arma::mat>(beta_npara+4, beta_npara+4); 
  arma::mat M2    = arma::zeros<arma::mat>(beta_npara+4, 2); 
  arma::mat OIMat = arma::zeros<arma::mat>(beta_npara+6, beta_npara+6); 
  arma::mat Score = arma::zeros<arma::mat>(1, 1);
  
  /*
   * Information matrix (TReC)
   */
  
  RcppT_compute_offset(z, RHO, KAPPA, ETA, GAMMA, tau1, tau2, offsets);
  RcppT_compute_expXbeta(X, BETA, expXbeta);
  mu = exp(offsets) % expXbeta;
  
  for(ii =0; ii<n; ii++){
    phi_mu_1    = 1.0 + phi*mu.at(ii);
    phi_y_1     = 1.0 + phi*y.at(ii);
    // dltrec_dmu  = y.at(ii)/mu.at(ii) - (1.0 + phi*y.at(ii))/phi_mu_1;
    d2ltrec_dmu = -1.0/(mu.at(ii)*phi_mu_1);
    
    if(z.at(ii) == 0){
      dmu_dkappa = (tau1.at(ii) + tau2.at(ii))*expXbeta.at(ii)*RHO.at(ii);
      dmu_deta   = 0;
      dmu_dgamma = 0;
      
      d2mu_dkappa_dgamma = 0;
      
    }else if(z.at(ii) == 1){
      
      dmu_dkappa = (tau1.at(ii) + tau2.at(ii)*GAMMA)*expXbeta.at(ii)*RHO.at(ii);
      dmu_deta   = (1-RHO.at(ii))*expXbeta.at(ii);
      dmu_dgamma = tau2.at(ii)*RHO.at(ii)*KAPPA*expXbeta.at(ii);
      
      d2mu_dkappa_dgamma = tau2.at(ii)*RHO.at(ii)*expXbeta.at(ii);
      
    }else if(z.at(ii) == 2){
      
      dmu_dkappa = (tau1.at(ii)*GAMMA + tau2.at(ii))*expXbeta.at(ii)*RHO.at(ii);
      dmu_deta   = (1-RHO.at(ii))*expXbeta.at(ii);
      dmu_dgamma = tau1.at(ii)*RHO.at(ii)*KAPPA*expXbeta.at(ii);
      
      d2mu_dkappa_dgamma = tau1.at(ii)*RHO.at(ii)*expXbeta.at(ii);
      
    }else{
      dmu_dkappa = (tau1.at(ii)+tau2.at(ii))*GAMMA*expXbeta.at(ii)*RHO.at(ii);
      dmu_deta   = 2*(1-RHO.at(ii))*expXbeta.at(ii);
      dmu_dgamma = (tau1.at(ii)+tau2.at(ii))*KAPPA*expXbeta.at(ii)*RHO.at(ii);
      
      d2mu_dkappa_dgamma = (tau1.at(ii)+tau2.at(ii))*RHO.at(ii)*expXbeta.at(ii);
      
    }
    dmu_keg.at(0,0) = dmu_dkappa;
    dmu_keg.at(0,1) = dmu_deta;
    dmu_keg.at(0,2) = dmu_dgamma;
    
    Iee.at(0,0) += d2ltrec_dmu*pow(dmu_dkappa, 2.0);  // kappa^2
    Iee.at(1,1) += d2ltrec_dmu*pow(dmu_deta, 2.0);    // eta^2
    Iee.at(2,2) += d2ltrec_dmu*pow(dmu_dgamma, 2.0);  // gamma^2
    
    Iee.at(0,1) += d2ltrec_dmu*dmu_dkappa*dmu_deta; // kappa*eta
    Iee.at(0,2) += d2ltrec_dmu*dmu_dkappa*dmu_dgamma;              //kappa*gamma
    Iee.at(1,2) += d2ltrec_dmu*dmu_deta*dmu_dgamma; //eta*gamma
    
    // Ipp += 2.0*pow(vphi, 3.0)*(R::digamma(y.at(ii)+vphi)-divphi-log(phi_mu_1)) +
    //   pow(vphi, 4.0)*(R::trigamma(y.at(ii)+vphi)-trvphi) +
    //   pow(vphi, 2.0)*(2.0*mu.at(ii)/phi_mu_1 - y.at(ii)) +
    //   (vphi+y.at(ii))*pow(mu.at(ii),2.0)/pow(phi_mu_1,2.0);
    
    Ibb -= (mu.at(ii)/phi_mu_1)*X.row(ii).t()*X.row(ii);
    
    // Iep -= (y.at(ii)-mu.at(ii))/pow(phi_mu_1,2.0)*dmu_keg.t();
    
    Ibe -= (1/phi_mu_1)*X.row(ii).t()*dmu_keg;
    
    // Ibp -= ((y.at(ii)-mu.at(ii))*mu.at(ii)/pow(phi_mu_1,2.0))*X.row(ii).t();
  }
  
  // printR_obj(Iee);
  // printR_obj(Ibb);
  // printR_obj(Ibp);
  // printR_obj(Ibe);
  
  /*
   * Information matrix (ASE)
   */
  
  RcppT_compite_pi(z_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau, pis);
  
  for(ii =0; ii<n_AS; ii++){
    double aa = pis.at(ii) * vtheta; 
    double bb = vtheta - aa;
    double diaa = R::digamma(aa);
    double dibb =  R::digamma(bb);
    double diaa_ni0 = R::digamma(aa + ni0.at(ii));
    double dibb_ni1 = R::digamma(bb + ni.at(ii) - ni0.at(ii));
    double triaa = R::trigamma(aa);
    double tribb =  R::trigamma(bb);
    // double triaa_ni0 = R::trigamma(aa + ni0.at(ii));
    // double tribb_ni1 = R::trigamma(bb + ni.at(ii) - ni0.at(ii));
    Expvec.zeros();
    RcppT_ASE_ExpFunc(ni.at(ii), pis.at(ii), vtheta, Expvec);
    
    dlASE_dpi  = vtheta*(diaa_ni0 - dibb_ni1 - diaa + dibb);
    d2lASE_dpi = pow(vtheta, 2.0)*(Expvec.at(2) + Expvec.at(3) - 
      triaa - tribb);
    if(z_AS.at(ii) == 0 | z_AS.at(ii) == 3){
      tmp1 = RHO_AS.at(ii)*tau.at(ii)*KAPPA + 2*(1 - RHO_AS.at(ii));
      tmp2 = RHO_AS.at(ii)*tauB.at(ii)/tmp1;
      tmp3 = RHO_AS.at(ii)*tau.at(ii)*(RHO_AS.at(ii)*tauB.at(ii)*KAPPA + 
        (1 - RHO_AS.at(ii)))/pow(tmp1, 2.0);
      dpi_dkappa = tmp2 - tmp3;
      dpi_deta   = 0;
      dpi_dgamma = 0;
      
      d2pi_dkappa = -(2.0*pow(RHO_AS.at(ii),2.0)*(1.0-RHO_AS.at(ii))*
        (2.0*tauB.at(ii) -tau.at(ii))*tau.at(ii))/pow(tmp1, 3.0);
      d2pi_deta = d2pi_dgamma  = d2pi_dkappa_eta = 
        d2pi_dkappa_gamma = d2pi_deta_gamma = 0;
      
    }else{
      tmp1 = RHO_AS.at(ii)*tauB.at(ii)*GAMMA*KAPPA + (1 - RHO_AS.at(ii))*ETA;
      tmp2 = RHO_AS.at(ii)*(tau.at(ii) - tauB.at(ii))*KAPPA +
        (1 - RHO_AS.at(ii)) + tmp1;
      tmp3 = tmp1/pow(tmp2, 2.0);
      tmp4 = 1/tmp2 -tmp3; 
      dpi_dkappa = RHO_AS.at(ii)*tauB.at(ii)*GAMMA*tmp4 -
        RHO_AS.at(ii)*(tau.at(ii) - tauB.at(ii))*tmp3;
      dpi_deta   = (1 - RHO_AS.at(ii))*tmp4;
      dpi_dgamma = RHO_AS.at(ii)*tauB.at(ii)*KAPPA*tmp4;
      
      tmp3 = tmp1/pow(tmp2, 3.0); // changed tmp3
      tmp_kappa = RHO_AS.at(ii)*(tau.at(ii) - tauB.at(ii)) + 
        RHO_AS.at(ii)* tauB.at(ii)*GAMMA;
      d2pi_dkappa = -2.0*GAMMA*RHO_AS.at(ii)*tauB.at(ii)*tmp_kappa/pow(tmp2,2.0) +
        2.0*pow(tmp_kappa,2.0)*tmp3;
      d2pi_deta   = -2.0*pow(1-RHO_AS.at(ii), 2.0)*(1/pow(tmp2, 2.0) - tmp3);
      d2pi_dgamma = -2.0*pow(RHO_AS.at(ii)*tauB.at(ii)*KAPPA,2.0)*
        (1/pow(tmp2, 2.0) - tmp3);
      d2pi_dkappa_eta  = (-RHO_AS.at(ii)*tauB.at(ii)*GAMMA*(1-RHO_AS.at(ii)) -
        tmp_kappa*(1-RHO_AS.at(ii)))/pow(tmp2, 2.0) +
        2.0*(1-RHO_AS.at(ii))*tmp_kappa*tmp3;
      d2pi_dkappa_gamma = RHO_AS.at(ii)*tauB.at(ii)/tmp2 -
        (pow(RHO_AS.at(ii)*tauB.at(ii), 2.0)*KAPPA*GAMMA + 
        tmp_kappa*RHO_AS.at(ii)*tauB.at(ii)*KAPPA)/pow(tmp2, 2.0) -
        RHO_AS.at(ii)*tauB.at(ii)*tmp3*tmp2 +
        2.0*tmp_kappa*RHO_AS.at(ii)*tauB.at(ii)*KAPPA*tmp3;
      d2pi_deta_gamma   = -2.0*KAPPA*(1-RHO_AS.at(ii))*RHO_AS.at(ii)*tauB.at(ii)/pow(tmp2,2.0) +
        2.0*KAPPA*(1.0-RHO_AS.at(ii))*RHO_AS.at(ii)*tauB.at(ii)*tmp3;
    }
    
    Iee.at(0,0) += d2lASE_dpi*pow(dpi_dkappa, 2.0)  + dlASE_dpi*d2pi_dkappa;
    Iae.at(0,1) += d2lASE_dpi*pow(dpi_deta, 2.0)    + dlASE_dpi*d2pi_deta;
    Iae.at(1,2) += d2lASE_dpi*pow(dpi_dgamma, 2.0)  + dlASE_dpi*d2pi_dgamma; 
    Iae.at(0,0) += d2lASE_dpi*dpi_dkappa*dpi_deta   + dlASE_dpi*d2pi_dkappa_eta; 
    Iae.at(1,0) += d2lASE_dpi*dpi_dkappa*dpi_dgamma + dlASE_dpi*d2pi_dkappa_gamma;
    Iae.at(0,2) += d2lASE_dpi*dpi_deta*dpi_dgamma   + dlASE_dpi*d2pi_deta_gamma; 
    
    score_alpha.at(0) += dlASE_dpi*dpi_deta;
    score_alpha.at(1) += dlASE_dpi*dpi_dgamma;
    
    W = -(R::digamma(ni.at(ii)+vtheta)+diaa*pis.at(ii) +dibb*(1.0-pis.at(ii)) -
      Expvec.at(0)*pis.at(ii)-Expvec.at(1)*(1.0-pis.at(ii)) - divtheta);
    dW_dTHETA = -pow(vtheta, 2.0)*(pow(pis.at(ii), 2.0)*(Expvec.at(2) - triaa) +
      pow(1-pis.at(ii), 2.0)*(Expvec.at(3) - tribb) +
      trvtheta - R::trigamma(vtheta+ni.at(ii)));
    Itt += 2.0*pow(vtheta,3.0)*W - pow(vtheta, 2.0)*dW_dTHETA;
    
    
    d2lASE_dTHETA_dpi = - pow(vtheta, 2.0)*(Expvec.at(0) - diaa -
      Expvec.at(1) + dibb) -
      pow(vtheta,3.0)*pis.at(ii)*(Expvec.at(2) - triaa) +
      pow(vtheta,3.0)*(1.0-pis.at(ii))*(Expvec.at(3) - tribb);
    Iet.at(0,0) += d2lASE_dTHETA_dpi*dpi_dkappa;
    Iet.at(1,0) += d2lASE_dTHETA_dpi*dpi_deta;
    Iet.at(2,0) += d2lASE_dTHETA_dpi*dpi_dgamma;
  }
  
  Iae.at(1,1) = Iae.at(0,2);
  Iaa = Iae.submat(0,1,1,2);
  
  Iee.at(0,1) += Iae.at(0,0);
  Iee.at(1,1) += Iae.at(0,1);
  Iee.at(1,2) += Iae.at(0,2);
  Iee.at(0,2) += Iae.at(1,0);
  Iee.at(2,2) += Iae.at(1,2);
  
  Iee.at(1,0) = Iee.at(0,1);
  Iee.at(2,0) = Iee.at(0,2);
  Iee.at(2,1) = Iee.at(1,2);
  
  Iat = Iet.submat(1,0,2,0);
  
  /*
   * Final information matrix
   */
  
  arma::uword np = beta_npara -1;
  M1.submat(0,0,np,np)             = Ibb;
  M1.submat(0,np+1,np,np+3)        = Ibe;
  M1.submat(np+1,0,np+3,np)        = Ibe.t();
  M1.submat(np+1,np+1,np+3,np+3)   = Iee;
  M1.submat(np+1,np+4,np+3,np+4)   = Iet;
  M1.submat(np+4,np+1,np+4,np+3)   = Iet.t();
  M1.at((np+4),(np+4))             = Itt;
  M2.rows(np+1,np+3) = Iae.t();
  M2.row(4+np)       = Iat.t();
  
  
  M1 = -M1;
  M2 = -M2;
  Iaa = -Iaa;
  
  
  OIMat.submat(0,0,(np+4),(np+4))           = M1;
  OIMat.submat(0,(np+5),(np+4),(np+6))      = M2;
  OIMat.submat((np+5),0,(np+6),(np+4))      = M2.t();
  OIMat.submat((np+5),(np+5),(np+6),(np+6)) = Iaa;
  
  // printR_obj(score_alpha);
  // printR_obj(Iee);
  // printR_obj(Iae);
  // printR_obj(Iaa);
  // printR_obj(Itt);
  // printR_obj(Iet);
  // printR_obj(OIMat);
  if(!OIMat.is_sympd()){
    Score.at(0,0) = -6.0;
    pval  = 1.0;
  }else{
    Score = score_alpha.t()*((Iaa - M2.t()*M1.i()*M2).i())*score_alpha;
    pval  = R::pchisq(Score.at(0,0),2,0,0);
  }
  
  if(std::isnan(Score.at(0,0))){
    Score.at(0,0) = -6.0;
    pval = 1.0;
  }
  return Rcpp::List::create(
    Rcpp::Named("Score", Score.at(0,0)),
    Rcpp::Named("pval", pval)
  );
}

/* ---------------------------
 * SNP-GENE pair
 ---------------------------*/

// [[Rcpp::export]]
void RcppT_trecase_mtest(const arma::mat& Y, const arma::mat& Y1,
                         const arma::mat& Y2, const arma::mat& Z,
                         const arma::mat& XX, const arma::vec& RHO,
                         const arma::mat& CNV1, const arma::mat& CNV2,
                         const arma::vec& SNP_pos,
                         const arma::uvec& sChr,
                         const arma::vec& gene_start,
                         const arma::vec& gene_end,
                         const arma::uvec& gChr,
                         const List& GeneSnpList,
                         const char* file_trec = "trecT.txt",
                         const char* file_trecase = "trecaseT.txt",
                         const bool& useLRT = false,
                         const double& transTestP = 0.01,
                         const double& cis_window=1e5,
                         const bool& useASE = 1, const int& min_ASE_total=8,
                         const int& min_nASE= 5, const int& min_nASE_het=5,
                         const double& eps=5e-5,
                         const arma::uword& max_iter=400L,
                         const bool& show=false){
  arma::uword gg, ii, xi, h1, h0, z0;
  arma::uword ssBegin = 0, ss = 0;
  double nSam = Y.n_rows;
  double ctcode = -1.0;
  double CisTrans_Chisq, CisTrans_Pval, log_theta, phi;
  int converge =0;
  // double pp = XX.n_cols+1;;
  arma::vec _y      = arma::zeros<arma::vec>(nSam);
  arma::vec y       = arma::zeros<arma::vec>(nSam);
  arma::vec z       = arma::zeros<arma::vec>(nSam);
  
  arma::vec y1      = arma::zeros<arma::vec>(nSam);
  arma::vec y2      = arma::zeros<arma::vec>(nSam);
  arma::vec RHO1    = arma::zeros<arma::vec>(nSam);
  
  arma::vec RHO_AS  = arma::zeros<arma::vec>(nSam);
  
  arma::vec _tau1   = arma::zeros<arma::vec>(nSam);
  arma::vec _tau2   = arma::zeros<arma::vec>(nSam);
  
  arma::vec tau1    = arma::zeros<arma::vec>(nSam);
  arma::vec tau2    = arma::zeros<arma::vec>(nSam);
  arma::vec tauB    = arma::zeros<arma::vec>(nSam);
  arma::vec tau     = arma::zeros<arma::vec>(nSam);
  arma::mat X       = arma::zeros<arma::mat>(nSam, XX.n_cols);
  
  arma::vec ni   = arma::zeros<arma::vec>(nSam);
  arma::vec lbc  = arma::zeros<arma::vec>(nSam);
  arma::vec ni0  = arma::zeros<arma::vec>(nSam);
  arma::vec z_AS = arma::zeros<arma::vec>(nSam);
  
  // arma::vec offsets = arma::zeros<arma::vec>(nSam);
  
  Rcpp::List res_trec, res_trecase, CT_score, res_ase, res_trec_ase;
  arma::vec PAR = arma::zeros<arma::vec>(3);
  arma::vec reg_par = arma::zeros<arma::vec>(XX.n_cols);
  
  //create files for TReC and TReCASE results
  FILE * f1, * f2;
  
  if(useASE){
    
    f2 = fopen(file_trecase, "w");
    
    fprintf(f2, "GeneRowID\tMarkerRowID\tTReCASE_kappa\tTReCASE_eta\tTReCASE_gamma\t");
    fprintf(f2, "TReCASE_LL.full\tTReCASE_pEta\tTReCASE_pGamma\t");
    for(xi=0;xi<XX.n_cols;xi++){
      fprintf(f2, "TReCASE_beta%d\t", xi);
    }
    fprintf(f2, "TReCASE_phi\tTReCASE_theta\t");
    fprintf(f2, "Converge\tCisTrans_Chisq\t");
    fprintf(f2, "CisTrans_Pvalue\tnSam\tnHet\n");
  }else{
    f1 = fopen(file_trec, "w");
    fprintf(f1, "GeneRowID\tMarkerRowID\tTReC_kappa\tTReC_eta\tTReC_gamma\t");
    fprintf(f1, "TReC_LL.full\tTReC_pEta\tTReC_pGamma\tTReC_Conv\t");
    for(xi=0;xi<XX.n_cols;xi++){
      fprintf(f1, "beta%d\t", xi);
    }
    fprintf(f1, "phi\tnSam\n");
  }
  
  if(GeneSnpList.length() > 0){
    
    for(gg=0; gg<GeneSnpList.length(); gg++){
      
      if(gg % 100 == 0){
        Rprintf("Begin analysis for Gene %d \n", gg+1);
      }
      
      if(GeneSnpList[gg]==R_NilValue){
        continue;
      }
      _y    = Y.col(gg);
      _tau1 = CNV1.col(gg);
      _tau2 = CNV2.col(gg);
      if(useASE){
        y1 = Y1.col(gg);
        y2 = Y2.col(gg);
      }
      
      //loop through the matrix to exclude NA value
      //organize tau1 tau2 tauB tau 
      arma::vec ssVec = GeneSnpList[gg]; 
      
      for(ssBegin = 0; ssBegin < ssVec.n_elem; ssBegin++){
        ss = ssVec.at(ssBegin)-1;
        // Rprintf("ss %d \n", ss+1);
        
        arma::vec zz2 = Z.col(ss);
        h1 = 0, h0 = 0, z0 = 0;
        X.zeros();
        z.zeros();
        y.zeros();
        z_AS.zeros();
        ni.zeros();
        ni0.zeros();
        tau1.zeros();
        tau2.zeros();
        tauB.zeros();
        tau.zeros();
        RHO1.zeros();
        RHO_AS.zeros();
        lbc.zeros();
        
        for(ii=0;ii<nSam;ii++){
          
          if(zz2.at(ii) != -9 & _tau1.at(ii) != -9 & _tau2.at(ii) != -9){
            z.at(z0)    = zz2.at(ii);
            X.row(z0)   = XX.row(ii);
            y.at(z0)    = _y.at(ii);
            tau1.at(z0) = _tau1.at(ii);
            tau2.at(z0) = _tau2.at(ii);
            RHO1.at(z0) = RHO.at(ii);
            // if(zz2.at(ii)==2){
            //   z.at(z0) = 1;
            // }else if(zz2.at(ii)==3){
            //   z.at(z0) = 2;
            // }
            // 
            if(useASE){
              if(y1.at(ii) + y2.at(ii) >= min_ASE_total){
                
                z_AS.at(h0) = zz2.at(ii);
                ni.at(h0)   = y1.at(ii) + y2.at(ii);
                tau.at(h0)  = _tau1.at(ii) + _tau2.at(ii);
                RHO_AS.at(h0) = RHO.at(ii);
                
                if(zz2.at(ii)==1){
                  ni0.at(h0) = y2.at(ii);
                  tauB.at(h0) = _tau2.at(ii); 
                  h1++;
                }else if(zz2.at(ii)==2){
                  ni0.at(h0) = y1.at(ii);
                  tauB.at(h0) = _tau1.at(ii);
                  h1++;
                }else if(zz2.at(ii)==0){
                  tauB.at(h0) = _tau2.at(ii); 
                  ni0.at(h0)  = y2.at(ii);
                }else if(zz2.at(ii)==3){
                  tauB.at(h0) = _tau1.at(ii);
                  ni0.at(h0)  = y1.at(ii);
                }
                lbc.at(h0) = R::lchoose(ni.at(h0), ni0.at(h0));
                h0++;
              }
            }
            z0++;
          }
        }        
        arma::vec lgy1 = Rcpp_lgy_add_1(y.subvec(0, z0-1)); //lgamma(y + 1)
        
        // begin tracase
        if(useASE & h1 >= min_nASE_het & h0 >= min_nASE){
          res_trecase = RcppT_trecase(y.subvec(0, z0-1), z.subvec(0, z0-1), 
                                      z_AS.subvec(0, h0-1), RHO1.subvec(0, z0-1), 
                                      RHO_AS.subvec(0, h0-1), X.rows(0, z0-1),
                                      tau1.subvec(0, z0-1), tau2.subvec(0, z0-1), 
                                      lgy1,
                                      ni0.subvec(0, h0-1), ni.subvec(0, h0-1), 
                                      tauB.subvec(0, h0-1), tau.subvec(0, h0-1), 
                                      lbc.subvec(0, h0-1),
                                      max_iter, eps, show);
          PAR = as<arma::vec>(res_trecase["PAR"]); 
          reg_par = as<arma::vec>(res_trecase["reg_par"]);
          arma::vec BETA = reg_par.subvec(0,X.n_cols-1);
          phi = exp(reg_par.at(X.n_cols));
          log_theta = res_trecase["log_theta"];
          converge = res_trecase["converge"];
          // Rprintf("phi %.2f\n",phi);
          // Rprintf("log_theta %.2f\n",log_theta);
          // Rprintf("trecase converge %d\n",converge);
          
          if(converge){
            if(!useLRT){
              // Rprintf("OBS");
              CT_score = RcppT_CisTrans_ScoreObs(PAR, y.subvec(0, z0-1),
                                                 z.subvec(0, z0-1), 
                                                 z_AS.subvec(0, h0-1), RHO1.subvec(0, z0-1), 
                                                 RHO_AS.subvec(0, h0-1), X.rows(0, z0-1), 
                                                 BETA, phi,
                                                 tau1.subvec(0, z0-1), tau2.subvec(0, z0-1),
                                                 lgy1,
                                                 ni0.subvec(0, h0-1), ni.subvec(0, h0-1), 
                                                 log_theta,
                                                 tauB.subvec(0, h0-1), tau.subvec(0, h0-1), 
                                                 lbc.subvec(0, h0-1));
              CisTrans_Chisq = as<double>(CT_score["Score"]);
              CisTrans_Pval = as<double>(CT_score["pval"]);
              
              if(CisTrans_Chisq < 0.0){
                // Rprintf("EXP");
                CT_score = RcppT_CisTrans_Score(PAR, y.subvec(0, z0-1),
                                                z.subvec(0, z0-1), 
                                                z_AS.subvec(0, h0-1), RHO1.subvec(0, z0-1), 
                                                RHO_AS.subvec(0, h0-1), X.rows(0, z0-1), 
                                                BETA, phi,
                                                tau1.subvec(0, z0-1), tau2.subvec(0, z0-1),
                                                lgy1,
                                                ni0.subvec(0, h0-1), ni.subvec(0, h0-1), 
                                                log_theta,
                                                tauB.subvec(0, h0-1), tau.subvec(0, h0-1), 
                                                lbc.subvec(0, h0-1));
                
              }
              CisTrans_Chisq = as<double>(CT_score["Score"]);
              CisTrans_Pval = as<double>(CT_score["pval"]);
              ctcode = CisTrans_Chisq; 
              if(CisTrans_Pval < transTestP){
                ctcode = -2.0;
              }
              // Rprintf("pval %.8f\n",CisTrans_Pval);
              // Rprintf("ct code %.2f\n",ctcode);
              
            }
            
            if(useLRT | CisTrans_Chisq < 0){ 
              // Rprintf("lrt\n");
              arma::vec para0 = arma::zeros<arma::vec>(5);
              res_trec_ase = RcppT_trec_ase(para0, y.subvec(0, z0-1), z.subvec(0, z0-1),
                                            z_AS.subvec(0, h0-1), RHO1.subvec(0, z0-1),
                                            RHO_AS.subvec(0, h0-1), X.rows(0, z0-1),
                                            tau1.subvec(0, z0-1), tau2.subvec(0, z0-1), 
                                            lgy1,
                                            ni0.subvec(0, h0-1), ni.subvec(0, h0-1), 
                                            tauB.subvec(0, h0-1), tau.subvec(0, h0-1), 
                                            lbc.subvec(0, h0-1), max_iter, eps, show);
              if(as<bool>(res_trec_ase["converge"])){
                CisTrans_Chisq = -2.0*(as<double>(res_trecase["LL"]) - as<double>(res_trec_ase["LL"]));
                CisTrans_Pval =  R::pchisq(CisTrans_Chisq,2,0,0);
                ctcode = CisTrans_Chisq; 
              }else{
                CisTrans_Chisq  = -6.0;
                CisTrans_Pval =  1.0;
                ctcode = -3.0;
              }
              
              // Rprintf("lrt pval %.8f\n",CisTrans_Pval);
              
              if(CisTrans_Pval < transTestP){
                ctcode = -2.0;
              }
            }
          }else{
            CisTrans_Chisq  = -5.0;
            CisTrans_Pval  = R_NaN;
          }
          
        }
        
        if(useASE == 0 | ctcode < 0.0 | 
           h1 < min_nASE_het | h0 < min_nASE | !converge){
          // Rprintf("trec\n");
          
          res_trec = RcppT_trec(y.subvec(0, z0-1), z.subvec(0, z0-1), 
                                RHO1.subvec(0, z0-1), X.rows(0, z0-1),
                                tau1.subvec(0, z0-1), tau2.subvec(0, z0-1), lgy1,
                                max_iter, eps, show); 
          PAR = as<arma::vec>(res_trec["PAR"]); 
          reg_par = as<arma::vec>(res_trec["reg_par"]);
          converge = res_trec["converge"];
          // Rprintf("trec converge %d\n",converge);
          if(useASE){
            
            // fprintf(f2, "GeneRowID\tMarkerRowID\tTReCASE_kappa\tTReCASE_eta\tTReCASE_gamma\t");
            fprintf(f2, "%d\t%d\t%.2e\t%.2e\t%.2e\t",
                    gg+1,ss+1,exp(PAR.at(0)),
                    exp(PAR.at(1)), exp(PAR.at(2)));
            // fprintf(f2, "TReCASE_LL.full\tTReCASE_pEta\tTReCASE_pGamma\t");
            fprintf(f2, "%.2e\t%.4e\t%.4e\t", 
                    as<double>(res_trec["LL"]),
                    as<double>(res_trec["p_eta"]), 
                    as<double>(res_trec["p_gamma"]));
            for(xi=0;xi<XX.n_cols;xi++){
              fprintf(f2, "%.2e\t", reg_par.at(xi));
            }
            // fprintf(f2, "TReCASE_phi\tTReCASE_theta\t");
            fprintf(f2, "%.2e\tNA\t",exp(reg_par.at(XX.n_cols)));
            // fprintf(f2, "Converge\tCisTrans_Chisq\t");
            fprintf(f2, "%d\t%.2f\t", converge, CisTrans_Chisq);
            // fprintf(f2, "CisTrans_Pvalue\tnSam\tnHetn");
            fprintf(f2, "%.4e\t%d\t%d\n",CisTrans_Pval, z0, h0);
            
          }else{
            // fprintf(f1, "GeneRowID\tMarkerRowID\tTReC_kappa\tTReC_eta\tTReC_gamma\t");
            fprintf(f1, "%d\t%d\t%.2e\t%.2e\t%.2e\t",
                    gg+1,ss+1,exp(PAR.at(0)),
                    exp(PAR.at(1)), exp(PAR.at(2)));
            fprintf(f1, "%.2e\t%.4e\t%.4e\t%d\t", 
                    as<double>(res_trec["LL"]),
                    as<double>(res_trec["p_eta"]), 
                    as<double>(res_trec["p_gamma"]),
                    converge);
            for(xi=0;xi<XX.n_cols;xi++){
              fprintf(f1, "%.2e\t", reg_par.at(xi));
            }
            fprintf(f1, "%.2e\t%d\n", exp(reg_par.at(xi)), z0);
          }
          
          
        }else{
          
          // fprintf(f2, "GeneRowID\tMarkerRowID\tTReCASE_kappa\tTReCASE_eta\tTReCASE_gamma\t");
          fprintf(f2, "%d\t%d\t%.2e\t%.2e\t%.2e\t",
                  gg+1,ss+1,exp(PAR.at(0)),
                  exp(PAR.at(1)), exp(PAR.at(2)));
          // fprintf(f2, "TReCASE_LL.full\tTReCASE_pEta\tTReCASE_pGamma\t");
          fprintf(f2, "%.2e\t%.4e\t%.4e\t", 
                  as<double>(res_trecase["LL"]),
                  as<double>(res_trecase["p_eta"]), 
                  as<double>(res_trecase["p_gamma"]));
          for(xi=0;xi<XX.n_cols;xi++){
            fprintf(f2, "%.2e\t", reg_par.at(xi));
          }
          // fprintf(f2, "TReCASE_phi\tTReCASE_theta\t");
          fprintf(f2, "%.2e\t%.2e\t",phi, 
                  exp(log_theta));
          // fprintf(f2, "Converge\tCisTrans_Chisq\t");
          fprintf(f2, "%d\t%.2f\t", 
                  converge,
                  CisTrans_Chisq);
          // fprintf(f2, "CisTrans_Pvalue\tnSam\tnHetn");
          fprintf(f2, "%.4e\t%d\t%d\n",CisTrans_Pval, z0, h0);
          
        }
        
        
        
      }
    }
    
  }else{
    for(gg=0; gg<Y.n_cols; gg++){
      
      if(gg % 100 == 0){
        Rprintf("Begin analysis for Gene %d \n", gg+1);
      }
      
      _y    = Y.col(gg);
      _tau1 = CNV1.col(gg);
      _tau2 = CNV2.col(gg);
      if(useASE){
        y1 = Y1.col(gg);
        y2 = Y2.col(gg);
      }
      
      //loop through the matrix to exclude NA value
      //organize tau1 tau2 tauB tau
      
      for(ss = ssBegin; ss < Z.n_cols; ss++){
        
        if(gChr.at(gg) != sChr.at(ss)){
          ssBegin = ss;
          break;
        }
        
        if(ss % 5000 == 0 & show){
          Rprintf("Begin analysis for SNP %d  \n", ss+1);
        }
        
        if(SNP_pos.at(ss) > gene_start.at(gg) - cis_window &&
           SNP_pos.at(ss) < gene_end.at(gg)   + cis_window){
          
          arma::vec zz2 = Z.col(ss);
          h1 = 0, h0 = 0, z0 = 0;
          X.zeros();
          z.zeros();
          y.zeros();
          z_AS.zeros();
          ni.zeros();
          ni0.zeros();
          tau1.zeros();
          tau2.zeros();
          tauB.zeros();
          tau.zeros();
          RHO1.zeros();
          RHO_AS.zeros();
          lbc.zeros();
          
          for(ii=0;ii<nSam;ii++){
            
            if(zz2.at(ii) != -9 & _tau1.at(ii) != -9 & _tau2.at(ii) != -9){
              z.at(z0)    = zz2.at(ii);
              X.row(z0)   = XX.row(ii);
              y.at(z0)    = _y.at(ii);
              tau1.at(z0) = _tau1.at(ii);
              tau2.at(z0) = _tau2.at(ii);
              RHO1.at(z0) = RHO.at(ii);
              // if(zz2.at(ii)==2){
              //   z.at(z0) = 1;
              // }else if(zz2.at(ii)==3){
              //   z.at(z0) = 2;
              // }
              
              if(useASE){
                if(y1.at(ii) + y2.at(ii) >= min_ASE_total){
                  
                  z_AS.at(h0) = zz2.at(ii);
                  ni.at(h0)   = y1.at(ii) + y2.at(ii);
                  tau.at(h0)  = _tau1.at(ii) + _tau2.at(ii);
                  RHO_AS.at(h0) = RHO.at(ii);
                  
                  if(zz2.at(ii)==1){
                    ni0.at(h0) = y2.at(ii);
                    tauB.at(h0) = _tau2.at(ii);
                    h1++;
                  }else if(zz2.at(ii)==2){
                    ni0.at(h0) = y1.at(ii);
                    tauB.at(h0) = _tau1.at(ii);
                    h1++;
                  }else if(zz2.at(ii)==0){
                    tauB.at(h0) = _tau2.at(ii);
                    ni0.at(h0)  = y2.at(ii);
                  }else if(zz2.at(ii)==3){
                    tauB.at(h0) = _tau1.at(ii);
                    ni0.at(h0)  = y1.at(ii);
                  }
                  lbc.at(h0) = R::lchoose(ni.at(h0), ni0.at(h0));
                  h0++;
                }
              }
              z0++;
            }
          }
          arma::vec lgy1 = Rcpp_lgy_add_1(y.subvec(0, z0-1)); //lgamma(y + 1)
          
          // begin tracase
          if(useASE & h1 >= min_nASE_het & h0 >= min_nASE){
            res_trecase = RcppT_trecase(y.subvec(0, z0-1), z.subvec(0, z0-1),
                                        z_AS.subvec(0, h0-1), RHO1.subvec(0, z0-1),
                                        RHO_AS.subvec(0, h0-1), X.rows(0, z0-1),
                                        tau1.subvec(0, z0-1), tau2.subvec(0, z0-1),
                                        lgy1,
                                        ni0.subvec(0, h0-1), ni.subvec(0, h0-1),
                                        tauB.subvec(0, h0-1), tau.subvec(0, h0-1),
                                        lbc.subvec(0, h0-1),
                                        max_iter, eps, show);
            PAR = as<arma::vec>(res_trecase["PAR"]);
            reg_par = as<arma::vec>(res_trecase["reg_par"]);
            arma::vec BETA = reg_par.subvec(0,X.n_cols-1);
            phi = exp(reg_par.at(X.n_cols));
            log_theta = res_trecase["log_theta"];
            converge = res_trecase["converge"];
            // Rprintf("phi %.2f\n",phi);
            // Rprintf("log_theta %.2f\n",log_theta);
            // Rprintf("trecase converge %d\n",converge);
            
            if(converge){
              if(!useLRT){
                // Rprintf("OBS");
                CT_score = RcppT_CisTrans_ScoreObs(PAR, y.subvec(0, z0-1),
                                                   z.subvec(0, z0-1),
                                                   z_AS.subvec(0, h0-1), RHO1.subvec(0, z0-1),
                                                   RHO_AS.subvec(0, h0-1), X.rows(0, z0-1),
                                                   BETA, phi,
                                                   tau1.subvec(0, z0-1), tau2.subvec(0, z0-1),
                                                   lgy1,
                                                   ni0.subvec(0, h0-1), ni.subvec(0, h0-1),
                                                   log_theta,
                                                   tauB.subvec(0, h0-1), tau.subvec(0, h0-1),
                                                   lbc.subvec(0, h0-1));
                CisTrans_Chisq = as<double>(CT_score["Score"]);
                CisTrans_Pval = as<double>(CT_score["pval"]);
                
                if(CisTrans_Chisq < 0.0){
                  // Rprintf("EXP");
                  CT_score = RcppT_CisTrans_Score(PAR, y.subvec(0, z0-1),
                                                  z.subvec(0, z0-1),
                                                  z_AS.subvec(0, h0-1), RHO1.subvec(0, z0-1),
                                                  RHO_AS.subvec(0, h0-1), X.rows(0, z0-1),
                                                  BETA, phi,
                                                  tau1.subvec(0, z0-1), tau2.subvec(0, z0-1),
                                                  lgy1,
                                                  ni0.subvec(0, h0-1), ni.subvec(0, h0-1),
                                                  log_theta,
                                                  tauB.subvec(0, h0-1), tau.subvec(0, h0-1),
                                                  lbc.subvec(0, h0-1));
                  
                  
                }
                CisTrans_Chisq = as<double>(CT_score["Score"]);
                CisTrans_Pval = as<double>(CT_score["pval"]);
                ctcode = CisTrans_Chisq;
                if(CisTrans_Pval < transTestP){
                  ctcode = -2.0;
                }
                // Rprintf("pval %.8f\n",CisTrans_Pval);
                // Rprintf("ct code %.2f\n",ctcode);
                
              }
              
              if(useLRT | CisTrans_Chisq < 0){
                // Rprintf("lrt\n");
                arma::vec para0 = arma::zeros<arma::vec>(5);
                res_trec_ase = RcppT_trec_ase(para0, y.subvec(0, z0-1), z.subvec(0, z0-1),
                                              z_AS.subvec(0, h0-1), RHO1.subvec(0, z0-1),
                                              RHO_AS.subvec(0, h0-1), X.rows(0, z0-1),
                                              tau1.subvec(0, z0-1), tau2.subvec(0, z0-1),
                                              lgy1,
                                              ni0.subvec(0, h0-1), ni.subvec(0, h0-1),
                                              tauB.subvec(0, h0-1), tau.subvec(0, h0-1),
                                              lbc.subvec(0, h0-1), max_iter, eps, show);
                if(as<bool>(res_trec_ase["converge"])){
                  CisTrans_Chisq = -2.0*(as<double>(res_trecase["LL"]) - as<double>(res_trec_ase["LL"]));
                  CisTrans_Pval =  R::pchisq(CisTrans_Chisq,2,0,0);
                  ctcode = CisTrans_Chisq;
                }else{
                  CisTrans_Chisq  = -6.0;
                  CisTrans_Pval =  1.0;
                  ctcode = -3.0;
                }
                
                // Rprintf("lrt pval %.8f\n",CisTrans_Pval);
                
                if(CisTrans_Pval < transTestP){
                  ctcode = -2.0;
                }
              }
            }else{
              CisTrans_Chisq  = -5.0;
              CisTrans_Pval  = R_NaN;
            }
            
          }
          
          if(useASE == 0 | ctcode < 0.0 |
             h1 < min_nASE_het | h0 < min_nASE | !converge){
            // Rprintf("trec\n");
            
            res_trec = RcppT_trec(y.subvec(0, z0-1), z.subvec(0, z0-1),
                                  RHO1.subvec(0, z0-1), X.rows(0, z0-1),
                                  tau1.subvec(0, z0-1), tau2.subvec(0, z0-1), lgy1,
                                  max_iter, eps, show);
            PAR = as<arma::vec>(res_trec["PAR"]);
            reg_par = as<arma::vec>(res_trec["reg_par"]);
            converge = res_trec["converge"];
            // Rprintf("trec converge %d\n",converge);
            if(useASE){
              
              // fprintf(f2, "GeneRowID\tMarkerRowID\tTReCASE_kappa\tTReCASE_eta\tTReCASE_gamma\t");
              fprintf(f2, "%d\t%d\t%.2e\t%.2e\t%.2e\t",
                      gg+1,ss+1,exp(PAR.at(0)),
                      exp(PAR.at(1)), exp(PAR.at(2)));
              // fprintf(f2, "TReCASE_LL.full\tTReCASE_pEta\tTReCASE_pGamma\t");
              fprintf(f2, "%.2e\t%.4e\t%.4e\t",
                      as<double>(res_trec["LL"]),
                      as<double>(res_trec["p_eta"]),
                      as<double>(res_trec["p_gamma"]));
              for(xi=0;xi<XX.n_cols;xi++){
                fprintf(f2, "%.2e\t", reg_par.at(xi));
              }
              // fprintf(f2, "TReCASE_phi\tTReCASE_theta\t");
              fprintf(f2, "%.2e\tNA\t",exp(reg_par.at(XX.n_cols)));
              // fprintf(f2, "Converge\tCisTrans_Chisq\t");
              fprintf(f2, "%d\t%.2f\t", converge, CisTrans_Chisq);
              // fprintf(f2, "CisTrans_Pvalue\tnSam\tnHetn");
              fprintf(f2, "%.4e\t%d\t%d\n",CisTrans_Pval, z0, h0);
              
            }else{
              // fprintf(f1, "GeneRowID\tMarkerRowID\tTReC_kappa\tTReC_eta\tTReC_gamma\t");
              fprintf(f1, "%d\t%d\t%.2e\t%.2e\t%.2e\t",
                      gg+1,ss+1,exp(PAR.at(0)),
                      exp(PAR.at(1)), exp(PAR.at(2)));
              fprintf(f1, "%.2e\t%.4e\t%.4e\t%d\t",
                      as<double>(res_trec["LL"]),
                      as<double>(res_trec["p_eta"]),
                      as<double>(res_trec["p_gamma"]),
                      converge);
              for(xi=0;xi<XX.n_cols;xi++){
                fprintf(f1, "%.2e\t", reg_par.at(xi));
              }
              fprintf(f1, "%.2e\t%d\n", exp(reg_par.at(xi)), z0);
            }
            
            
          }else{
            
            // fprintf(f2, "GeneRowID\tMarkerRowID\tTReCASE_kappa\tTReCASE_eta\tTReCASE_gamma\t");
            fprintf(f2, "%d\t%d\t%.2e\t%.2e\t%.2e\t",
                    gg+1,ss+1,exp(PAR.at(0)),
                    exp(PAR.at(1)), exp(PAR.at(2)));
            // fprintf(f2, "TReCASE_LL.full\tTReCASE_pEta\tTReCASE_pGamma\t");
            fprintf(f2, "%.2e\t%.4e\t%.4e\t",
                    as<double>(res_trecase["LL"]),
                    as<double>(res_trecase["p_eta"]),
                    as<double>(res_trecase["p_gamma"]));
            for(xi=0;xi<XX.n_cols;xi++){
              fprintf(f2, "%.2e\t", reg_par.at(xi));
            }
            // fprintf(f2, "TReCASE_phi\tTReCASE_theta\t");
            fprintf(f2, "%.2e\t%.2e\t",phi,
                    exp(log_theta));
            // fprintf(f2, "Converge\tCisTrans_Chisq\t");
            fprintf(f2, "%d\t%.2f\t",
                    converge,
                    CisTrans_Chisq);
            // fprintf(f2, "CisTrans_Pvalue\tnSam\tnHetn");
            fprintf(f2, "%.4e\t%d\t%d\n",CisTrans_Pval, z0, h0);
            
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

