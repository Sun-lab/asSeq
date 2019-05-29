#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <stdio.h>


// [[Rcpp::depends("RcppArmadillo")]]

template<typename T>
void printR_obj(const T& obj){
  Rcpp::Rcout << obj << std::endl;
}

const double LOG2PI = log(2*arma::datum::pi);

// --------------------
// Intermediate Functions
// --------------------

// [[Rcpp::export]]
double Rcpp_norm(const arma::vec& a){
  return arma::norm(a);
}

// [[Rcpp::export]]
double Rcpp_logSumExp(const arma::vec& log_x){

  if( log_x.n_elem == 1 ){
    return log_x.at(0);
  } else {
    double max_val = max(log_x);

    arma::vec log_x_2 = log_x - max_val;

    return log(arma::sum(arma::exp(log_x_2))) + max_val;
  }
}

// [[Rcpp::export]]
double Rcpp_min_diff(const arma::vec& x){
  return arma::min(arma::diff(arma::sort(x)));
}

// [[Rcpp::export]]
arma::vec Rcpp_lgy_add_1(const arma::vec& y){
  return arma::lgamma(y + 1);
}

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
 * NB regression
 */
// [[Rcpp::export]]
double RcppT_reg_LL(const arma::vec& y, const arma::mat& X,
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
double RcppT_loglikNB(const arma::vec& y, const double& phi, 
                     const arma::vec& lgy1, const arma::vec& mu){
  //lgy1 = std::lgamma(y + 1)
  double loglik = 0.0;
  arma::uword ii;
  double vphi = 1.0/phi;
  double lgvphi = lgamma(vphi);
  
  for(ii = 0; ii<lgy1.size(); ii++){
    if(y.at(ii) > 0){
      loglik += lgamma(y.at(ii) + vphi) -lgvphi - lgy1.at(ii) +
        y.at(ii) * std::log(mu.at(ii));
    }
    loglik += vphi*std::log(vphi) -
      (vphi+y.at(ii)) * std::log(vphi+mu.at(ii));
  }
  
  return loglik;
}

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
  
  old_LL = RcppT_loglikNB(y, phi, lgy1, mu);
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
  double new_LL, curr_LL, LL0;
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
      printR_obj(new_LL);
      printR_obj(curr_LL);
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
  double lrt = (new_LL-LL0)*2.0;
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


// ----------------------------------------------------------------------
// TRECASE
// ----------------------------------------------------------------------

// [[Rcpp::export]]
double RcppT_TReCASE_LL_KEG(const arma::vec& para, const int& H0,
                        const arma::vec& y, const arma::vec& z,
                        const arma::vec& z_AS, const arma::vec& RHO,
                        const arma::vec& RHO_AS, const arma::vec& tau1,
                        const arma::vec& tau2, const arma::vec& ni0,
                        const arma::vec& ni, const double& THETA,
                        const arma::vec& tauB, const arma::vec& tau,
                        const arma::vec& lgy1, const arma::vec& lbc,
                        arma::vec& pis, arma::vec& mu){
  double LL_NB, LL_BB;
  double log_theta = log(THETA);
  LL_NB = RcppT_loglikNB_KEG(para, H0, y, z, phi, RHO, tau1, 
                             tau2, lgy1, expXbeta, offsets, mu);
  LL_BB = RcppT_loglikBB_KEG(para, H0, z_AS, RHO_AS, ni, ni0, log_theta,
                            lbc, tauB, tau, pis);
  return LL_NB + LL_BB;
}

double RcppT_TReCASE_LL(const arma::vec& y, ){
  
  LL_NB = RcppT_loglikNB(y, phi, lgy1, mu) 
  LL_BB = RcppT_loglikBB_THETA(ni, ni0, log_theta, pis, lbc);

}
