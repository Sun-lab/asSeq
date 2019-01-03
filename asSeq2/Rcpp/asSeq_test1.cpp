#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>

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
double Rcpp_loglikNB(const double& phi,const arma::vec& mu1,
                const arma::vec& y, const arma::vec& lgy1){
  //lgy1 = std::lgamma(y + 1)
  double loglik = 0.0;
  arma::uword ii;
  double vphi = 1.0/phi;
  
  for(ii = 0; ii< y.n_elem; ii++){
    if(y.at(ii) > 0){
      loglik += lgamma(y.at(ii) + vphi) -lgamma(vphi) - lgy1.at(ii) +
        y.at(ii) * std::log(mu1.at(ii));
    }
    loglik += vphi*std::log(vphi) - 
      (vphi+y.at(ii)) * std::log(vphi+mu1.at(ii));
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
    loglik += y.at(ii)*std::log(mu1.at(ii)) - mu1.at(ii) - lgy1.at(ii);
    //loglik += R::dpois(y.at(ii), mu1.at(ii), true);
  }
  return loglik;
}

// [[Rcpp::export]] 
double Rcpp_logLTReC(const double& bxj, const arma::vec& y, 
                     const arma::mat& X, const arma::vec& z,
                     const arma::vec& BETA, const double& phi, const bool& fam_nb,
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
    return Rcpp_loglikNB(phi, mu, y, lgy1);
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
  double grad = 0.0, df_dmu =0.0, dmu_db = 0.0, df_dmu2 =0.0;
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
Rcpp::List Rcpp_trec_bxj(const arma::vec& y, const arma::mat& X,
                         double& bxj, const arma::vec& z,
                         const arma::vec& BETA, double& phi,const bool& fam_nb,
                         const arma::vec& lgy1, const arma::uword& max_iter = 4e3,
                         const double& eps = 1e-7,const bool& show = true){
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
    if(show){
      if(uu == 0){
        printR_obj("Failed update");
      } else if(uu == 1){
        printR_obj("Newton-Raphson update");
      } else {
        printR_obj("Gradient-Descent update");
      }
    }
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
    Rcpp::Named("bxj", old_bxj)
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
  
  for(ii = 0; ii< y.n_elem; ii++){
    mu.at(ii) = std::exp( arma::dot(X.row(ii).t(),BETA) + offsets.at(ii) );
    if(y.at(ii) > 0){
      LL += lgamma(y.at(ii) + vphi) -lgamma(vphi) - lgy1.at(ii) +
        y.at(ii) * std::log(mu.at(ii));
    }
    LL += vphi*std::log(vphi) - 
      (vphi+y.at(ii)) * std::log(vphi+mu.at(ii));
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

  for(ii = 0; ii < y.n_elem; ii++){
    //mu = std::exp( arma::dot(X.row(ii).t(),BETA) + offsets.at(ii) );
    phi_mu1 = (1.0 + phi*mu.at(ii));
    // Part BETA
    GRAD.subvec(0,pp - 1) += (y.at(ii) - mu.at(ii))/phi_mu1 * X.row(ii).t();
    // Part PHI*
    GRAD.at(pp) += (-1.0*vphi*vphi)*(R::digamma(y.at(ii) + vphi) - digam_vphi -
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
  arma::uword ii;
  arma::mat HESS = arma::zeros<arma::mat>(X.n_cols, X.n_cols);
  
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
  arma::uword ii;
  arma::mat HESS = arma::zeros<arma::mat>(X.n_cols, X.n_cols);
  
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
    
    if(show){
      if(uu == 0){
        printR_obj("Failed update");
      } else if(uu == 1){
        printR_obj("Newton-Raphson update");
      } else {
        printR_obj("Gradient-Descent update");
      }
    }
    
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
    for(jj=0; jj<30; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4, jj);
      new_xk    = xk + tmp_alpha * p_k;
      new_LL    = fnscale * Rcpp_reg_LL(y, X, offsets, new_xk, fam_nb, lgy1, mu);
      
      if(new_LL < old_LL){ //minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * Rcpp_reg_grad(y, X, mu, new_xk,fam_nb) - gr_k;
        ys  = arma::dot(y_k, s_k);
        
        if(ys > 0.0){
          if(show) printR_obj("Update xk and inv_Bk");
          ISYT   = I_num_params - (s_k * y_k.t()) /ys;
          inv_Bk = ISYT * inv_Bk * ISYT.t() + s_k * s_k.t() / ys;
        }else{
          if(show) printR_obj("Update xk only");
        }
        xk = new_xk; 
        old_LL = new_LL;
        uu = 1;
        break;
      }
    }
    
    if(uu==0){
      if(Rcpp_norm(gr_k) > 1.0){
        if(show) printR_obj("Reset inv_Bk");
        inv_Bk = I_num_params;
      }else{
        if(show) printR_obj("Failed in search");
        break;
      }
    }
    
    //check convergence 
    if(iter > 0){
      if(std::abs(curr_LL - old_LL) < eps && 
         Rcpp_norm(curr_xk - xk) < eps){
        gr_k = Rcpp_reg_grad(y, X, mu, xk, fam_nb);
        if(Rcpp_norm(gr_k) < eps){
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


// [[Rcpp::export]]
Rcpp::List Rcpp_trec(const arma::vec& y, const arma::mat& X,  
                     const arma::vec& z, const bool& fam_nb, 
                     const arma::vec& lgy1,const arma::uword& max_iter = 4e3,
                     const double& eps = 1e-7, const bool& show = false){
    arma::uword iter;
    arma::uword converge = 0;
    arma::uword pp = X.n_cols+fam_nb;
    arma::vec offsets = arma::zeros<arma::vec>(z.n_elem);
    Rcpp::List new_reg, new_bxj_fit;
    double new_LL, new_bxj, phi, curr_LL, LL0;
    double curr_bxj = 0.0;
    arma::vec curr_reg_par = arma::zeros<arma::vec>(pp);
    arma::vec new_reg_par = arma::zeros<arma::vec>(pp);
    arma::vec BETA = arma::zeros<arma::vec>(pp-1);
    
    // initial regression fit 
    new_reg = Rcpp_reg_BFGS(y, X, offsets, curr_reg_par, fam_nb, lgy1, 
                                max_iter, eps, show);
    curr_reg_par = as<arma::vec>(new_reg["PAR"]);
    curr_LL = as<double>(new_reg["LL"]);
    LL0     = as<double>(new_reg["LL"]);
    BETA    = curr_reg_par.subvec(0, X.n_cols-1);
    phi     = std::exp(curr_reg_par.at(pp-1));
    
  while(iter < max_iter){
    
    //update bxj
    new_bxj_fit = Rcpp_trec_bxj(y, X, curr_bxj, z, BETA, phi, fam_nb, lgy1,
                                max_iter, eps, show);
    new_bxj     = as<double>(new_bxj_fit["bxj"]);
    new_LL      = as<double>(new_bxj_fit["LL"]);

    if(new_LL < curr_LL) printR_obj("likelihood decreased for bxj");
    //printR_obj(curr_LL);
    
    //update BETA, phi 
    compute_offset(new_bxj, z, offsets);
    new_reg     = Rcpp_reg_BFGS(y, X, offsets, curr_reg_par, fam_nb, lgy1, 
                                max_iter, eps, show);
    
    new_LL      = as<double>(new_reg["LL"]);
    new_reg_par = as<arma::vec>(new_reg["PAR"]);
    
    BETA = new_reg_par.subvec(0, X.n_cols-1);
    phi  = std::exp(new_reg_par.at(pp-1));
    
    if(new_LL < curr_LL) printR_obj("likelihood decreased for betas");
    //printR_obj(curr_LL);
    
    if(iter > 0){
      if( std::abs(curr_LL - new_LL) < eps && 
          Rcpp_norm(curr_reg_par - new_reg_par) < eps &&
          std::abs(curr_bxj - new_bxj) < eps){
          converge = 1;
          break;
      }
    }
    
    curr_reg_par = new_reg_par;
    curr_bxj     = new_bxj;
    curr_LL      = new_LL;
    iter++;
    
  }
  return Rcpp::List::create(
      Rcpp::Named("bxj", new_bxj),
      Rcpp::Named("reg_par", Rcpp::NumericVector(new_reg_par.begin(), new_reg_par.end())),
      Rcpp::Named("LL", new_LL),
      Rcpp::Named("lrt", (new_LL-LL0)*2),
      Rcpp::Named("converge", converge)
    );
}

/* ---------------------------
 * ASE
 ---------------------------*/

// [[Rcpp::export]]
double Rcpp_loglikBB(const arma::vec& ni, const arma::vec& ni0, 
                const double& Pi, const double& theta, 
                const arma::vec& lbc){
  
  arma::uword ii;  
  double loglik = 0.0;
  double vtheta = 1/theta;
  double aa = Pi*vtheta;
  double bb = vtheta - aa;
  double lgvab = - lgamma(aa) - lgamma(bb) + lgamma(vtheta);
  
  if(Pi == 0.5){
    for(ii=0;ii<ni.n_elem;ii++){
      printR_obj(loglik);
      loglik += lbc.at(ii) + lgamma(aa + ni0.at(ii)) +
        lgamma(bb + ni.at(ii) - ni0.at(ii)) + lgvab - 
        lgamma(vtheta + ni.at(ii));
    } 
  }else{
    for(ii=0;ii<ni.n_elem;ii++){
    printR_obj(loglik);
    loglik += lbc.at(ii) + lgamma(aa + ni0.at(ii)) +
      lgamma(bb + ni.at(ii) - ni0.at(ii)) - lgamma(vtheta + ni.at(ii)) + lgvab;
  }
    }
  
  return loglik;
}

// [[Rcpp::export]]
double Rcpp_log_BB(const arma::uword& x,
                   const arma::uword& n,const double& pp,
                   const double& psi,const double& lbc){
  
  double aa = pp / psi;
  double bb = (1.0 - pp) / psi;
  
  return lbc + lgamma(x + aa) +
    lgamma(n - x + bb) +
    lgamma(aa + bb) -
    lgamma(n + aa + bb) -
    lgamma(aa) - lgamma(bb);
}

// [[Rcpp::export]]
double Rcpp_vec_log_BB(const arma::uvec& x,
                       const arma::uvec& n,const double& pp,
                       const double& psi,const arma::vec& lbc){
  
  double out = 0.0;
  arma::uword ii;
  for(ii = 0; ii < x.n_elem; ii++){
    out += Rcpp_log_BB(x.at(ii),n.at(ii),
                       pp,psi,lbc.at(ii));
  }
  return out;
}

// [[Rcpp::export]]
double Rcpp_ase_bxj(const arma::vec& ni, const arma::vec& ni0,
                    const double& Pi, const double& theta){

  arma::uword ii;
  double grad = 0.0;
  double vtheta = 1/theta;
  double aa = Pi*vtheta;
  double bb = vtheta - aa;

  for(ii=0;ii<ni.n_elem;ii++){
    grad += R::digamma(aa + ni0.at(ii)) -
      R::digamma(bb + ni.at(ii) - ni0.at(ii)) - R::digamma(aa) +
      R::digamma(bb);
  }


}
