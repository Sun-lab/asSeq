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
double Rcpp_logLTReC(const double& bxj, const arma::vec& lgy1,
                const arma::vec& y, const arma::vec& z,
                const arma::vec& mu, const double& phi,
                const bool& fam_nb, arma::vec& mu1, arma::vec& offsets){
  // z is the genotype vector take value 0,1,2 (same as x in the R code)
  // lgy1 = lgamma(y+1)
  arma::uword ii; 
  //arma::vec mu1 = arma::zeros<arma::vec>(y.n_elem);
  
  for(ii =0; ii<y.n_elem; ii++){
    if(z.at(ii) == 2){
      offsets.at(ii) = bxj;
      //mu1.at(ii) = mu.at(ii)*std::exp(bxj - b0);
      mu1.at(ii) = mu.at(ii)*std::exp(bxj);
    }else if(z.at(ii) == 1){
      offsets.at(ii) = std::log((1+std::exp(bxj))/2);
      //mu1.at(ii) = mu.at(ii)*(1+std::exp(bxj))/(1+std::exp(b0));
      mu1.at(ii) = mu.at(ii)*(1+std::exp(bxj))/2;
    }else{
      mu1.at(ii) = mu.at(ii);
    }
  }
  
  if(fam_nb){
    return Rcpp_loglikNB(phi, mu1, y, lgy1);
  }else{
    return Rcpp_loglik_pois(mu1, y, lgy1);
  }
}

// [[Rcpp::export]]
double Rcpp_grad_bxj_trec(const double& bxj, const arma::vec& y,
                     const arma::vec& z, const arma::vec& mu1,
                     const double& phi, const bool& fam_nb ){
  double grad = 0.0, dg_dmu =0.0, dmu_db = 0.0;
  arma::uword ii;
  
  for(ii =0; ii<y.n_elem; ii++){
    if(fam_nb){
      dg_dmu = y.at(ii)/mu1.at(ii) - (1.0+phi*y.at(ii))/(1.0+phi*mu1.at(ii));
    }else{
      dg_dmu = y.at(ii)/mu1.at(ii) - 1.0;
    }
    if(z.at(ii) == 2){
      dmu_db = mu1.at(ii);
    }else if(z.at(ii) == 1){
      dmu_db = mu1.at(ii)*std::exp(bxj)/(1.0+std::exp(bxj));
    }else{
      dmu_db = 0.0;
    }
    
    grad += dg_dmu*dmu_db;
  }
  return grad;
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

// [[Rcpp::export]]
Rcpp::List Rcpp_NB_reg(const arma::vec& y,
                       const arma::mat& X,const arma::vec& offsets,
                       const arma::vec& params0, const arma::vec& lgy1,
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
    old_LL = Rcpp_NB_reg_LL(y,X,offsets,old_PARAMS, lgy1, mu);
    old_grad = Rcpp_NB_reg_grad(y,X,mu,old_PARAMS);
    old_hess_grad = -1.0 * arma::inv(Rcpp_NB_reg_Hess(y,X,mu,old_PARAMS)) * old_grad;
    old_grad = old_grad / std::max(1.0,Rcpp_norm(old_grad));
    old_hess_grad = old_hess_grad / std::max(1.0,Rcpp_norm(old_hess_grad));
    uu = 0;
    for(jj = 0; jj <= 15; jj++){
      new_PARAMS = old_PARAMS + old_hess_grad / std::pow(4.0,jj);
      new_LL = Rcpp_NB_reg_LL(y,X,offsets,new_PARAMS, lgy1,mu);
      if( new_LL > old_LL ){
        old_PARAMS = new_PARAMS;
        old_LL = new_LL;
        uu = 1;
        break;
      } else {
        new_PARAMS = old_PARAMS + old_grad / std::pow(4.0,jj);
        new_LL = Rcpp_NB_reg_LL(y,X,offsets,new_PARAMS, lgy1,mu);
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
      if( std::abs(curr_LL - old_LL) < eps && Rcpp_norm(curr_PARAMS - old_PARAMS) < eps ){
        old_grad = Rcpp_NB_reg_grad(y,X,mu,old_PARAMS);
        old_hess_grad = -1.0 * arma::inv(Rcpp_NB_reg_Hess(y,X,mu,old_PARAMS)) * old_grad;
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
Rcpp::List Rcpp_NB_reg_BFGS(const arma::vec& y,
                            const arma::mat& X,const arma::vec& offsets,
                            const arma::vec& params0, const arma::vec& lgy1,
                            const arma::uword& max_iter = 4e3,
                            const double& eps = 1e-7,const bool& show = true){
  
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
    old_LL = fnscale * Rcpp_NB_reg_LL(y, X, offsets, xk, lgy1, mu);
    gr_k = fnscale * Rcpp_NB_reg_grad(y, X, mu, xk);
    p_k = -1.0 * inv_Bk * gr_k;
    inv_norm_p_k =  1.0 / std::max(1.0, Rcpp_norm(p_k));
    
    //line search for new xk
    for(jj=0; jj<30; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4, jj);
      new_xk = xk + tmp_alpha * p_k;
      new_LL = fnscale * Rcpp_NB_reg_LL(y, X, offsets, new_xk, lgy1, mu);
      if(new_LL < old_LL){ //minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * Rcpp_NB_reg_grad(y, X, mu, new_xk);
        ys = arma::dot(y_k, s_k);
        if(ys > 0.0){
          if(show) printR_obj("Update xk and inv_Bk");
          ISYT = I_num_params - (s_k * y_k.t()) /ys;
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
        gr_k = Rcpp_NB_reg_grad(y, X, mu, xk);
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
  
  old_LL = Rcpp_NB_reg_LL(y, X, offsets, xk, lgy1, mu);
  //gr_k = Rcpp_NB_reg_grad(y, X, mu, xk);
  return Rcpp::List::create(
    Rcpp::Named("converge", converge),
    Rcpp::Named("LL", old_LL),
    Rcpp::Named("iter", iter),
    Rcpp::Named("norm_GRAD", Rcpp_norm(gr_k)),
    Rcpp::Named("PAR", Rcpp::NumericVector(xk.begin(), xk.end()))
  );
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
  
// [[Rcpp::export]]
Rcpp::List Rcpp_pois_reg(const arma::vec& y,
                         const arma::mat& X,const arma::vec& offsets,
                         const arma::vec& params0, const arma::vec& lgy1,
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
    old_LL = Rcpp_pois_reg_LL(y,X,offsets,old_PARAMS, lgy1, mu);
    old_grad = Rcpp_pois_reg_grad(y,X,mu,old_PARAMS);
    old_hess_grad = -1.0 * arma::inv(Rcpp_pois_reg_Hess(y,X,mu,old_PARAMS)) * old_grad;
    old_grad = old_grad / std::max(1.0,Rcpp_norm(old_grad));
    old_hess_grad = old_hess_grad / std::max(1.0,Rcpp_norm(old_hess_grad));
    uu = 0;
    for(jj = 0; jj <= 15; jj++){
      new_PARAMS = old_PARAMS + old_hess_grad / std::pow(4.0,jj);
      new_LL = Rcpp_pois_reg_LL(y,X,offsets,new_PARAMS, lgy1,mu);
      if( new_LL > old_LL ){
        old_PARAMS = new_PARAMS;
        old_LL = new_LL;
        uu = 1;
        break;
      } else {
        new_PARAMS = old_PARAMS + old_grad / std::pow(4.0,jj);
        new_LL = Rcpp_pois_reg_LL(y,X,offsets,new_PARAMS, lgy1,mu);
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
      if( std::abs(curr_LL - old_LL) < eps && Rcpp_norm(curr_PARAMS - old_PARAMS) < eps ){
        old_grad = Rcpp_pois_reg_grad(y,X,mu,old_PARAMS);
        old_hess_grad = -1.0 * arma::inv(Rcpp_pois_reg_Hess(y,X,mu,old_PARAMS)) * old_grad;
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
/*
// [[Rcpp::export]]
Rcpp::List Rcpp_trec_BFGS(const arma::vec& y, const arma::mat& X, 
                          const arma::vec& offsets, const arma::vec& params0, 
                          const arma::vec& lgy1, const arma::uword& max_iter = 4e3, 
                          const double& eps = 1e-7,const bool& show = true){
  
}

*/

/* ---------------------------
 * ASE
 ---------------------------*/
