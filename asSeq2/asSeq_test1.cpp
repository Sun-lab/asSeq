#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>

// [[Rcpp::depends("RcppArmadillo")]]
/* ---------------------------
 * 
 ---------------------------*/
/* ---------------------------
 * TReC
 ---------------------------*/

using namespace Rcpp;

// [[Rcpp::export]] 
double loglikNB(const double& phi, const int& n, const arma::vec& lgy1,
                const arma::vec& mu, const arma::vec& y){
  //lgy1 = std::lgamma(y + 1)
  double loglik = 0.0;
  arma::uword ii;
  double vphi = 1.0/phi;
  
  for(ii = 0; ii<n; ii++){
    if(y.at(ii) > 0){
      loglik += lgamma(y.at(ii) + vphi) -lgamma(vphi) - lgy1.at(ii) +
        y.at(ii) * std::log(mu.at(ii));
    }
    loglik += vphi*std::log(vphi) - 
      (vphi+y.at(ii)) * std::log(vphi+mu.at(ii));
  }
  
  return loglik;
}

// [[Rcpp::export]] 
double loglik_pois(const int& n, const arma::vec& lgy1,
                   const arma::vec& mu, const arma::vec& y){
  //lgy1 = std::lgamma(y + 1)
  double loglik = 0.0;
  arma::uword ii;
  for(ii =0; ii<n; ii++){
    loglik += y.at(ii)*std::log(mu.at(ii)) - mu.at(ii) - lgamma(y.at(ii)+1);
  }
  return loglik;
}

// [[Rcpp::export]] 
double logLTReC(const double& bxj, const int& n, const arma::vec& lgy1,
                const arma::vec& y, const arma::vec& z,
                const arma::vec& mu, const double b0, const double& phi,
                const bool& fam_nb, arma::vec& mu1, arma::vec& offsets){
  // z is the genotype vector take value 0,1,2 (same as x in the R code)
  // lgy1 = lgamma(y+1)
  arma::uword ii; 
  //arma::vec mu1 = arma::zeros<arma::vec>(n);
  
  for(ii =0; ii<n; ii++){
    if(z.at(ii) == 2){
      offsets.at(ii) = bxj;
      mu1.at(ii) = mu.at(ii)*std::exp(bxj - b0);
    }else if(z.at(ii) == 1){
      offsets.at(ii) = (1+std::exp(bxj))/2;
      mu1.at(ii) = mu.at(ii)*(1+std::exp(bxj))/(1+std::exp(b0));
    }else{
      mu1.at(ii) = mu.at(ii);
    }
  }
  
  if(fam_nb){
    return loglikNB(phi,n,lgy1, mu1, y);
  }else{
    return loglik_pois(n,lgy1, mu1, y);
  }
}

// [[Rcpp::export]]
double grad_bxj_trec(const double& bxj, const int& n, const arma::vec& y,
                     const arma::vec& z, const arma::vec& mu1, const double b0,
                     const double& phi, const bool& fam_nb ){
  double grad = 0.0, dg_dmu =0.0, dmu_db = 0.0;
  arma::uword ii;
  
  for(ii =0; ii<n; ii++){
    if(fam_nb){
      dg_dmu = y.at(ii)/mu1.at(ii) - (1.0+phi*y.at(ii))/(1.0+phi*mu1.at(ii));
    }else{
      dg_dmu = y.at(ii)/mu1.at(ii) - 1.0;
    }
    if(z.at(ii) == 2){
      dmu_db = mu1.at(ii);
    }else if(z.at(ii) == 1){
      dmu_db = mu1.at(ii)*std::exp(bxj)/(1.0+std::exp(bxj));
    }
    
    grad += dg_dmu*dmu_db;
  }
  return grad;
}
