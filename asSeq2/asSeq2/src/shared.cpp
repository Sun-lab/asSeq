#include "shared.h"


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


// [[Rcpp::export]]
double Rcpp_loglikNB(const arma::vec& y, const double& phi, 
                      const arma::vec& lgy1, const arma::vec& mu){ //duplicated 
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

