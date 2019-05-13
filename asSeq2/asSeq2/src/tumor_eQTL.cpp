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
void compute_offset(const arma::vec& z, const arma::vec& RHO, const int& n,
                    double& KAPPA, double& ETA, double& GAMMA,
                    const arma::vec& tau1, const arma::vec& tau2,
                    arma::vec& offsets){
  arma::uword ii;
  for(ii = 0; ii<n; ii++){
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
void compute_mu(const arma::mat& X, const arma::vec& offsets,
                const arma::vec& betas, const int& n,
                arma::vec& mu){
  arma::uword ii;
  for(ii = 0; ii<n; ii++){
    mu.at(ii) = std::exp(arma::dot(X.row(ii).t(), betas + offsets.at(ii)));
  }
}

// [[Rcpp::export]]
double Rcpp_loglikNB(const double& phi, const int& n, const arma::vec& lgy1,
                const arma::vec& mu, const arma::vec& y){
  //lgy1 = std::lgamma(y + 1)
  double loglik = 0.0;
  arma::uword ii;
  double vphi = 1.0/phi;
  double lgvphi = lgamma(vphi);

  for(ii = 0; ii<n; ii++){
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
arma::vec Rcpp_grad_NB(double& KAPPA, double& ETA, double& GAMMA, const int& H0,
                       const int& n, const arma::vec& y, const arma::vec& z,
                       const arma::vec& RHO, const arma::mat& X,
                       const arma::vec& betas, double& phi,
                       const arma::vec& tau1, const arma::vec& tau2,
                       const arma::vec& offsets, const arma::vec& mu){

  double dltrec_dmu = 0, dmu_dkappa = 0, dmu_deta = 0, dmu_dgamma = 0,
    expXbeta;
  arma::uword ii;
  arma::vec grad = arma::zeros<arma::vec>(3);

  for(ii =0; ii<n; ii++){
    expXbeta   = std::exp(arma::dot(X.row(ii).t(), betas));
    dltrec_dmu = y.at(ii)/mu.at(ii) - (1 + phi * y.at(ii))/(1 + phi * mu.at(ii));

      if(z.at(ii) == 0){

        dmu_dkappa = (tau1.at(ii) + tau2.at(ii))*expXbeta*RHO.at(ii);

      }else if(z.at(ii) == 1){

        dmu_dkappa = (tau1.at(ii) + tau2.at(ii)*GAMMA)*expXbeta*RHO.at(ii);
        dmu_deta   = (1-RHO.at(ii))*expXbeta;
        dmu_dgamma = tau2.at(ii)*RHO.at(ii)*KAPPA*expXbeta;

      }else if(z.at(ii) == 2){

        dmu_dkappa = (tau1.at(ii)*GAMMA + tau2.at(ii))*expXbeta*RHO.at(ii);
        dmu_deta   = (1-RHO.at(ii))*expXbeta;
        dmu_dgamma = tau1.at(ii)*RHO.at(ii)*KAPPA*expXbeta;

      }else{
        dmu_dkappa = (tau1.at(ii)+tau2.at(ii))*GAMMA*expXbeta*RHO.at(ii);
        dmu_deta   = 2*(1-RHO.at(ii))*expXbeta;
        dmu_dgamma = (tau1.at(ii)+tau2.at(ii))*KAPPA*expXbeta*RHO.at(ii);
      }

    grad.at(0) += dltrec_dmu*dmu_dkappa;
    grad.at(1) += dltrec_dmu*dmu_deta;
    grad.at(2) += dltrec_dmu*dmu_dgamma;

  }
  grad.at(0) *= KAPPA;
  grad.at(1) *= ETA;
  grad.at(2) *= GAMMA;

  return grad;
}
