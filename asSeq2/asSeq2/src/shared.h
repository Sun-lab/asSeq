#ifndef _SHARED_H_
#define _SHARED_H_

#include <RcppArmadillo.h>
#include <Rmath.h>

// [[Rcpp::depends("RcppArmadillo")]]

// --------------------
// Intermediate Functions
// --------------------

double Rcpp_norm(const arma::vec& a);

double Rcpp_logSumExp(const arma::vec& log_x);

double Rcpp_min_diff(const arma::vec& x);

arma::vec Rcpp_lgy_add_1(const arma::vec& y);


/* ---------------------------
* NB regression (duplicated function)
*/
double Rcpp_reg_LL(const arma::vec& y, const arma::mat& X,
                   const arma::vec& offsets, const arma::vec& PARAMS,
                   const arma::vec& lgy1, arma::vec& mu);

arma::vec Rcpp_reg_grad(const arma::vec& y, const arma::mat& X,
                        const arma::vec& mu, const arma::vec& PARAMS);

Rcpp::List Rcpp_reg_BFGS(const arma::vec& y, const arma::mat& X,
                         const arma::vec& offsets, const arma::vec& params0,
                         const arma::vec& lgy1,
                         const arma::uword& max_iter = 4e3,
                         const double& eps = 1e-7, const bool& show = true);

double Rcpp_loglikNB(const arma::vec& y, const double& phi, 
                     const arma::vec& lgy1, const arma::vec& mu);




#endif