#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <omp.h>

// [[Rcpp::depends("RcppArmadillo")]]

// Get genotype quicker

// [[Rcpp::export]]
arma::uvec Rcpp_dose2geno(const arma::mat& doses,
	const double& geno_thres = 0.8){
	
	arma::uword num_loci = doses.n_rows, ii, tmp_index;
	arma::uvec geno = arma::zeros<arma::uvec>(num_loci);
	arma::vec tmp_dose = arma::zeros<arma::vec>(3);
	
	for(ii = 0; ii < num_loci; ii++){
		tmp_dose = doses.row(ii).t();
		tmp_index = arma::index_max(tmp_dose);
		if( tmp_dose.at(tmp_index) >= geno_thres ){
			geno.at(ii) = tmp_index; // 0,1,2
		} else {
			geno.at(ii) = 5; // aka NA, missing genotype
		}
	}
	
	return geno;
}