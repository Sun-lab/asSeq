\name{get_reduced_boot}
\alias{get_reduced_boot}
\title{
  A bootstrap of the gene expression with reduced effect size.
}
\description{
  Producing a bootstrap of the gene expression with reduced effect size. 
  The procedure is setup to take MatrixEQTL style results extended with the 
  number of effective tests and produce a set of the reduced bootstraps 
  with target minimum p-values.
}
\usage{
get_reduced_boot(ppi, target.min.ps=NULL, target.perm.ps, genenm, mQTL.fit, 
                   expr.mat, min.SNP, covars, nsam, use.norm=F)
}
\arguments{
  \item{ppi}{Which bootstrap from the grid to process.}
  \item{target.min.ps}{A grid of minimum p-values (by default is null, and is 
     populated using target permutation p-values and effective number of SNPs}
  \item{target.perm.ps}{A grid of target permutation p-values (requires a good
     guess of effective number of SNPs}
  \item{genenm}{gene name}
  \item{mQTL.fit}{MatrixEQTL results extended with effective number of tests guess}
  \item{expr.mat}{MatrixEQTL style gene expression of the target gene}
  \item{min.SNP}{MatrixEQTL style minimum SNP information for the target gene}
  \item{covars}{MatrixEQTL style covariates in TReC part of the model needed to produce a parametric bootstrap}
  \item{nsam}{number of samples used in the analysis}
  \item{use.norm}{instead of residuals from the fit use errors simulated 
  from normal distribution (by default FALSE)}
}
  
\details{
  The function produces a bootstrap of the gene expression with reduced effect size. 
  The procedure is setup to take MatrixEQTL style results extended with the 
  number of effective tests and produce a set of the reduced bootstraps 
  with target minimum p-values. 
  Fitting linear model versus minimum SNP we obtain effect size estimate X, 
  which is adjusted according to the target minimum p-value.
}
\value{
 The function returns a vector with a reduced effect size bootstrap
}

%\examples{
%redboot = get_reduced_boot(1, target.perm.ps=1e-2, i=1,
%mQTL.fit=eigenMT,
%expr.mat = permEst$geneM,
%min.SNP=gen.sub,
%covars=permEst$cvrtM,
%nsam=ncol(gen.sub))
%}

\author{ Vasyl Zhabotynsky, Wei Sun weisun@email.unc.edu }
\keyword{ methods }

