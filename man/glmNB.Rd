\name{glmNB}
\alias{glmNB}
\title{
Fit a generalized linear model for negative binomial family.
}
\description{
Fit a generalized linear model for negative binomial family, and return 
the fitted values as well as Pearson residuals.
}
\usage{
glmNB(y, X, link, offset = NULL, naPercent = 0.4, maxit = 50, 
  conv = 1e-08, fitted = NULL, scoreTestP=0.05,  trace = 1)
}
\arguments{
  \item{y}{a vector of the response variable. }
  
  \item{X}{ a matrix of covariate data, each column is one covariate and 
    each row is one observation. Do not include intercept in 
    this covariate matrix}
  
  \item{link}{
	a character string specifying the link function, the allowed values 
	include "logit", "log", "identity", and "inverse".}
  
  \item{offset}{a vector specifying offset in the linear predictor scale, 
    eta = Xb + offset.}
  
  \item{naPercent}{the maximum percent of missing that are allowed in the 
  data. }
  
  \item{maxit}{the maximum number of iteration}
  
  \item{conv}{convergence tolerance. The algorithm converges if the change of  
    likelihood and over dispersion parameter are both smaller than conv}
  
  \item{fitted}{a vector of the fitted value, should be the same length as y.}
  
  \item{scoreTestP}{ p-value cutoff for score test of the over dispersion in 
    Poisson model. For each association regression, we first fit a Poisson 
    model, and move to negative binomial if the score test p-value smaller 
    than scoreTestP}

  \item{trace}{ Print out tracking information if trace > 1.}
}
\value{
A list of the following components:
  \item{phi }{over dispersion parameter, phi equals to 0 when there is no over dispersion}
  
  \item{prior }{prior weight}
  
  \item{offset }{offset}
  
  \item{loglik }{2 times log likelihood}
  
  \item{rank }{rank of the covariate matrix}
  
  \item{fitted }{a vector of fitted values, which equals NA if response or 
  any covariate is NA}
  
  \item{resid }{a vector of residual values}
  
  \item{weights }{a vector of weights}
  
  \item{scale }{scale parameter of exponential family}
  
  \item{dfResid }{degree of freedom of the residuals}
  
  \item{nIter }{number of iterations of the IRLS algorithm}
  
  \item{X }{covariate matrix}
}
\author{
Wei Sun weisun@email.unc.edu
}
\note{
This is implemented by a C version of the R code in MASS/glm.nb
}

\keyword{ methods }


