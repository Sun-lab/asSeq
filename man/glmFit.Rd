\name{glmFit}
\alias{glmFit}
\title{
Fit a generalized linear model.
}
\description{
Fit a generalized linear model, and return the fitted values 
as well as Pearson residuals. Now only support faimliy 
"binomial", "poisson", "gaussian", and "gamma".
}
\usage{
glmFit(family, link, y, X, offset = NULL, naPercent = 0.4, nTotal = NULL, 
       maxit = 50, conv = 1e-05, fitted = NULL, trace = 1)
}
\arguments{
  \item{family}{
	a character string specifying the distribution family, the allowed values 
	include "binomial", "poisson", "gaussian", and "gamma".}
  
  \item{link}{
	a character string specifying the link function, the allowed values 
	include "logit", "log", "identity", and "inverse".}
  
  \item{y}{
	a vector of the response variable. }
  
  \item{X}{
	a matrix of covariate data, each column is one covariate and 
	each row is one observation. Do not include intercept in 
	this covariate matrix}
  
  \item{offset}{
	a vector specifying offset in the linear predictor scale, 
	eta = Xb + offset.}
  
  \item{naPercent}{
	the maximum percent of missing that are allowed in the data. }
  
  \item{nTotal}{
	if faimly is "binomial", nTotal specify the total number of trails 
	in each subject}
  
  \item{maxit}{
	the maximum number of iteration}
  
  \item{conv}{
	convergence tolerance. Let wss be the weighted sum of perason 
	residual square, the algorithm converges if 
	|wss-wss_last|/wss_last < conv}
  
  \item{fitted}{
	a vector of the fitted value, should be the same length as y}
  
  \item{trace}{ Print out tracking information if trace > 1.}
}

\value{
A list of the following components:
  \item{family }{distribution family}
  \item{link }{link function}
  \item{offset }{offset}
  \item{rank }{rank of the covariate matrix}
  \item{fitted }{a vector of fitted values, which equals NA if response or 
  any covariate is NA}
  \item{resid }{a vector of residual values}
  \item{weights }{a vector of weights}
  \item{scale }{scale parameter of exponential family}
  \item{dfResid }{degree of freedom of the residuals}
  \item{nIter }{number of iterations of the IRLS algorithm}
  \item{X }{covariate matrix}
  \item{beta }{the coefficient of the last covariate in the model}
}
\author{
Wei Sun weisun@email.unc.edu
}
\note{
Part of codes is adopted from CNVtools/glm_test.cpp
}

\keyword{ methods }

