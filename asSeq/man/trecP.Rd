\name{trecP}
\alias{trecP}
\title{
  Permutation test for TReC model
}

\description{
  eQTL mapping of total read count by generalized linear model. The total 
  read count is in contrast to the allele specific read count. 
}

\usage{
trecP(Y, X, Z, link="log", weight = NULL, offset = NULL, adjZ=TRUE, 
  local.only = TRUE, local.distance = 2e+05, eChr = NULL, ePos = NULL, 
  mChr = NULL, mPos = NULL, converge = 1e-08, scoreTestP = 0.05, 
  np.max=5000, np=c(20, 100, 500, 1000, 2500),
  aim.p=c(0.5, 0.2, 0.1, 0.05, 0.02), 
  confidence.p = 0.01, trace = 0, maxit = 100)
}

\arguments{
\item{Y}{ matrix of gene expression (total read count per gene per sample). 
    Each row corresponds to one sample and each column corresponds to one gene. 
  }
    
  \item{X}{ matrix of cofounding covariates. 
    Each row corresponds to one sample and each column corresponds to one 
    gene. Usually this include log(total reads per sample) 
    and PCs obtained from log transformed expression matrix.}
    
  \item{Z}{ matrix of genotype data, Each row corresponds to one sample 
    and each column corresponds to one SNP. }
  
  \item{link}{ link function, only 'log' and 'identity' are allowed. 
    A 'identity' link is not recommended since it may lead to errors 
    in model fitting because mu can be negative. 
  }

  \item{weight}{ A vector of weight of each sample.}
    
  \item{offset}{ A vector of offset of each sample in the model 
    log(mu) = offset + Xb.}
  
  \item{adjZ}{
   whether to adjust Z, if there is no adjustment, we assume a linear 
   relation linkfun(mu) ~ Z; otherwise a linear relation mu ~ Z;
  }

  \item{local.only}{ whether examine local eQTL only. }
  
  \item{local.distance}{ if the distance between a gene and a marker
    is smaller than local.distance, we call it a local eQTL. }
   
  \item{eChr}{ chromosomes of gene expression traits, must be
    a vector of integers, which means sex chromosomes (X, Y) need
    to be coded as numbers }
   
  \item{ePos}{ positions of gene expression traits, must be
    a vector of integers}
   
  \item{mChr}{ chromosomes of markers, must be a vector of integers,
    which means sex chromosomes (X, Y) need to be coded as numbers }
   
  \item{mPos}{ positions of markers, must be a vector of integers}
  
  \item{converge}{ convergence tolerance to be passed to function glmNB. 
    The iterative algorithm in glmNB converges if the change of likelihood  
    and over dispersion parameter are both smaller than converge }
   
  \item{scoreTestP}{ p-value cutoff for score test of the over dispersion in 
    Poisson model. For each association regression, we first fit a 
    Poisson model, and move to negative binomial if the score test p-value 
    smaller than scoreTestP}

  \item{np.max}{ the maximum number of permutations }
   
  \item{np}{ number of permutations, a vector in ascending order, 
    see the explanation for \code{aim.p} }
   
  \item{aim.p}{ the permutation p-value threshold, corresponding
    to the number of permutations in \code{np}. We stop the 
    permutations if we can confident to say the permutation p-value 
    is larger than aim.p[i] after np[i] permutations. }
   
  \item{confidence.p}{ A binomial p-value determining whether
    the permutation p-value estimated from the np[i] permutations
    is larger than aim.p[i]. }

  \item{trace}{ If trace=0, do not output tracing information.  
    if trace = 1, output brief tracing information, and 
    if trace > 1. output detailed tracing information. }

  \item{maxit}{ maximum number of iterations used in glmNB. }
  
}
\details{
  Different amount of permutations are performed for the each gene  
}
\value{
  A data frame with columns gene ID, marker ID,
  nominal p-value, permutation p-value and number of permutations
}
\author{ Wei Sun wsun@bios.unc.edu }
\keyword{ methods }

