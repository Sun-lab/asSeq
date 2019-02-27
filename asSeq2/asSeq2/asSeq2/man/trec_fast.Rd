\name{trec_fast}
\alias{trec_fast}
\title{
  Total Read Count (TReC) model
}
\description{
  eQTL mapping of total read count (from both alleles of one gene) 
  by generalized linear model.
}
\usage{
trec(Y, X, Z, output.tag, p.cut, link="log", offset = NULL, 
  adjZ=TRUE, local.only = TRUE, local.distance = 2e+05, eChr = NULL, ePos = NULL, 
  mChr = NULL, mPos = NULL, converge=5e-5, convergeGLM=1e-8, scoreTestP=0.05,
  trace = 1, maxit = 100)
}
\arguments{
  \item{Y}{ matrix of gene expression (total read count per gene per sample). 
    Each row corresponds to one sample and each column corresponds to one gene. 
  }
  \item{X}{ matrix of cofounding covariates. 
    Each row corresponds to one sample and each column corresponds to one gene.
    Usually this include log(total reads per sample) 
    and PCs obtained from log transformed expression matrix.
  }
  \item{Z}{ matrix of genotype data, Each row corresponds to one sample 
    and each column corresponds to one SNP. 
  }
  \item{output.tag}{
    tag of output file. The actual output files will be tag_eqtl.txt and 
    tag_freq.txt.
  }
  \item{p.cut}{ p-value cut-off, only the associations with p-value smaller 
  	than p.cut are saved. 
  }
  \item{link}{ link function, only 'log' and 'identity' are allowed. 
  A 'identity' link is not recommended since it may lead to errors 
  in model fitting because mu can be negative. 
  }

  \item{offset}{
    A vector of offset of each sample in the model linkfun(mu) = offset + Xb.
  }
  
  \item{adjZ}{
   whether to adjust Z, if there is no adjustment, we assume a linear relation
   linkfun(mu) ~ Z; otherwise a linear relation mu ~ Z;
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
  
  \item{converge}{ convergence tolerance. The algorithm converges if changes 
  of all the parameters are smaller than value. }
   
  \item{convergeGLM}{ convergence tolerance to be passed to function glmNB.}

  \item{scoreTestP}{ p-value cutoff for score test of the over dispersion in 
  Poisson model. For each association regression, we first fit a Poisson model, 
  and move to negative binomial if the score test p-value smaller than 
  scoreTestP}

  \item{trace}{ If trace=0, do not output tracing information.  
  if trace = 1, output brief tracing information, and 
  if trace > 1. output detailed tracing information. }

  \item{maxit}{ maximum number of iterations used in glmNB  }
}
\details{
  For each gene, first a base-line model only including the covariates in X 
  is fitted using negative binomial regression. If the over dispersion parameter
  is not significantly different from zero, a Poisson model is used for all
  the eQTL mapping of this gene.  
}
\value{
The function return 1 if no error occurs.
 The results are written into two output files: 
 one file (output.tag_eqtl.txt) stores the eQTL mapping result that
 pass the p-value threshold p.cut, including gene ID,
 marker ID, Family (Poisson and Negative Binomial), 
 Chi-square statistic, and p-value; 
 the other file (output.tag_freq.txt) includes a frequency table 
 of all the p-values.}
\author{ Wei Sun weisun@email.unc.edu }
\keyword{ methods }

