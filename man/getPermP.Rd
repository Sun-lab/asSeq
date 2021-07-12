\name{getPermP}
\alias{getPermP}
\title{
  Estimating permutation p-value for a given gene. 
}
\description{
  Estimating permutation p-value for a given gene. The procedure is setup to 
  take MatrixEQTL compatible results for one gene at a time.
}
\usage{
getPermP(permEst, num.grid = 100, slice=2000, n.perm=1000, ini.perm=100, 
                  nsub=NA, write.res=F)
}
\arguments{
  \item{permEst}{A list containing gene-relevant inforomation in MatrixEQTL compatible format.
  \itemize{
    \item snpM - matrix with genotyping data
    \item geneM - gene expression data
    \item cvrtM - covariates
    \item snpspos - SNP names and positions
    \item genepos - gene name and position
    \item outpf - output file name for the results,
    \item pvOutputThreshold pvalue cutoff for trans tests (default 1e-300),
    \item pvOutputThreshold.csv - pvalue cutoff for cis tests (defualt 1),
    \item cisDist - distance to inform MatrixEQTL which SNPs are cis (default 1e9 since it is more efficient to preprocess SNPs to submit only cis SNPs to the procedure)
    \item effNumGuess - initial guess of effective number of SNPs 
    (can be obtained from some fast external tool such as eigenMT 
    or we found using 1/4 of the SNPs as a guess). Procedure starts with the 
    initial guess and adjusts it depending on the observed results in first 100 permutations.
    \item verbose (default FALSE) - parameter informing MatrixEQTL to be less verbose,
    \item pvalue.hist (default FALSE) - parameter informing MatrixEQTL to just run the main analysis
    \item min.pv.by.genesnp (default FALSE)
    \item noFDRsaveMemory (default FALSE)
    \item outdir - supply output directory
  } 
  }
  
  \item{num.grid}{target number of grid points for permutation p-value estimation 
  (default 100) 
  }
  
  \item{slice}{Technical MatrixEQTL parameter for number of genes/snps per slice 
  (default 2000 as not too mermory intensive)
  }

  \item{ini.perm}{First N permutations used to evaluate appropriateness of the 
  initial number of effective number of SNPs and populating the grid.
  }
  
  \item{nsub}{If you want to run permutations on a subset of samples 
  (default NA letting the program to run estimation on all samples)
  }
  
  \item{write.res}{
    informs whether the results should be written to the hard drive as opposed to just returning an object
  }
  
}
\details{
  Estimating permutation p-value for a given gene. The procedure is setup to 
  take MatrixEQTL compatible results for one gene at a time.
  For a given gene, multiple bootstraps with adjusted effect sizes are produced 
  and used to draw a relationship between minimum p-values and permutation p-values.
  Since, the relationship follows logistic distribution we use glm to fit the results.
  We also include linear model fit for comparizon.
}
\value{
 The function returns a list with several objects:
 \itemize{
   \item{summ - eigenMT style object for the minimum SNP, which in addition to 
   typical MatrixEQTL object contains effective number of SNPs we get after adjusting inital guess.
   This is still not the final number of effective tests, that can be computed comparing estimated
   permutation p-value eigenMT$pred.permGLM or eigenMT$pred.permLM}
   \item{tim - time it took to run estimation}
   \item{vals - a matrix with corresponding minimum p-values and permutation p-values 
   for the adjusted bootstraps, can be used for diagnostic showing that relationship 
   fits expected pattern and to fit other types of regressions}
   \item{mEQTL - MatrixEQTL object for original gene fit}
   \item{min.snp - min SNP genotypes}
 }
}

%\examples{
%permEst = list(snpM=as.matrix(g.ini[,subs]),
%geneM=as.matrix(exprj[,subs]),
%cvrtM=as.matrix(covar[subs,]),
%snpspos=snpspos,
%genepos=genepos,
%outpf=sprintf("%s_mEQTL_unr.txt", rownames(exprj)[1]),
%pvOutputThreshold=1e-300,
%pvOutputThreshold.csv=1,
%cisDist=1e9,#we have already preprocessed SNPs
%effNumGuess=nrow(g.ini)/4,
%verbose=FALSE, pvalue.hist=FALSE,
%min.pv.by.genesnp = FALSE,
%noFDRsaveMemory=FALSE,
%outdir="unreduced")
%#updNtests=sprintf("%s_updtests.csv", rownames(exprj)[1])
%me = getPermP(permEst, n.perm=100, ini.perm=25)
%}

\author{ Vasyl Zhabotynsky, Wei Sun weisun@email.unc.edu }
\keyword{ methods }

