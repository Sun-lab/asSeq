\name{asSeq2-package}
\alias{asSeq2-package}
\alias{asSeq2}
\docType{package}
\title{
\packageTitle{asSeq2}
}
\description{  
  eQTL mapping using both Total Read Count (TReC) 
  and/or Allele Specific Expression (ASE).
}
\usage{
trecase(Y, Y1 = NULL, Y2 = NULL, Z, XX, SNPloc, geneloc, fam_nb = T,
    file_trec = "trec.txt", file_trecase = "trecase.txt", cis_window = 100000L,
    useASE = 1L, min_ASE_total = 8L, min_nASE = 10L, eps = 5e-05,
    max_iter = 400L, show = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  A matrix of size N-by-G of TReC values for G genes and N subjects. Data is stored as one column per gene and one row per subject.
}
  \item{Y1}{
  A matrix of size N-by-G of read counts mapping to arbitrarily defined allele 1. Data is stored as one column per gene and one row per subject. Y1 is set to NULL when useASE = 0 (only test eQTL effect of total read count).
}
  \item{Y2}{
  A matrix of size N-by-G of read counts mapping to arbitrarily defined allele 1. Data is stored as one column per gene and one row per subject. Y2 is set to NULL when useASE = 0 (only test eQTL effect of total read count).
}
  \item{Z}{
  A matrix of size N-by-S of genotype data for N subjects at S potential SNPs. Data is stored as one row per subject, one column per SNP. A single item of Z takes values 0,1,2, or 3 corresponding to genotypes AA, AB, BA, and BB.
}
  \item{XX}{
  A matrix of size N-by-M of covariates for N subjects and M measured covariates. Data is stored as one column per covariate, one row per subject.
}
  \item{SNPloc}{
  A data.frame of size S-by-3 SNP information for S SNPs and colnames have to be c('snp', chr, 'pos'). 
}
  \item{geneloc}{
  A data.frame of size G-by-3 SNP information for G genes and colnames have to be c('gene', chr, 'start','end'). 
} 
\item{fam_nb}{
  Which regression family to use. 
  0: Poisson regression.
  1: means run negtive binomial regression.
  
}
\item{file_trec}{
  Path to where TReC output file need to be saved.
}
\item{file_trecase}{
  Path to where TReCASE output file need to be saved.
}
  \item{cis_window}{
  Value determining local eQTL. A SNP is considered local to a gene if and only if it is located within the gene body or within cis_window bp from either the start or end of the gene.
}
  \item{useASE}{
  Indicator of whether or not to use ASE data for any gene (1=Use ASE,0=TReC only). Default is set to 1.
}
  \item{min_ASE_Total}{
  The minimum number of Allele Specific reads a subject must have to be included in the ASE analysis. Recommended value of 8.
}
  \item{min_nASE}{
  The minimum number of subjects that must exceed the min_ASE_Total threshold in order for ASE data to be included.
}
  \item{eps}{
  Specifies whether TReCASE model has converged. Default of 5e-5.
}
 \item{max_iter}{
  Maximum number of iterations for TReC or TReCASE model fits. Default value of 400.
}
  \item{show}{ If show=0, do not output tracing information.  
  if show = 1, output tracing information. }

}
\details{
  Joint model of both total expression and allele-specific expression will only. If useASE = 0, only model of total expression will be fitted. When useASE = 1, we only run TReCASE model for the SNP with smallest p-value in each gene.
}
\value{
All output is stored in the file specified in file_trec and/or file_trecase in tab-delimited format. 

file_trec:

GeneRowID: the number corresponds to nth gene in geneloc 

MarkerRowID: the number corresponds to nth SNP in SNPloc  

TReC_b: Estimates of eQTL effect.

TReC_Chisq: Chi-square statistic.

TReC_Pvalue: p-value.

TreC_Conv: 0-TReC model does not coverged, 1-TReC model coverged.

beta: The coeffecients of Negtive Binomial regression.

phi: The overdispersion parameter in Negtive Binomial regression.



file_trecase:

TReC_* means values for TReC model.

Joint_* means values for TReCASE model.

Joint_theta: The overdispersion parameter in Beta-Binomial regression.

Converge: 0: TReCASE model does not coverged, 1: TReCASE model coverged.

CisTrans_Chisq: The value of the test statistic for the Cis-Trans Score Test

CisTrans_Pvalue: Pvalue for the Cis-Trans Score Test

}
\author{ Wei Sun weisun@email.unc.edu}

\examples{
{

  }
}
