\name{KLK1}
\alias{KLK1}
\docType{data}
\title{
  Data for the TReCASE model of one gene KLK1 (ENSG00000167748) 
  vs. one SNP (rs1054713), and together with covaraites data.
}
\description{
  This data object is a list that is to be used to demonstrate function trecase.
}
\usage{data(KLK1)}
\details{
  This list includes 9 components. 
  \itemize{
	\item{Y:}{ The total number of reads of gene KLK1 in the 65 samples. }
	
	\item{Y1:}{ The number of allele-sepecific reads of gene KLK1 in 
	the 65 samples, from the 1st haplotype. }
	
	\item{Y2:}{ The number of allele-sepecific reads of gene KLK1 in 
	the 65 samples, from the 2nd haplotype. }
	
	\item{X:}{ matrix of cofounding covariates. Each row corresponds 
	to one sample and each column corresponds to one gene. }
	
	\item{Z:}{ Genotype of SNP rs1054713. It takes values of 0 (AA), 1 (AB), 
	3 (BA) or 4(BB), where A and B are two alleles of the SNP. AB and BA are 
	different since AB means the 1st/2nd haplotype harbors the A/B allele, 
	respectively, while BA means the 1st/2nd haplotype harbors the B/A allele, 
	respectively. 
	}
  
	\item{eChr:}{ chromosome of KLK1 }
	
	\item{ePos:}{ position of KLK1}
	
	\item{mChr:}{ chromosome SNP rs1054713 }
	
	\item{mPos:}{ positions of SNP rs1054713 }
  }
}
\source{
The genotye data is from HapMap website and the RNA-seq data  
is from Pickrell et al. (2010). Understanding mechanisms
underlying human gene expression variation with RNA sequencing. 
Nature. 464(7289):768-72
}
\references{
	Wei Sun (2011), A Statistical Framework for eQTL Mapping Using RNA-seq Data, 
	Biometrics, in press
}
\examples{
data(KLK1)
}
\keyword{datasets}
