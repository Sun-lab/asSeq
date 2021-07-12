\name{ase}
\alias{ase}
\title{
Allele Specific Expression (ASE) model
}
\description{
  eQTL mapping using imbalance of allele-specific read counts. 
  }
\usage{
ase(Y1, Y2, Z, output.tag, p.cut, min.AS.reads=5, min.AS.sample=5, 
  min.n.het=5, local.only = TRUE, local.distance = 2e+05, 
  eChr = NULL, ePos = NULL, mChr = NULL, mPos = NULL, trace = 0)
}
\arguments{
  \item{Y1}{ matrix of gene expression from one allele (haplotype). 
    Each row corresponds to one sample and each column corresponds to one gene. 
  }
  
  \item{Y2}{ matrix of gene expression from the other allele (haplotype). 
    Each row corresponds to one sample and each column corresponds to one gene. 
    Must be the same dimension as Y1.
  }
  
\item{Z}{ matrix of phased genotype data, Each row corresponds to one sample 
    and each column corresponds to one SNP. Must take value of 0 (AA), 1 (AB), 
    3 (BA) or 4(BB), where A and B are two alleles of the SNP. AB and BA are 
    different since AB means the 1st/2nd haplotype harbors the A/B allele, 
    respectively, while BA means the 1st/2nd haplotype harbors the B/A allele, 
    respectively. 
  }
  
  \item{output.tag}{
    tag of output file. The actual output files will be tag_eqtl.txt and 
    tag_freq.txt.
  }

  \item{p.cut}{ p-value cut-off, only the associations with p-value smaller than 
    p.cut are saved. 
  }

  \item{min.AS.reads, min.AS.sample}{ one gene will be considered for 
  allele-specific association analysis only if at least min.AS.sample samples
  have at least min.AS.reads allele-specific reads.
  }
  
  \item{min.n.het}{ an allele-specific eQTL will be evaluated only if the 
  at least min.n.het samples have heterozygous genotype at the query SNP. 
  }
  
  \item{local.only}{ whether examine local eQTL only. }
  
  \item{local.distance}{ if the distance between a gene and a genetic marker
   is smaller than local.distance, we call it as a local eQTL. }
   
  \item{eChr}{ chromosomes of gene expression traits, must be
   a vector of integers, which means sex chromosomes (X, Y) need
   to be coded as numbers }
   
  \item{ePos}{ positions of gene expression traits, must be
   a vector of integers}
   
  \item{mChr}{ chromosomes of markers, must be a vector of integers,
   which means sex chromosomes (X, Y) need to be coded as numbers }
   
  \item{mPos}{ positions of markers, must be a vector of integers}

  \item{trace}{ If trace=0, do not output tracing information.  
  if trace = 1, output brief tracing information, and 
  if trace > 1. output detailed tracing information. }

}

\value{
The function return 1 if no error occurs.

 The results are written into two output files: 
 
 one file (output.tag_eqtl.txt) stores the eQTL mapping result that pass the p-value 
 threshold p.cut, including gene ID, marker ID, Chi-square statistic, 
 and p-value; 
 
 the other file (output.tag_freq.txt) includes a frequency table of all the p-values.
 }
\author{ Wei Sun weisun@email.unc.edu }
\keyword{ methods }

