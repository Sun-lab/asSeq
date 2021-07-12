\name{quantifyAsReads}
\alias{quantifyAsReads}
\title{ Quantifies number of allele-specific reads along with positions of SNPs used  }
\description{
  A sequence read is counted as allele specific only if it harbors one or more
  heterozygous genetic markers, and their genotypes are consistent with one 
  and only one of the two haplotypes specified in file snpList.
}
\usage{
quantifyAsReads(input, snpList, outputTag=NULL, flag2kp=0, flag2rm=3844, 
         prop.cut=.5, min.avgQ=10, min.snpQ=10, phred=33, skip=0)
}
\arguments{
  \item{input}{ input file name, this should be file with BAM format, 
  e.g., the output file of Tophat }
  
  \item{snpList}{ name of the file including a SNP list, must be tab-delimited 
  file with four columns, chromosome, position, allele 1 and allele 2, without 
  header. }

  \item{outputTag}{ output files will be named as outputTag_hap1.bam,  
  outputTag_hap2.bam and outputTag_hapN.bam }
  
  \item{flag2kp}{ a bitwise flag that specifies the reads to keep. By default, 
  flag2kp=0, means to keep all the reads, except those excluded by flagrm. 
  For example, if flag2kp = 67 = 1 + 2 + 64, we keep the reads that is paired in 
  sequencing (1), mapped in a proper pair (2), and only the first read in 
  the pair (64). }
  
  \item{flag2rm}{ By default, flag2rm = 3844 = 2048 + 1024 + 512 + 256 + 4, 
  i.e., the query sequence itself is unmapped (4) or the alignment is not 
  primary (256), or it fails quality checks (512), or it is PCR duplicates 
  (1024), or it is supplementary alignment (2048). }
  
  \item{prop.cut}{ assign a read to one haplotype if the proportion of 
  informative SNPs suggesting it is from one haplotype is larger than prop.cut.
  }

  \item{min.avgQ}{ Minimum of average sequencing quality score, the default value 
  is 10}
  
  \item{min.snpQ}{ Minimum sequencing quality allowed for a nucleotide to be used 
  as a SNP, the default value is 10}

  \item{phred}{
  baseline phred score. Use 33 for Illumina 1.8+, and use 64 for sequences from 
  previous Illumina pipeline. see http://en.wikipedia.org/wiki/FASTQ_format for 
  more details. }
  
  \item{skip}{ The number of lines of the data file to skip before beginning to
 read data }
 
}

\value{
The output will be named as outputTag, which links to the file in which quantification is stored.
}
\author{ Wei Sun weisun@email.unc.edu }
\keyword{ utilities }

