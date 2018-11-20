\name{prepareBAM}
\alias{prepareBAM}
\title{
  Prepare BAM files for the following counting steps
}
\description{
 Three steps are involved in the preparation. Sorting the bam file by read name, 
 delete those reads that are mapped to multiple places in the genome, and filter 
 out those reads of low quality. 
}
\usage{
prepareBAM(input, outputTag, sortIt = TRUE, getUniqMapping = TRUE, fixmate = FALSE, 
filterIt = TRUE, sort.nmax = NULL, sort.mem = NULL, min.avgQ = 10, 
min.mapQ = 10, flag2kp = 0, flag2rm = 1796, phred = 33, 
keepSortedBAM = FALSE, keepUniqBAM = FALSE)
}
\arguments{
  \item{input}{
  input file name, should be a bam file
}
  \item{outputTag}{
  a tag of the output file name. For example, if outoutTag is "sample1", 
  the output file after sorting is "sample1_sorted_by_name.bam", and the output 
  file after sorting, taking unique reads, and filtering is 
  "sample1_sorted_by_name_uniq_filtered.bam". 
}
  \item{sortIt}{
  let prepareBAM to sort the bam file by read name or not.
  WARNING: Only change the default if you are certain that it is already sorted by read name.
  If the input bam file is not
  sorted by name, this step is required, otherwise the counting of number of 
  reads per exon or per gene will not be correct. 
}
  \item{fixmate}{
  to fix the isPaired and isMateMapped flags for paired reads with same qnames.
  Only works when \code{getUniqMapping=TRUE}.
}
  \item{getUniqMapping}{
  to delete those non-uniquely mapped reads or not. It is strongly recommended 
  to delete those reads otherwise the counting may not be accurate.
}
  \item{filterIt}{
  to filter out low quality reads or not, based on the criterion specified by 
  parameters \code{min.avgQ}, \code{min.mapQ}, \code{flag2kp}, \code{flag2rm}, 
  and \code{phred}. If filterIt is FALSE, these parameters will not be used. 
}
  \item{sort.nmax}{
  sorting is done by splitting the bam files into a few smaller files and 
  sorting those smaller files separately. sort.nmax specify the maximum number 
  of such smaller files. 
}
  \item{sort.mem}{
  sort.mem specify the maximum amount of memory (Mb) used by the sorting step.
  The default value is 1024 (Mb). 
}
  \item{min.avgQ}{
  Each base pair of a read has a sequencing quality score. A read is 
  discarded if the average of those base-pair level quality scores 
  is smaller than min.avgQ
}
  \item{min.mapQ}{
  if the mapping quality score of a is smaller than min.mapQ, this 
  read is disgarded
}
\item{flag2kp}{ a bitwise flag that specifies the reads to keep. By default, 
  flag2kp=0, means to keep all the reads, except those excluded by flag2rm. 
  For example, if flag2kp = 67 = 1 + 2 + 64, we keep the reads that is paired 
  in sequencing (1), mapped in a proper pair (2), and only the first read in 
  the pair (64). }
  
\item{flag2rm}{ By default, flag2rm = 1796 = 1024 + 512 + 256 + 4, 
  i.e., the query sequence itself is unmapped (4) or the alignment is not 
  primary (256), or it fails quality checks (512), or it is PCR duplicates 
  (1024). }

  \item{phred}{
  baseline phred score. Use 33 for Illumina 1.8+, and use 64 for sequences from 
  previous Illumina pipeline. see http://en.wikipedia.org/wiki/FASTQ_format for 
  more details. }
  
  \item{keepSortedBAM}{
  to keep the sorted (by name) bam files or not. By default the sorted bam files
  will be deleted if there are following steps o obtaining unique reads and/or
  filtering. 
}
  \item{keepUniqBAM}{
  to keep the bam files after removing non-uniquely mapped reads or not. By 
  default, such bam files are removed if there is a following step of filtering}
}

\value{
the results will be output as a bam file. 
}

\author{ Wei Sun weisun@email.unc.edu and Vasyl Zhabotynsky}
\keyword{ utilities }
