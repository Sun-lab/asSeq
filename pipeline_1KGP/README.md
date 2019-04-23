# pipeline for eQTL analysis using 1000 Genome Project (1KGP) data. 
For real data analysis I provided adding real data we provide data snippets for two genes that can go throught the analysis

1. Using mapped files and data from TRECASE_MLE process according to RASQUAL approach using GATK

2. Recreate VCF file with these GATK counts

3. Fit RASQUAL

permuted data analysis would have one fix correction to the second step: - permuting SNP information in VCF file

Details to each step:

1. This step requires several steps fixing the flags and adding tags for GATK to be used
a. ensure that VCF file is tab delimited, ordered and has unique snp positions
b. fig mate status of the flags
c. reorder bam file chromosom info to match reference FASTA file
d. add group tag to each read - would be sample ID
e. perform counting with options 
-U ALLOW_N_CIGAR_READS - allowing N's
--minDepth 1 - outputing all the locations with at least one count (can filter out later)
--minMappingQuality 10 - mapping quality at least 10
--minBaseQuality 20 - only consider bases with quality at least 20