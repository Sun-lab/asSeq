permuted data analysis would have one fix correction to the second step: - permuting SNP information in VCF file

Details to each step:

1. This step requires several steps fixing the flags and adding tags for GATK to be used

a. ensure that VCF file is tab delimited, ordered and has unique snp positionsb. fig mate status of the flags

c. reorder bam file chromosom info to match reference FASTA file

d. add group tag to each read - would be sample ID

e. perform counting with options 

-U ALLOW_N_CIGAR_READS - allowing N's

--minDepth 1 - outputing all the locations with at least one count (can filter out later)

--minMappingQuality 10 - mapping quality at least 10

--minBaseQuality 20 - only consider bases with quality at least 20


2. create VCF file with allele specific counts

3. run RASQUAL