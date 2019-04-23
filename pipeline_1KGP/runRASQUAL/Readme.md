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

note, java version may be lead to gatk not recognizing ASEReadCounter

code worked with OpenJDK Runtime Environment (build 1.8.0) and didn't work with Java(TM) SE Runtime Environment 18.3 (build 10.0.2) giving a message

"Invalid command line: Malformed walker argument: Could not find walker with name: ASEReadCounter"

2. create VCF file with allele specific counts

a. count number of SNPs within a window - will be used as input for RASQUAL

b. add snp level allele-specific counts from individual level files produced in step 1 to a vcf (stored in a separate folder)
(for permutation data the modification is - petrurbing genotype data at each SNP)

c. save in binary format other required files if they don't exist yet: design matrix (X_ex.bin), matrix of offsets (K_ex.bin) and total read counts (Tcnt_ex.bin)

3. run RASQUAL

Includes several mandatory options: 

a. it takes part of a VCF file as stdin (using tabix with a certain range)

b. -y - total read counts (Y matrix)

c. -k - offsets (K matrix)

d. -x - design (X matrix)

e. -p - number of columns in X-matrix

f. -n - number of samples

g. -j - row in Y-matrix corresponding to the jene that is to be tested

h. -l, -m - technical parameters to ensure that enough space for fSNPs and rSNPs is dedicated (can just put number of rSNPs)

i. -s, -e - starting and ending points for exons of the gene (in order for a program to learn which SNPs to consider fSNPs with allele-specific counts)

j. -z - convert genome imputation quality score (R^2 or I^2) into allelic probability (AP)

k. --force - force running the code even if it takes many SNPs.

l. --posterior-genotype - do posterior genotype update
  rasBASE = sprintf("rasqual -y %s -k %s -x %s -p %s -n %s", totb, kbin, xbin, np, nsam)
  rasPOSj = sprintf("-j %s -l %s -m %s -s %s -e %s", gnj, nl, nm, stj, enj)
  rasOPT = "-z --force --posterior-genotype"# --population-only --min-coverage-depth 0.0 --as-only 
