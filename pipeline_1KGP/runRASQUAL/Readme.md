permuted data analysis would have one fix correction to the second step: - permuting SNP information in VCF file

Details to each step:

1. step1_use_GATK.R. This step requires several steps fixing the flags and adding tags for GATK to be used
+ a. ensure that VCF file is tab delimited, ordered and has unique snp positions
+ b. fig mate status of the flags
+ c. reorder bam file chromosom info to match reference FASTA file
+ d. add group tag to each read - would be sample ID
+ e. perform counting with options 
```
-U ALLOW_N_CIGAR_READS - allowing N's

--minDepth 1 - outputing all the locations with at least one count (can filter out later)

--minMappingQuality 10 - mapping quality at least 10

--minBaseQuality 20 - only consider bases with quality at least 20
```
note, java version may be lead to gatk not recognizing ASEReadCounter

code worked with OpenJDK Runtime Environment (build 1.8.0) and didn't work with Java(TM) SE Runtime Environment 18.3 (build 10.0.2) giving a message

"Invalid command line: Malformed walker argument: Could not find walker with name: ASEReadCounter"

2. step2_update_VCF.R/step2_permute_VCF.R create VCF file with allele specific counts
+ a. count number of SNPs within a window - will be used as input for RASQUAL
+ b. add snp level allele-specific counts from individual level files produced in step 1 to a vcf (stored in a separate folder)
(for permutation data the modification is - petrurbing genotype data at each SNP)
+ c. save in binary format other required files if they don't exist yet: design matrix (X_ex.bin), matrix of offsets (K_ex.bin) and total read counts (Tcnt_ex.bin)

3. run RASQUAL

Includes several mandatory options: 
+ RASQUAL takes part of a VCF file as stdin (using tabix with a certain range)
```
-y - total read counts (Y matrix)
-k - offsets (K matrix)
-x - design (X matrix)
-p - number of columns in X-matrix
-n - number of samples
-j - row in Y-matrix corresponding to the jene that is to be tested
-l, -m - parameters suggesting to RASQUAL number of fSNPs and rSNPs (can put number of rSNPs for each of parameters)
-s, -e - starting and ending points for exons of the gene (in order for a program to learn which SNPs to consider fSNPs with allele-specific counts)
```
and a few extra options:
```
-z - convert genome imputation quality score (R^2 or I^2) into allelic probability (AP)
--force - force running the code even if it takes many SNPs.
--posterior-genotype - do posterior genotype update
other options to consider (were used in simulations): --min-coverage-depth 0.0, --as-only, --population-only
```

note. initial VCF can be created form imputed data according to pipeline presented in step8 and step9 in:
https://github.com/Sun-lab/eQTL-PoO-analysis-pipeline/tree/master/parents/block2_phasining_imputation
