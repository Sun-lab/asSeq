Note you will need to ensure that tabix and path to RASQUAL are addedto the environment
module add tabix;module add r;
export PATH=/nas/longleaf/home/zhabotyn/progs/rasqual-master/bin/:$PATH


We simulate the data in format acceptable by TReCASE style - one line corresponds to one gene.
If there is n fSNPs, then number of columns of genotype classes and allele-specific counts would have n*#Individuals
RASQUAL, however, has #Individual columns for total and allele-specific counts, with several fSNPs recorded at appropriate positions.
One gene, thus, can use several linesfrom VCF file and one line from total read counts file.
Consequently, the first step, is to reformat input into VCF format by taking appropriate parts of allele-specific count matrix and write it in VCF style.

Then, we save several files that are needed by RASQUAL for input:

X.bin - design matrix for covariates (X_suffix.bin)

K.bin - offset matrix (K_suffix.bin)

Y.bin - total read counts matrix (tot_suffix.bin)

and 

VCF file - containing fSNPs and tSNPs (hap_suffix.vcf)

we use suffix because we may run simulations in parallel, so simulation for each profile should have distinctive file.

X, K and Y are provided using options -x, -k and -y

whereas VCF file is provided as stdin - hence we use tabix with selecting a range of candidate SNPs (to model appropriate window size)

Additionally among required options are the options:

a. -s and -e giving lists of starting and ending positions of the exons (in this simulation for simplicity we just use one simulated exon)
According to these parameters RASQUAL will choose which SNPs are fSNPs.

b. Option -j is required to point which row of Y(total read counts) matrix to start processing.

c. -l and -m - provide number of tSNPs and fSNPs for this gene (for simplicity we provide an upper bound which uses somewhat more memory, but otherwise is not very important)

d. -n - #Individuals (sample size)

In this simulation we also ask to fit only 1 SNP using option --rsnp (and added a parameter that forces returning only one SNP test -z) to improve timing.

Finally, we add three options: --fix-delta --fix-phi (since we don't have reference and mapping bias in this simulations and --min-coverage-depth 0.0 to use all fSNPs in the simulation.

We provide three blocks:

1. fit joint TReC+ASReC model
2. fit only TReC part of the model
3. fit only ASReC part of the model

Resulting values will be stored in 
rasq/fin_suffix.csv, rasq/fin_suffix_tot.csv and rasq/fin_suffix_ase.csv