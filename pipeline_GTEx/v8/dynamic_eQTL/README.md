#  Detect dynamic eQTL using ASE

Check whether found allelic imbalance depends on the age, TP53 or CTCF expression

## Models: 

long - model including PEER factors and PCs

short - model that only include the variable of interest, including  age, or the gene expression of TP53 or CTCF


##  ASE counts 

We converted original hap1 and hap2 counts in the following fashion:

taking the most significant SNP classified into 4 groups according to TReCASE model AA, AB, BA and BB

if SNP was from the group 0, 1 or 3 - keep the values as they are

for the group 2 - swap hap1 and hap2 values.

Call these processed values as A1 and A2.

Also we classified whether a SNP for a given individual and gene is heterozygous: 0 - homozygous, 1 - heterozygous. Call these SNPInd

The output file is written in a format: A1|A2|SNPInd. 

Since the most significant SNP could be different in the short and the long models we record a file for each the model.