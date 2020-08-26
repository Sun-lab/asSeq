# pipeline for eQTL analysis using GTEx data.
Step by step code pipeline is in code subfolder In short: 

step 0 - trims suspect counts

step 1 - the data and reformat it to separate genes format in order to paralellize computatoins

step 2 -  basic eQTL fit using TReCASE (not permutation corrected)

step 4 - estimating effective number of SNPs in order to estimate permutation p-value using bootstraps of MatrixEQTL

step 5 - collecting gene level results and recoding counts to be used for dynamic eQTL study

step 5a - Principal Component analysis on allele-specific counts 

step 6 - analysis of dynamic eQTLs for three conditions of interest: age, CTCF and TP53 expression
- nocov - model with no other covariates
- gPC - model including top 2 genotype Principal Components
- gPC_PF - model including top 2 genotype PCs and top 5 PEER factors

step 8 - pathway analsys for the findings from the step 6

step 9 - summarizing singificant findings for all tissues

example - a specific example that can be used to run steps 2 and 4 to produce permuted p-values for a particular gene
