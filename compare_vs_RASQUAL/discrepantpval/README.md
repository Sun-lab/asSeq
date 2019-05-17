1. This folder compares discrepant results from TReCASE and RASQUAL model.
significant_combined_spline.png show dependency of fraction of significant genes (at a cutoff listed in the legend) given that in other model they are significant (panels (a) and (b)) and given that in other model they are not significant (panels (c) and (d))
- The results are plotted vs number of fSNPs, spline-smoothed.

2. added markdown for model fit (and required input file): linear_model_cov.csv and markdown.Rmd
or plain r code in fit_lm.R

Columns of linear_model_cov.csv file are:

- fID - gene id
- sID - snp id
- SNPchr - snp chormosome
- SNPpos - snp position
- Ref - reference allele
- Alt - alternative allele
- AF - allele frequency
- HWEChisq - Hardy Weinberg Equilibrium Chi2
- impQual - imputation quality
- qval - log10 RASQUAL q-value
- Chi2 - Chi2 value
- Pi - RASQUAL estimate of eQTL effect
- DeltaErr - mapping error estimate
- PhiBias - reference bias estimate
- OD - RASQUAL over-dispersion (larger value - smaller over-dispersion)
- corfSNP - correlation of original fSNP and RASQUAL corrected
- corrSNP - same for rSNPs
- qvalT - TReCASE q-value
- pvalT - TReCASE p-value
- pvalR - RASQUAL p-value
- odNB - TReCASE Negative-Binomial over-dispersion
- odBB - TReCASE Beta-Binomial over-dispersion     
- mpp - median p-value for RASQUAL model refitted on permuted SNPs
- asT - allele-specific expression for 280 samples using TReCASE approach
- asR - allele-specific expression for 280 samples using RASQUAL approach


3. updated markdown file: 
- try to replace the smoothScatter with scatter plot with approprate point size and color. 
- Draw some scatter plot of covariates versus the response variable for those ~220 genes. 
- added plots for covariates, 
- added model fits at several scanario cutoffs (all genes, more discrepant and less discrepant genes) to see the robustness - NB OD always remains significantly associated with stronger TReCASE p-values
- cleaned up the code and split into smaller chunks