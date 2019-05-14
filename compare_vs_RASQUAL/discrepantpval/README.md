1. This folder compares discrepant results from TReCASE and RASQUAL model.
significant_combined_spline.png show dependency of fraction of significant genes (at a cutoff listed in the legend) given that in other model they are significant (panels (a) and (b)) and given that in other model they are not significant (panels (c) and (d))
- The results are plotted vs number of fSNPs, spline-smoothed. 

2. added markdown for model fit (and required input file): linear_model_cov.csv and markdown.Rmd
or plain r code in fit_lm.R

3. updated markdown file: 
- added plots for covariates, 
- added model fits at several scanario cutoffs (all genes, more discrepant and less discrepant genes) to see the robustness - NB OD always remains significantly associated with stronger TReCASE p-values
- cleaned up the code and split into smaller chunks