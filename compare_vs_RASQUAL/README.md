# compare with RASQUAL
compare the methods of TReCASE vs. RASQUAL  

the files in this folder is under active updating and please do not use them for any anlaysis yet. 

## to-do's 

1. Add code to run RASQUAL for simulaiton 

- added this code to rasqual_simu folder

2. Add simulation code to simulate data with beta-binomial variation within individual. 

- updated the package to allow simuX function to simulate within individual over-dispersion provided as an extra parameter (old code should still work)

3. Add the code to real data anlaysis, such as orginal data and permutated data in folder asSeq/pipeline_1KGP/

- finishing uploading real data analysis using RASQUAL

4. Add description of simulaiton setup in Supplementary materails. Refine the writting in simulation sections, particular the last two sections. 

- added some simulation setup and extended next sections, will look at it again and add one more iteration

5. For figures of Fisher information matrix, add the code to generate data, modify the figures to replace “as-cnt” with ASReC, and for the panel with OD as title while you are actually looking at rho, you can use OD (rho) instead of OD.
- changed, also changed "hessian" label to "information"

6. Compare matrixEQTL vs. TReCASE in asSeq2 in real data. Use 200kb is fine.  For real data anlaysis pipeline, we may consider a screening step using matrixEQTL. We do not want to run TReCASE for every gene and every SNP since it is time consuming. So first run matrixEQTL and only keep the asossocaitons with p-value < 0.1, and run TReCASE for those cases. To evaluate whehter such screening step makes sense, we will actualy run matirxEQTL and TReCASE for all gene,SNP pairs and compare results. When running matrixEQTL, first quantile normalize gene expression gene by gene, and thus matrixEQTL results is somewhat a rubust version of linear regression that use data quantiles rather than the raw data. Use the code in asSeq/pipeline_1KGP/quantile_transformation.R. We will see sometime TReCASE give significant association, but matrixEQTL does not. This is often the case with outlier. So we should discard those associaitons. 

## to-do's (mostly finished)

1. Compare the density of NB and BB used by TReCASE and RASQUAL. 

2. Added fisher directory in GitHub: it contains the code to calculate Fisher's information matrix and summarizes SE and SEE of \pi and \rho estimate. (\pi is the proportion of success and \rho is a version of over-dispersion parameter connected to \theta=\rho/(1-\rho)). 
 Paul et al. paper was used for Fisher information matirx. 
- In this simulation we concentrate only on allele-specific counts since inflation of type 1 error was observed in splitting individuals into multiple-SNPs.
- Add one more code to reproduce the results. 

2.1 I added another look at hessian by plotting how hessian component for both OD (parameterized by \rho) and eQTL effect are behaving depending on number of allele-specific counts in a null hypothesis scenario (\pi=0.5), plotted in two scenarios of \theta=0.1 and 0.5

3. Simulate data
- added a clean simulation that allows to simulate panel by panel for panels (a)-(c) of Figure 4 (using 500 simulation 1 panel can be simulated in 15 minutes)
- added a simulation the same simulation setup, but with a subset of profiles matching figure 5

Simulation setup:
- (a) mean total read count 100
- (b) gene level proportion of allele-specific counts 10%
- (c) between sample NB and BB over-dispersion are provided (within sample is assumed to be binomial at simulation step)
- (d) for multiple SNPs allele specific reads are uniformely split

4. Regenerate Figure 4b for type I error evaluation using n=64, and a small number of replicates. 
- regenerated panels (a), (b) and (c) for figure 4 for 500 simulations and 

5. Using the same data, regenerate Figure 5(c) to evaluate the estimation of over-dispersion. 
- regenerated panel 5 (c) using RASQUAL style TReCASE model simulated above (generating such RASQUAL itself as shown in figure 5 (c) is quite longer) 

6. Using the same data, regenerate Figure 6(f) to evaluate the estimate of estimation of genetic effect, for the situation that true effect is 0, draw the distribution of estimated genetic effect. 
- added figure 6f
(regarding distribution of estimated over-dispersion and genetic effects I think the other plot with hessian matrix estimates along with distribution of estimated effects is more of interest)

7. Confirm the defintion of feature SNP. A SNP is a feature SNP if it is heterozygous in at least one sample? then the number of feature SNPs per gene should increase as sample size increases. 
- there are some cutoffs for fSNPs which kept during the analysis: --min-coverage-depth 0.05
- figures put in fsnps folder
- need to check whether further how it filters out some of fSNPs

8. Linear model fit to associate y=log10(TReCASE p-value)-log10(RASQUAL p-value) with all factors. Summarize the relation between y and each factor by scatter plot.
- this model fit is added to the subsection in supplement on google drive.



