# compare with RASQUAL
compare the methods of TReCASE vs. RASQUAL  

the files in this folder is under active updating and please do not use them for any anlaysis yet. 

## to-do's
# compare with RASQUAL
1. In the supplementary materials in google drive, added explicit definition of beta-binomialand negative binomials in TReCASE, including definition which parameter is over-dispersion parameter. Also explicitly defined the two distributions used by RASQUAL. Described similarities and discrepancies in over-dispersion assumptions (in the google drive supplement)
- added to the supplementary materials


2. Added fisher directory in GitHub: it contains the code to calculate Fisher's information matrix and summarizes SE and SEE of \pi and \rho estimate. (\pi is the proportion of success and \rho is a version of over-dispersion parameter connected to \phi=\rho/(1-\rho)). 
 Paul et al. paper was used for Fisher information matirx. 
- In this simulation we concentrate only on allele-specific counts since inflation of type 1 error was observed in splitting individuals into multiple-SNPs.
- SSE_vs_modelSE_an_100.png - observed vs model based standard error for simulation assuming on average 10 allele-specific reads per individual
(I added two more for 20 and 40 reads - files anding with _200 and _400)
- It can be seen that for eQTL effect model based variation decreases when reads are split in multiple SNPs
(for OD parameter estimate the picture is similar as long as counts get high enough - for low counts direction can be opposite)
- SSE_vs_modelSE_an_100.png with a brief description is added to the supplement in google drive

2.1 I added another look at hessian by plotting how hessian component for both OD and eQTL effect are behaving depending on number of allele-specific counts in a null hypothesis scenario (\pi=0.5), plotted in two scenarios of \phi=0.1 and 0.5
a) hess_relation1_0.1.png - plotting contribution to hessian
hess_logscale1_0.1.png - plotting contribution to hessian on log scale
It confirms observed above that splitting in multiple SNPs always increases for eQTL parameter and for large enough AS counts in OD parameter.
b) Since you mentioned plotting it in terms of standard deviation I also added the versions of the plot in terms of SD:
sd_by_asc_1ind_0.1.png
sd_by_asc_1ind_0.5.png
I also added relative SD(multi-SNP vs 1-SNP) to illustrate that the ratio even gets bigger for larger counts
relSD_0.1.png & relSD_0.1.png
c) These plots are added with a brief description to the supplement in google drive

3. Simulate data
- added a clean simulation that allows to simulate panel by panel for panels (a)-(c) of Figure 4 (using 500 simulation 1 panel can be simulated in 15 minutes)
- added a simulation the same simulation setup, but with a subset of profiles matching figure 5

Simulation setup:
(a) mean total read count 100
(b) gene level proportion of allele-specific counts 10%
(c) between sample NB and BB over-dispersion are provided (within sample is assumed to be binomial at simulation step)
(d) for multiple SNPs allele specific reads are uniformely split

4. Regenerate Figure 4b for type I error evaluation using n=64, and a small number of replicates. 
- regenerated panels (a), (b) and (c) for figure 4 for 500 simulations and 

5. Using the same data, regenerate Figure 5(c) to evaluate the estimation of over-dispersion. 
- regenerated panel 5 (c) using RASQUAL style TReCASE model simulated above (generating such RASQUAL itself as shown in figure 5 (c) is quite longer) Also plot the distribution of estimated over-dispersion. 

6. Using the same data, regenerate Figure 6(f) to evaluate the estimate of estimation of genetic effect, for the situation that true effect is 0, draw the distribution of estimated genetic effect. 
- added figure 6f
(regarding distribution of estimated gentetic effects I think the other plot with hessian matrix estimates is more of interest)

7. Use the resutls of 500 replicates, to study what factors is associated with small p-values under the null. For example, a scatter plot of over-dispersion estiamte  vs. -log10(p-value). 
figures put in figures subfolder

8. Confirm the defintion of feature SNP. A SNP is a feature SNP if it is heterozygous in at least one sample? then the number of feature SNPs per gene should increase as sample size increases. 
- there are some cutoffs for fSNPs which kept during the analysis: --min-coverage-depth 0.05
figures put in fsnps folder
need to check whether further how it filters out some of fSNPs

7. Linear model fit to associate y=log10(TReCASE p-value)-log10(RASQUAL p-value) with all factors. Summarize the relation between y and each factor by scatter plot.
this model fit is added to the subsection in supplement on google drive.


## to-do's with lower priority

1. share the code for RQSQUAL (and TReCASE) analysis in the GitHub.

2. compare resutls of matrix eQTL and TReCASE in asSeq2 in real data. Use 200kb is fine. 

