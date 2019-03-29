# compare with RASQUAL
compare the methods of TReCASE vs. RASQUAL  

the files in this folder is under active updating and please do not use them for any anlaysis yet. 

## to-do's
# compare with RASQUAL
1. compare the methods of TReCASE vs. RASQUAL. In the supplementary materials in google drive, add explicit definition of beta-binomialand negative binomials in TReCASE, i.e., define which parameter is over-dispersion parameter. Also explicitly define the two distributions used by RASQUAL. Make sure their assumption of over-dispersion parameters of NB and BB are the same is equivalent to our definition. 

2. Finish the code check_beta_binomial_MLE.R in GitHub. Calculate Fisher's information matrix and summarize SE and SEE of \pi estimate. here \pi is the proportion of scucess. You can refer to page 3 of Paul et al. paper in google drive for Fisher information matirx. SE is standard error from empirical distribution of \pi estimates. SEE is an estimate based on Fisher information matrix. 


## to-do's

1. Simulate data
- added a clean simulation that allows to simulate panel by panel for panels (a)-(c) of Figure 4 (using 500 simulation 1 panel can be simulated in 15 minutes)
- added a simulation the same simulation setup, but with a subset of profiles matching figure 5
Simulation setup:
mean total read count 100

gene level proportion of allele-specific counts 10%

between sample NB and BB over-dispersion are provided (within sample is assumed to be binomial at simulation step)

for multiple SNPs allele specific reads are uniformely split

2. Regenerate Figure 4b for type I error evaluation using n=64, and a small number of replicates. 
- regenerated panels (a), (b) and (c) for figure 4 for 500 simulations and 

3. Using the same data, regenerate Figure 5(c) to evaluate the estimation of over-dispersion. 
- regenerated panel 5 (c) using RASQUAL style TReCASE model simulated above (generating such RASQUAL itself as shown in figure 5 (c) is quite longer) Also plot the distribution of estimated over-dispersion. 

4. Using the same data, regenerate Figure 6(f) to evaluate the estimate of estimation of genetic effect, for the situation that true effect is 0, draw the distribution of estimated genetic effect. 

5. Use the resutls of 500 replicates, to study what factors is associated with small p-values under the null. For example, a scatter plot of over-dispersion estiamte vs. -log10(p-value). 

6. Confirm the defintion of feature SNP. A SNP is a feature SNP if it is heterozygous in at least one sample? then the number of feature SNPs per gene should increase as sample size increases. 

7. Linear model fit to associate y=log10(TReCASE p-value)-log10(RASQUAL p-value) with all factors. Summarize the relation between y and each factor by scatter plot.

## to-do's with lower priority

1. share the code for RQSQUAL (and TReCASE) analysis in the GitHub.

2. compare resutls of matrix eQTL and TReCASE in asSeq2 in real data. Use 200kb is fine. 

