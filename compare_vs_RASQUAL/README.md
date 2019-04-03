# compare with RASQUAL
compare the methods of TReCASE vs. RASQUAL  

the files in this folder is under active updating and please do not use them for any anlaysis yet. 

## to-do's
# compare with RASQUAL
1. compare the methods of TReCASE vs. RASQUAL. In the supplementary materials in google drive, add explicit definition of beta-binomialand negative binomials in TReCASE, i.e., define which parameter is over-dispersion parameter. Also explicitly define the two distributions used by RASQUAL. Make sure their assumption of over-dispersion parameters of NB and BB are the same is equivalent to our definition. 
Negative-Binomial is defined in Supplement Figure 1 (eQTL_2019) and top of page 40 of Supplement (Kumasaka paper). 
- in TReCASE for each sample over-dispersion is set to be 
1/theta_i=theta
- in RASQUAL, assuming a set of sample level offsets K_i and genotype-dependent function Q_i which is 2(1-pi), 0.5 and 2pi for G_i=0, 1 and 2
theta_i=theta*K_i*G_i
Beta-Binomial described in Supplement Figure 2 (eQTL_2019) and bottom of page 41 of Supplement (Kumasaka paper).
- In TReCASE for each sample over-dispersion theta_i is alpha+beta=1/theta for all samples with alpha=pi/theta
- In RASQUAL (alpha+beta)=theta*h*K_i*G_i (and on page 42 they set h=1 for allele-specific counts) = theta*K_i*G_i

So, in RASQUAL theta*K_i*G_i is set to be the same in total and allele specific counts. 
It is defined somewhat different from TReCASE though, where theta is defined independently of offset K_i and G_i


2. Finish the code check_beta_binomial_MLE.R in GitHub. Calculate Fisher's information matrix and summarize SE and SEE of \pi estimate. here \pi is the proportion of scucess. You can refer to page 3 of Paul et al. paper in google drive for Fisher information matirx. SE is standard error from empirical distribution of \pi estimates. SEE is an estimate based on Fisher information matrix. 

Note: didn't get a good closed form example why using a splitted sample gets wrong hessian.
a) Got two simulations (for mean number allele-specific reads 10 and 50 in unsplit data) showing the differnce in observed and model based variation
b) Got a small illustration for behavior of hessian coefficients (under pi=0.5 and od=0.5 - behaved similarly at a few other combinations tried):
at pi=0.5 covariance off-diagonal parts of hessian are 0, so we can purely concentrate on diagonal values corresponding to OD and eQTL, for them assuming that reads are equivalently split per SNP. We can see that for a sufficiently large number of reads per gene splitting them into multiple SNPs produces higher diagonal values of hessian matrix. For over-dispersion 0.5 it is true as long as there is more than 1 allele-specific read, for smaller over-dispersion (such as 0.1) about 5 allele-specific reads produce higher value for OD (and any value for eQTL)

## to-do's

1. Simulate data
- added a clean simulation that allows to simulate panel by panel for panels (a)-(c) of Figure 4 (using 500 simulation 1 panel can be simulated in 15 minutes)
- added a simulation the same simulation setup, but with a subset of profiles matching figure 5

Simulation setup:
(a) mean total read count 100
(b) gene level proportion of allele-specific counts 10%
(c) between sample NB and BB over-dispersion are provided (within sample is assumed to be binomial at simulation step)
(d) for multiple SNPs allele specific reads are uniformely split

2. Regenerate Figure 4b for type I error evaluation using n=64, and a small number of replicates. 
- regenerated panels (a), (b) and (c) for figure 4 for 500 simulations and 

3. Using the same data, regenerate Figure 5(c) to evaluate the estimation of over-dispersion. 
- regenerated panel 5 (c) using RASQUAL style TReCASE model simulated above (generating such RASQUAL itself as shown in figure 5 (c) is quite longer) Also plot the distribution of estimated over-dispersion. 

4. Using the same data, regenerate Figure 6(f) to evaluate the estimate of estimation of genetic effect, for the situation that true effect is 0, draw the distribution of estimated genetic effect. 
- added figure 6f
(regarding distribution of estimated gentetic effects I think the other plot with hessian matrix estimates is more of interest)

5. Use the resutls of 500 replicates, to study what factors is associated with small p-values under the null. For example, a scatter plot of over-dispersion estiamte  vs. -log10(p-value). 
plotted vs 
- OD(NB)
- OD(BB)
- b0

6. Confirm the defintion of feature SNP. A SNP is a feature SNP if it is heterozygous in at least one sample? then the number of feature SNPs per gene should increase as sample size increases. 
- there are some cutoffs for fSNPs which kept during the analysis: --min-coverage-depth 0.05

7. Linear model fit to associate y=log10(TReCASE p-value)-log10(RASQUAL p-value) with all factors. Summarize the relation between y and each factor by scatter plot.

## to-do's with lower priority

1. share the code for RQSQUAL (and TReCASE) analysis in the GitHub.

2. compare resutls of matrix eQTL and TReCASE in asSeq2 in real data. Use 200kb is fine. 

