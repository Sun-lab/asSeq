# compare with RASQUAL
compare the methods of TReCASE vs. RASQUAL  

the files in this folder is under active updating and please do not use them for any anlaysis yet. 

## to-do's
# compare with RASQUAL
compare the methods of TReCASE vs. RASQUAL  

the files in this folder is under active updating and please do not use them for any anlaysis yet. 

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

<<<<<<< HEAD
1. Simulate data
- added a clean simulation that allows to simulate panel by panel for panels (a)-(c) of Figure 4 (using 500 simulation 1 panel can be simulated in 15 minutes)
- added a simulation the same simulation setup, but with a subset of profiles matching figure 5
Simulation setup:
mean total read count 100
gene level proportion of allele-specific counts 10%
between sample NB and BB over-dispersion are provided (within sample is assumed to be binomial at simulation step)
for multiple SNPs allele specific reads are uniformely split
 
=======
1. Simulate data, run analysis using two methods, TReCASE and TReCASE-RL (RASQUAL-like), which set the over-dispersion of NB and BB to be the same and model the ASReC of each fSNP by a beta-binomial distribution. 
>>>>>>> 09251af78797e63e5b309ec377a4dff1dc44c190
2. Regenerate Figure 4b for type I error evaluation using n=64, and a small number of replicates. 
- regenerated panels (a), (b) and (c) for figure 4 for 500 simulations and 

3. Using the same data, regenerate Figure 5(c) to evaluate the estimation of over-dispersion. 
- regenerated panel 5 (c) using RASQUAL style TReCASE model simulated above (generating such RASQUAL itself as shown in figure 5 (c) is quite longer)

4. Using the same data, regenerate Figure 6(d) to evaluate the estimate of estimation of genetic effect, only consider the situation that true effect is 0. 
?
