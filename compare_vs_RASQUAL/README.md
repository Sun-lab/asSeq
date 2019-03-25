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
- regenerated panel 5 (c) using RASQUAL style TReCASE model simulated above (generating such RASQUAL itself as shown in figure 5 (c) is quite longer)

4. Using the same data, regenerate Figure 6(d) to evaluate the estimate of estimation of genetic effect, only consider the situation that true effect is 0. 
?