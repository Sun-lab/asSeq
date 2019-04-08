This folder studies issues coming from splitting reads into multiple SNPs.
1. fishercl.R - code for calculating fisher matrix and likelihood, simush_hesscl.R - to run a simulation on 1, 2, 4 and 8 SNPs fitting allele-specific counts only
2. SSE_vs_modelSE_an_%s.png - results for 10, 20, and 40 allele-specific reads per individuals (in files with suffix 100, 200 and 400 respectively)
3. Hessian behaviour when split and individual into multiple SNPs:
- We calculate hessian for allele-specific model under null (pi=0.5) with two over-dispersion parameters 0.1 and 0.5 (in this setup offdiagonal coefficients are 0).
Assuming that in one individual reads are split evenly into multiple SNPs we calculate hessian starting with 1 read per SNP in 8 SNP scenario (which corresponds to 2 reads in 4 SNPS, 4 reads in 2 SNPs and 8 reads per individual) and going up to 20 reads per SNP in 8 SNP scenario.
- hess_relation1_0.1.png & hess_relation1_0.5.png - diagonal hessian components 
- hess_logscale1_0.1.png & hess_logscale1_0.5.png - diagonal hessian components on log scale
- sd_by_asc_1ind_0.1.png & sd_by_asc_1ind_0.5.pig - implied standard deviation
- relSD_0.1.png & relSD_0.5.png - SD(#SNP)/SD(1SNP)
