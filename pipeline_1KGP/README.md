# pipeline for eQTL analysis using 1000 Genome Project (1KGP) data. 
For real data analysis I provided adding real data we provide data snippets for two genes that can go throught the analysis

1. Using mapped files and data from TRECASE_MLE process according to RASQUAL approach using GATK

2. Recreate VCF file with these GATK counts

3. Fit RASQUAL

# To-dos



+ 0.a Check how much time does it take to run TReC or TReCASE for all genes and all SNPs. 

+ 0.b Re-run eQTL mapping using asSeq2 and compare the results between asSeq1 and asSeq2. 

+ 0.c include the code of eigenMT in the repository

+ 0.d Get 2 x 2 table for comparing the permtution p-value and estiamated permutation pvalue

+ 0.e Generate p-value histogram using TReC/TReCASE, RASQUAL after permuting SNP genotype data. 

1. Compare the results from MatrixEQTL, TReC, TReCASE and RASQUAL. For example, by two-way table for different p-value cutoffs, say 0.1, 0.01, 1e-3, 1e-6, 1e-8. Consider p-value per (gene, SNP) pair instead using minimum p-value. something like this. 

```
              (0,1e-06] (1e-06,0.001] (0.001,0.01] (0.01,0.1] (0.1,1]   total
(0,1e-06]         91833         10712          743        226     463  103977
(1e-06,0.001]     20825        131456        32326       6463    2302  193372
(0.001,0.01]       3044         57216       137675      75936   10333  284204
(0.01,0.1]         4954         32365       132427     645987  295067 1110800
(0.1,1]           15012         43271        69960     515532 5830884 6474659
total            135668        275020       373131    1244144 6139049 8167012
```

+ 1.0 Run MatrixEQTL after normal quantile transformation of gene expression gene by gene. 

+ 1.1 For those with MatrixEQTL p-value > 0.1 and TReC p-value < 1e-6, check whehter there are outliers by plotting some examples. 

+ 1.2 TReC p-value < 1e-6 and Matrix eQTL p-value > 0.1, check their TReC p-value after removing top 3 samples with largest read counts. 

+ 1.3 Generate the scatter plot of -log10(min MatrixEQTL p-value) vs. -log10(min TReC p-value), the plot may be messy, then the two-way table is better

+ 1.4 generate a histogram of -log10(min MatrixEQTL p-value)\

+ 1.5 Check for those cases with very small MatrixEQTL p-values but large TReC p-values. 


2. For each gene, we got multiple p-values across multiple SNPs, how to get a permutation p-value or an approximation. Let p_min be the minimum p-value of this gene across all SNPs, and let d_eff be effect number of independent tests. Then our p-value for this gene is p_gene = min(p_min * d_eff, 1). To define the number of independent tests, we want fit a linear model of the relation between minimum p-value versus permution p-value. We can perturb the expression of this gene a few times, for example, by parametric bootstrap given SNP genotype and the linear regression betwween gene expression and SNP genotype, and get more data points. Maybe we can use 10 data points and permute 1000 times for each data points. Thus this requires 10,000 permutation per gene in total. Here when you do parametric bootstrap, you want to choose the eQTL effects in a range so that the permutation p-value spread the range of 0.001 and 0.1. 

+ 2.0 For our estimate, if the estimated number of tests is larger than the number of SNPs, we will truncated it to be the number of SNPs. 

+ 2.1 We can do this for 100 genes first and compare the d_eff with some simple methods such as eigenMT: https://www.ncbi.nlm.nih.gov/pubmed/26749306 

+ 2.2 Summarize the number of SNPs per gene and the number of independent SNPs per gene and generate a scatter plot. 

+ 2.3 Use an empirical Bayes approach to re-estimae intercept and slope. Using the results across all genes to estimte distribution of intercept and slopes by a normal distribution. Use those distribution as prior to re-estiamte intercept and slope. Compare the MLE vs. posterior by scatter plot. 

3. Summarize the location of eQTLs. Find the strongest eQTL per gene, apply a threshold, and plot their locations by a density plot. 

4. Check whether GWAS signals are enriched in eQTL sites, for example, using LD regression. 
