# pipeline for eQTL analysis using 1000 Genome Project (1KGP) data. 
provide data snippets for two genes that can go throught the analysis

# To-dos

0. Add up the total amount of time for asSeq1/asSeq2, RASQUAL, and TReCASE score test for 5000 permutations, and our new appraoch with 1000 permuatations for 20 replicates, and add them in Supplementary Materials. 

+ 0.1 add the bamtools to asSeq2.

+ 0.2 check why TReC return NAs.

1. Setup the whole workflow

+ 1.1 Run DEseq2 to replace outliers using 4/n as cutoff. Need to include all the covariates (including sex) other than SNP genotype in this model. Make sure SNPs has been filtered by MAF 0.05. 

+ 1.2 Run Matrix eQTL (eigen MT) for all gene, SNP pairs, and choose to run TReCASE only for those pairs with Matrix eQTL p-value < 0.05

+ 1.3 Estimate the number of indepdent tests per gene 

+ 1.4 Obtain permutaion p-value of the minimum TReCASE p-value per gene

+ 1.5 Choose a permutation p-value cutoff across genes by controlling FDR. 

+ 1.6 Run the above steps 1.2-1.4 using permuted genotype data and check the distribution of permutation p-value. Compared the resutls when runing steps 1.3-1.4. 


2. Estimate the number of independent tests per gene

+ 2.1 using normal quantile transformed gene expression data to fit the MatrixEQTL and eigenMT. 

+ 2.1 Based on the results of eigenMT, estimate the original min p-value for permutation p-value in the range of 0.002 to 0.2. 

+ 2.2. Use the resisudals for each gene to bootstrap gene expression data after plugging in eQTL effects. 

+ 2.3. Run 1000 permutatinos to esitmate permuation p-value. 


3. Summarize the location of eQTLs. Find the strongest eQTL per gene, apply a threshold, and plot their locations by a density plot. 

+ 3.1 Make this plot for all three methods, MatrixEQTL, RASQUAL, and asSeq2, if overlapping of three density curves are not clear, you can just draw another one for two methods MatrixEQTL and asSeq2. May split the positve and negative strand. 

+3.2 Make this plot for the eQTLs identified only by Matrix EQTL, only by asSeq2 and by both. two eQTLs are the same if their R2 is > 0.8. 


4. Summarize the resutls by comparing annotations. 

Compare the reuslts of matrixEQTL vs. asSeq2, how frequent the resutls of one method overlap with some annotated regions, such as promoter, enhancer

+ 4.0 prepare a histogram of permutation p-values for MatrixEQTL and TReCASE, and a scatter plot of permutation p-values of the two methods. Check the p-value cutoffs for differnt FDR levels, such as 5% or 10%. 

+ 4.1 https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state
STATES FOR EACH 200bp BIN: 
https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/STATEBYLINE/

This should give you the state of each genomic location across many tissues and cell lines. find the cell line that match with our data, lymphoblast, if there is no such cell line, we should discuss which one is most similar. 


+4.2. https://egg2.wustl.edu/roadmap/web_portal/DNase_reg.html#delieation



6. Check whether GWAS signals are enriched in eQTL sites, for example, using LD regression. 
you can use the information from this page

https://github.com/bulik/ldsc/wiki/Partitioned-Heritability

https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses

LD regressoin use GWAS test statistic as response variable, covariates are the LD score, which calculate the summation of R2 of this SNP versus all nearby SNPs. Some example codes are here

#!/bin/bash
module load anaconda2
source activate ldsc
ml bedtools
ct="selected.Astro"
for chr in {1..22}
do
python make_annot.py \
		--gene-set-file ../data/${ct}.txt \
		--gene-coord-file ../data/geneList_all_hg19_final.txt \
		--windowsize 100000 \
		--bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
		--annot-file ../data/${ct}.${chr}.annot.gz
done

for chr in {1..22}
do
python ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
    --ld-wind-cm 1 --annot ../data/${ct}.${chr}.annot.gz --thin-annot \
    --out ../data/${ct}.${chr} --print-snps hapmap3_snps/hm.${chr}.snp
done
source deactivate


python munge_sumstats.py \
--sumstats ../data/clozuk_pgc2.meta.sumstats.txt.gz \
--out ../data/clozuk_pgc2.meta.sumstats \
--N 105318 --merge-alleles \
../data/clozuk_pgc2.meta.sumstats.info9.snplist.txt.gz



python ldsc.py \
--h2-cts ../data/clozuk_pgc2.meta.sumstats_with_rsid.gz \
--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
--out ../data/clozuk_pgc2.meta_ct6 \
--ref-ld-chr-cts ../data/ct6.ldcts \
--w-ld-chr weights_hm3_no_hla/weights.




(in progress)
