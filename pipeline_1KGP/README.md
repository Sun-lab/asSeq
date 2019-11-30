# pipeline for eQTL analysis using 1000 Genome Project (1KGP) data. 
provide data snippets for two genes that can go through the analysis

# To-dos

+ 0.1 add the bamtools to asSeq2.

+ 0.2 check why TReC return NAs.

1. Setup the whole workflow

+ 1.1 Run DEseq2 to replace outliers using 4/n as cutoff. Need to include all the covariates (including sex) other than SNP genotype in this model. Make sure SNPs has been filtered by MAF 0.05. 

+ 1.2 Run Matrix eQTL (eigen MT) for all gene, SNP pairs, and choose to run TReCASE only for those pairs with Matrix eQTL p-value < 0.05

+ 1.3 For analysis with a subset of samples, check the SNP has been filtered by MAF 0.05. 

2. Summarize the location of eQTLs. Find the strongest eQTL per gene, apply a threshold, and plot their locations by a density plot. 

+ 2.1 Make this plot for the eQTLs identified only by Matrix EQTL, only by asSeq2 and by both. two eQTLs are the same if their R2 is > 0.8. 

3. Summarize the results from GTEx whole blood. 

+ 3.1 make sure the gene location file is ok. 

+ 3.2 filter SNPs using MAF 0.05 for the permutations. 

+ 3.3 summarize the results using permutation p-value based FDR cutoff, probably using FDR 0.01 as cutoff. For example, compare the # of findings with respect to the number of fSNPs. 

4. Summarize the results by comparing annotations. 

Compare the results of matrixEQTL vs. asSeq2, how frequent the results of one method overlap with some annotated regions, such as promoter, enhancer

+ 4.1 https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state
STATES FOR EACH 200bp BIN: 
https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/STATEBYLINE/. This should give you the state of each genomic location across many tissues and cell lines. find the cell line that match with our data, lymphoblast, if there is no such cell line, we should discuss which one is most similar. 

+ 4.2. https://egg2.wustl.edu/roadmap/web_portal/DNase_reg.html#delieation

