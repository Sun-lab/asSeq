# pipeline for eQTL analysis using 1000 Genome Project (1KGP) data. 
provide data snippets for two genes that can go throught the analysis

# To-dos

+ 0.1 add the bamtools to asSeq2.

+ 0.2 check why TReC return NAs.

1. Setup the whole workflow

+ 1.1 Run DEseq2 to replace outliers using 4/n as cutoff. Need to include all the covariates (including sex) other than SNP genotype in this model. Make sure SNPs has been filtered by MAF 0.05. 

+ 1.2 Run Matrix eQTL (eigen MT) for all gene, SNP pairs, and choose to run TReCASE only for those pairs with Matrix eQTL p-value < 0.05

+ 1.3 For anlaysis with a subset of samples, check the SNP has been filtered by MAF 0.05. 

2. Summarize the location of eQTLs. Find the strongest eQTL per gene, apply a threshold, and plot their locations by a density plot. 

+2.1 Make this plot for the eQTLs identified only by Matrix EQTL, only by asSeq2 and by both. two eQTLs are the same if their R2 is > 0.8. 


4. Summarize the resutls by comparing annotations. 

Compare the reuslts of matrixEQTL vs. asSeq2, how frequent the resutls of one method overlap with some annotated regions, such as promoter, enhancer

+ 4.1 https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state
STATES FOR EACH 200bp BIN: 
https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/STATEBYLINE/. This should give you the state of each genomic location across many tissues and cell lines. find the cell line that match with our data, lymphoblast, if there is no such cell line, we should discuss which one is most similar. 

+ 4.2. https://egg2.wustl.edu/roadmap/web_portal/DNase_reg.html#delieation

