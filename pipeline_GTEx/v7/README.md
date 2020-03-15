# pipeline for eQTL analysis using GTEx data.

Genotype calls and RNA-seq data for whole blood sample are downloaded from dbGaP at https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v7.p2. 

### 1. data processing
The data preparation for genotype data is similar to that in pipeline_TCGA/doc/workflow_asSeq.Rmd except that:
  
  + we modified step1 in stepA - Indels, multiallelic SNPs and duplicates are removed. 
    - genotype_data_prep/step1_prepare_pedfile.sh
  
  + we skiped step 6 (imputation) in stepB. 

### 2. Run TReCASE

   R script: TReCASE/
   
   + step1_trecase.R -- TReCASE analysis with replacing outlier steps. 
   
### 3. Run RASQUAL: 
similar to pipeline_1KGP.

   R script: RASQUAL/
   
### 4. LD Score: 

   + step1_create_eQTL_geneList.R -- create control gene list, eQTL gene list, and gene coordinate file.
      
      - coordinate file: all the genes used in eQTL analysis, with its chromosome infomation, and start and end location to be the location of the most significant eQTL SNP
      
      - control gene list: genes whose most significant eQTL SNP has q-value > 0.05
      - eQTL gene list: genes whose most significant eQTL SNP has q-value <= 0.05
   
   + anno_ldsc.sh and control.sh -- annotating and computing LD score for eQTL gene list and control gene list respectively. 
   
   + ldscore_analysis.sh -- partitioning LD score
   

