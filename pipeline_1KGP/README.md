# pipeline for eQTL analysis using 1000 Genome Project (1KGP) data. 
For real data analysis I provided adding real data we provide data snippets for two genes that can go throught the analysis

1. Using mapped files and data from TRECASE_MLE process according to RASQUAL approach using GATK

2. Recreate VCF file with these GATK counts

3. Fit RASQUAL

permuted data analysis would have one correction to the second step: - permuting SNP information in VCF file
