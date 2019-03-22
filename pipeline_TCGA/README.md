# asSeq
allele-specific sequencing data analysis. 

the pipeline to obtain allele-specific sequencing reads. Including genotype imputation and phasing, extration of allele-specific reads and counting the allele-speciic reads. The files are organized as follows. 
```
├── code
│   ├── flow.txt
│   ├── stepA_prepare_input_data
│   ├── stepB_phasing_and_imputation
│   └── stepC_collect_read_count
├── data
│   ├── chr
│   ├── chr_flipped
│   ├── data_genotype
│   ├── data_snp
│   ├── imputed
│   ├── phased
│   └── sample_geno_call_22.txt
└── doc
    ├── workflow_asseq.Rmd
    └── workflow_asseq.html
```

In the folder ```doc```, you can see an R markdown file and the ouptut in html format, which provide the explanation and example coede for the whole pipeline. 

In the folder ```code```, R or shell scripts are divided into three steps (step A, B, and C) and are saved in three folders. 

The folder ```data``` contains example data. 


