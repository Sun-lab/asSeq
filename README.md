# asSeq
eQTL mapping using allele-specific gene expression data


This repository includes the package asSeq that implements the method TReCASE (Sun 2012). 
Please refer to another repository for all the pipelines for data anlaysis: 
https://github.com/Sun-lab/asSeq_pipelines


## Installation 
 To install this package in R, use 
 
 ```
    library("devtools");
    install_github("Sun-lab/asSeq")
 ```

There are many warning messages, and most of them are due to the codes we copied from bamtools. None of the warning message should affect the accuracy of our methods or the success of installation. 

The success of compiling of the c/c++ codes and installation may depend on the compiler and handling of temporary files for each machine/operating system. For example, in a Mac environment, the package installation in Rstudio may fail because of permission to access system folders. In such situations, download the source package (e.g., asSeq_0.99.501.tar.gz) and run the following command in terminal should work

 ```
R CMD install asSeq_0.99.501.tar.gz
 ```



## Reference

Sun, W. (2012). A statistical framework for eQTL mapping using RNA‐seq data. Biometrics, 68(1), 1-11.

Zhabotynsky et al. (2021) eQTL mapping using allele-specific geneexpression
