<div align="left">
<a href=""><img src="https://img.shields.io/badge/R-%23276DC3.svg?style=square&logo=r&logoColor=pink&label=asSeq" height="100" /></a>
</div>

<!-- badges: start -->
![R](https://img.shields.io/badge/R-%23276DC3.svg?style=square&logo=r&logoColor=pink)
![CRAN status](https://www.r-pkg.org/badges/version/asSeq)
[![DOI](https://zenodo.org/badge/DOI/10.1111/j.1541-0420.2011.01654.x.svg)](https://doi.org/10.1111/j.1541-0420.2011.01654.x)
[![DOI](https://zenodo.org/badge/DOI/10.1371/journal.pgen.1010076.svg)](https://doi.org/10.1371/journal.pgen.1010076)
<!-- badges: end -->

## Introduction

eQTL mapping using allele-specific gene expression data

This repository includes the package asSeq that implements the method **TReCASE** (Sun 2012). 
Please refer to another repository for all the pipelines for data anlaysis: 
https://github.com/Sun-lab/asSeq_pipelines

## Installation

To install this package in R, use 
 
```R
library("devtools");
install_github("Sun-lab/asSeq")
```

There are many warning messages, and most of them are due to the codes we copied from bamtools. None of the warning message should affect the accuracy of our methods or the success of installation. 

The success of compiling of the c/c++ codes and installation may depend on the compiler and handling of temporary files for each machine/operating system. For example, in a Mac environment, the package installation in Rstudio may fail because of permission to access system folders. In such situations, download the source package (e.g., [asSeq_0.99.501.tar.gz](https://github.com/Sun-lab/asSeq/raw/master/asSeq_0.99.501.tar.gz)) and run the following command in terminal should work

```Shell
R CMD install asSeq_0.99.501.tar.gz
```

## References

[Sun, W.](https://github.com/sunway1999) (2012). A statistical framework for eQTL mapping using RNA-seq data. *Biometrics*, 68(1), 1-11. [[HTML](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3218220/), [PDF](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3218220/pdf/nihms-307768.pdf)]

[Zhabotynsky, V.](https://github.com/yaceya), [Huang, L.](https://github.com/licaih), [Little, P.](https://github.com/pllittle), [Hu, Y. J.](https://sph.emory.edu/faculty/profile/index.php?FID=yijuan-hu-8694), [Pardo-Manuel de Villena, F.](https://www.med.unc.edu/genetics/directory/fernando-pardo-manuel-de-villena-phd/), [Zou, F.](https://sph.unc.edu/adv_profile/fei-zou-phd/), [Sun, W.](https://github.com/sunway1999) (2022). eQTL mapping using allele-specific count data is computationally feasible, powerful, and provides individual-specific estimates of genetic effects. *PLoS Genetics*, 18(3), e1010076. [[HTML](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010076)]
