
# ------------------------------------------------------------
# datM is a read count matrix, each row corresponds
# to a gene and each column corresponds to a sample
#
# tot is the total number of reads per sample
# s75 is the 75 percentile of gene expression
# ------------------------------------------------------------

tot = colSums(datM)
s75 = apply(datEA, 2, quantile, prob=0.75)

# ------------------------------------------------------------
# check the correlation between tot and s75
# ------------------------------------------------------------

cor(tot, s75)

# ------------------------------------------------------------
# normalize the gene expression data by 75 percentile
# ------------------------------------------------------------

nDat = t(log10(t((datEA + 1))/s75))
dim(nDat)

# ------------------------------------------------------------
# remove the mean value of each gene, and then set missing 
# value as 0. Certainly a better missing vlaue imputation 
# method is desirable if there are larger number of missings
# ------------------------------------------------------------

datR14Pr = nDat - rowMeans(nDat, na.rm=TRUE)
datR14Pr[is.na(datR14Pr)] = 0 

# ------------------------------------------------------------
# eigen value decomposion of the matrix with the number of 
# rows and the number of samples equal to sample size
# ------------------------------------------------------------

covdatR1 = t(datR14Pr) %*% datR14Pr / nrow(datR14Pr)
dim(covdatR1)
prdatR1  = eigen(covdatR1)

# ------------------------------------------------------------
# check the 20 eigen-values, and first 3 PCs
# ------------------------------------------------------------

prdatR1$values[1:20]

PC1 =  prdatR1$vectors[,1]
PC2 =  prdatR1$vectors[,2]
PC3 =  prdatR1$vectors[,3]

