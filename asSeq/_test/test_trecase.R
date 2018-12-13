
library(MASS)
source("~/research/R/asSeq/_test/aseR.R")
source("~/research/R/asSeq/_test/trecaseR.R")

library(asSeq)

# ---------------------------------------------------------------------
# test for cis-eQTL
# ---------------------------------------------------------------------

chr1 = 2

# ---------------------------------------------------------------------
# test one case where glmNBlog fails but glmNB works
#
# chr1 = 2
# ---------------------------------------------------------------------

setwd("~/research/eQTL_seq/result/YRI_Joint_Permute/")

output1 = sprintf("TReCASE_Permute_YRI_expression_chr%s_info.txt", chr1)
einfo   = read.table(output1, sep="\t", header=TRUE)
dim(einfo)
einfo[1:2,1:7]

X = read.table("Covariates.txt", sep = "\t", header = TRUE)
dim(X)
X[1:2,]

X = data.matrix(X)

cn = read.table("YRI_samples.txt", sep="\t", header=TRUE)
dim(cn)

chr2 = sprintf("chr%d", chr1)

st1 = "\n---------------------------------------------------\n"
message(sprintf("%s chr: %d %s", st1, chr1, st1))

output1 = sprintf("TReCASE_Permute_YRI_geno_chr%s_info.txt", chr1)
output2 = sprintf("TReCASE_Permute_YRI_geno_chr%s_data.txt", chr1)

genoInf = read.table(output1, sep="\t", header=TRUE)
genoDat = read.table(output2, sep="\t", header=TRUE)
dim(genoInf)
dim(genoDat)
genoInf[1:2,]
genoDat[1:2,1:5]

if(any(cn$name != names(genoDat))){
  stop("sample ID do not match\n")
}

X = cbind(X, as.numeric(cn$gender))

output2 = sprintf("TReCASE_Permute_YRI_expression_chr%s_data_TReC.txt", chr1)
output3 = sprintf("TReCASE_Permute_YRI_expression_chr%s_data_ASE1.txt", chr1)
output4 = sprintf("TReCASE_Permute_YRI_expression_chr%s_data_ASE2.txt", chr1)

Y       = read.table(output2, sep="\t")
Y1      = read.table(output3, sep="\t")
Y2      = read.table(output4, sep="\t")
Z       = t(data.matrix(genoDat))

dim(Y1)
dim(Y2)
dim(Y)
dim(Z)

Y1[1:2,1:8]
Y2[1:2,1:8]
Y[1:2,1:8]
Z[1:2,1:8]

chrI = as.integer(chr1)
eChr = rep(chrI, ncol(Y1))
ePos = einfo$txStart

mChr = rep(chrI, ncol(Z))
mPos = genoInf$position_b36

# ---------------------------------------------------------------------
# test one case where glmNBlog fails but glmNB works
#
# eID = 46
# mID = 39546
#---------------------------------------------------------------------

eID = 925
mID = 18082

y   = Y[,eID]
y1  = Y1[,eID]
y2  = Y2[,eID]
z1  = Z[,mID]

setwd("~/")
trecase(y, y1, y2, X, z1, output.tag="test_temp", p.cut=0.1, eChr=chr1,
   ePos=ePos[eID], mChr=mChr[mID], mPos=mPos[mID], trace=6, maxit=100)

rest = read.table("test_temp_eqtl.txt", header=TRUE, sep="\t")
rest

system("rm test_temp_eqtl.txt")
system("rm test_temp_freq.txt")

# ---------------------------------------------------------------------
# test TreCASE model
# ---------------------------------------------------------------------

eID = which(einfo$name2=="ENSG00000115934")
mID = which(genoInf$rsID=="rs2136600")

y   = Y[,eID]
y1  = Y1[,eID]
y2  = Y2[,eID]
z1  = Z[,mID]
z2  = z1
z1[z1 > 2] = z1[z1 > 2] - 2

table(z1)
table(z2)
ta = trecaseR(y, y1, y2, X, z1, z2, plotIt = FALSE)

ta

plot(ta$logLik)

trecase(y, y1, y2, X, Z[,mID], output.tag="test_temp", p.cut=0.1, eChr=chr1,
ePos=ePos[eID], mChr=mChr[mID], mPos=mPos[mID], trace=1, maxit=100)

rest = read.table("test_temp_eqtl.txt", header=TRUE, sep="\t")
rest

system("rm test_temp_eqtl.txt")
system("rm test_temp_freq.txt")
