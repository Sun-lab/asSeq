chr1   = 2
eStart = 925
eEnd   = 926

library(asSeq) 


# ---------------------------------------------------------------------
# test for cis-eQTL
# ---------------------------------------------------------------------

setwd("~/research/eQTL_seq/result/YRI_Joint_Permute/")

output1 = sprintf("TReCASE_Permute_YRI_expression_chr%s_info.txt", chr1)
einfo   = read.table(output1, sep="\t", header=TRUE)
dim(einfo)
einfo[1:2,1:7]

X = read.table("Covariates.txt", sep = "\t", header = TRUE)
dim(X)
X[1:2,1:5]

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


Y  = data.matrix(Y[,eStart:eEnd])
Y1 = data.matrix(Y1[,eStart:eEnd])
Y2 = data.matrix(Y2[,eStart:eEnd])
eChr = eChr[eStart:eEnd]
ePos = ePos[eStart:eEnd]

ta = trecaseP(Y, Y1, Y2, X, Z, min.nTotal=5, min.N=5, min.Nhet=5, 
              cis.only = TRUE, cis.distance = 2e+05, eChr = eChr, 
              ePos = ePos, mChr = mChr, mPos = mPos, trace=6, 
              np.max=10, np=c(5), aim.p=c(0.4))

