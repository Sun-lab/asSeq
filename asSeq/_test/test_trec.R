
# -------------------------------------------------------------------------
# read in data
# -------------------------------------------------------------------------
setwd("~/research/eQTL_seq/result/YRI_Joint_Permute/")

eD = read.table("TReCASE_Permute_YRI_expression_chr2_data_TReC.txt", 
  sep="\t")

X = read.table("Covariates.txt", sep = "\t", header = TRUE)
dim(X)
X[1:2,1:5]
X = data.matrix(X)


geno = read.table("TReCASE_Permute_YRI_geno_chr2_data.txt", sep="\t", header=TRUE)
dim(geno)

genoInfo = read.table("TReCASE_Permute_YRI_geno_chr2_info.txt", sep="\t", header=TRUE)
dim(genoInfo)
genoInfo[1:5,]

## eQTL pairs identified by TReCASE
eID = 46
mID = 39568

# -------------------------------------------------------------------------
# check Trec model
# -------------------------------------------------------------------------

y  = as.numeric(eD[,eID])
z1 = as.numeric(geno[mID,])
z2 = z1
z2[z1==3] = 1
z2[z1==4] = 2

library(asSeq)
trec(y, X, z2, output.tag="test_temp", p.cut=1.0, cis.only = FALSE)

rest = read.table("test_temp_eqtl.txt", header=TRUE, sep="\t")
rest

system("rm test_temp_eqtl.txt")
system("rm test_temp_freq.txt")

# -------------------------------------------------------------------------
# check linear model
# -------------------------------------------------------------------------

normscore = function(vec) {
  len  = length(na.omit(vec))+1
  rank = rank(na.omit(vec))
  ties = (rank - floor(rank)) > 0
  new.vec = vec[!is.na(vec)] 
  new.vec[!ties]=qnorm(rank[!ties]/len)
  new.vec[ties] =0.5*(qnorm((rank[ties]+0.5)/len)+qnorm((rank[ties]-0.5)/len))
  vec[!is.na(vec)] = new.vec
  vec
}

z3 = as.numeric(geno[which(genoInfo$rsID == "rs4359651"),])
z4 = z3
z4[z3==3] = 1
z4[z3==4] = 2

y2 = normscore(y)
l3 = lm(y2 ~ X + z4)
summary(l3)

# -------------------------------------------------------------------------
# check TReC model using the marker identified by lienar model
# -------------------------------------------------------------------------

library(asSeq)

trec(y, X, z4, output.tag="test_temp", p.cut=1.0, cis.only = FALSE)

rest = read.table("test_temp_eqtl.txt", header=TRUE, sep="\t")
rest

system("rm test_temp_eqtl.txt")
system("rm test_temp_freq.txt")

# -------------------------------------------------------------------------
# check TReC model using the marker identified by lienar model
# without adjusting Z, i.e., use glmNB instead of glmNBlog
# -------------------------------------------------------------------------

trec(y, X, z4, output.tag="test_temp", p.cut=1.0, cis.only = FALSE, adjZ=FALSE)

rest = read.table("test_temp_eqtl.txt", header=TRUE, sep="\t")
rest

system("rm test_temp_eqtl.txt")
system("rm test_temp_freq.txt")

# -------------------------------------------------------------------------
# check results of glm.nb
# -------------------------------------------------------------------------

library(MASS)
g1 = glm.nb(y ~ X + z4)
g0 = glm.nb(y ~ X)
anova(g0, g1)
