args=(commandArgs(TRUE))
#args = c("specifications_Brain_Caudate_basal_ganglia.txt", "5e5")
args

rt.dir = getwd()
specf = args[1]
cis_window = as.numeric(args[2])

specs = unlist(read.table(specf, as.is=T))
pref = specs[1]
wrk.dir = sprintf("%s/%s", rt.dir, pref)
setwd(wrk.dir)

specs = unlist(read.table("specifications.txt", as.is=T))
specs
pref = specs[1]
nsam = specs[2]
queue = specs[3]
days = specs[4]
#wrk.dir = sprintf("%s/%s", rt.dir, pref)
bas.dir = specs[16]
wrk.dir = specs[14]
lib.dir = specs[15]
mem = sprintf("%sg", specs[5])
mem

source(sprintf("%s/helpers.R", lib.dir))

library(VGAM)

#setwd(sprintf("/pine/scr/z/h/zhabotyn/R01/GTEx/v8/2020_01_20/%s", pref))

#cov.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8/Annotations/GTEx_Analysis_v8_eQTL_covariates"
#covs = read.table(sprintf("%s/Brain_Caudate_basal_ganglia.v8.covariates.txt", cov.dir), as.is=T)
#pheno.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8/Annotations"
#pheno = read.table(sprintf("%s/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", pheno.dir), as.is=T, sep="\t", header=T)
anno.dir = sprintf("%s/Annotations", bas.dir)
covs.dir = sprintf("%s/GTEx_Analysis_v8_eQTL_covariates", anno.dir)
pheno.fil = sprintf("%s/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", anno.dir)
pheno = read.table(pheno.fil, as.is=T, sep="\t", header=T)
pheno[1:4,]

pheno$agen = 25
pheno$agen[pheno$AGE=="20-29"] = 25
pheno$agen[pheno$AGE=="30-39"] = 35
pheno$agen[pheno$AGE=="40-49"] = 45
pheno$agen[pheno$AGE=="50-59"] = 55
pheno$agen[pheno$AGE=="60-69"] = 65
pheno$agen[pheno$AGE=="70-79"] = 75
pheno$sex = -1
pheno$sex[pheno$SEX==2] = 1
#covars = read.table(sprintf("/pine/scr/z/h/zhabotyn/R01/GTEx/v8/Annotations/GTEx_Analysis_v8_eQTL_covariates/%s.v8.covariates.txt",pref), header=T, as.is=T)
covars = read.table(sprintf("%s/%s.v8.covariates.txt",covs.dir, pref), header=T, as.is=T)
rownames(covars) = covars$ID; covars = covars[,-1]
colnames(covars) = gsub("\\.", "-", colnames(covars))
covars[1:2,1:5]

modelL = "long"
modelS = "short"
resL = read.csv(sprintf("%s_%s_5e+05_long_TReCASE_trimmed.csv", pref, nsam), as.is=T)
resS = read.csv(sprintf("%s_%s_5e+05_short_TReCASE_trimmed.csv", pref, nsam), as.is=T)
resL[1:4,]
resS[1:4,]

cntL = read.csv(sprintf("counts_with_min_snp_%s_%s_5e+05_long.csv", pref, nsam), as.is=T)
cntS = read.csv(sprintf("counts_with_min_snp_%s_%s_5e+05_short.csv", pref, nsam), as.is=T)
colnames(cntL) = gsub("\\.", "-", colnames(cntL))
colnames(cntS) = gsub("\\.", "-", colnames(cntS))
rownames(cntL) = cntL[,1]; cntL = cntL[,-1]; cntL = as.matrix(cntL)
rownames(cntS) = cntS[,1]; cntS = cntS[,-1]; cntS = as.matrix(cntS)
cntL[1:4,1:4]
cntS[1:4,1:4]

phenoL = pheno[match(colnames(cntL),pheno$SUBJID),]
m = match(phenoL$SUBJID, colnames(covars));m
covars = covars[,m]
table(colnames(covars) == phenoL$SUBJID)

eqLi = sprintf("TReCASE_%s_%s_5e+05_long.csv", pref, nsam, as.is=T)
eqSi = sprintf("TReCASE_%s_%s_5e+05_short.csv", pref, nsam, as.is=T)
eqLi = read.csv(eqLi)
eqSi = read.csv(eqSi)

cand = c(grep("ENSG00000141510", rownames(cntL)), grep("ENSG00000102974", rownames(cntL)))
totL = apply(cntL, 1:2, splitting, block=1)
crv = colSums(totL); crv = median(crv)/crv
nrmL = totL%*%diag(crv);summary(colSums(nrmL))
phenoL$tp53 = log10(nrmL[cand[1],])
phenoL$ctcf = log10(nrmL[cand[2],])

alpha1 = 0.01
kp = which(eqLi$permp<alpha1 & eqLi$trans_Pvalue>=.1); length(kp)
cntLi = cntL[kp,]
kp = which(eqSi$permp<alpha1 & eqSi$trans_Pvalue>=.1); length(kp)
cntSi = cntS[kp,]


minsL = apply(cntLi, 1:2, splitting)
minsS = apply(cntSi, 1:2, splitting)



hap1L = apply(cntLi, 1:2, splitting, block=2)
hap2L = apply(cntLi, 1:2, splitting, block=3)
hap1S = apply(cntSi, 1:2, splitting, block=2)
hap2S = apply(cntSi, 1:2, splitting, block=3)


calc_prop = function(i, hap1, hap2, msnp, cutoff=10, off=.5){
  aA = hap1[,i]
  aB = hap2[,i]  
  table(is.na(aA))
  aA[msnp[,i] %in% c(0,4)] = NA
  aB[msnp[,i] %in% c(0,4)] = NA
  table(is.na(aA))
  flip = which(msnp[,i]==3)
  tmp = aA[flip]  
  aA[flip] = aB[flip]
  aB[flip] = tmp 
  aC = aA+aB<cutoff
  aA[aC] = NA
  aB[aC] = NA
  table(is.na(aA))
  eff = log(aB+off)-log(aA+off)
  eff
}
fracL = sapply(1:ncol(hap1L), calc_prop, hap1=hap1L, hap2=hap2L, msnp=minsL)
fracS = sapply(1:ncol(hap1S), calc_prop, hap1=hap1S, hap2=hap2S, msnp=minsS)
rownames(fracL) = rownames(hap1L)
rownames(fracS) = rownames(hap1S)
colnames(fracL) = colnames(fracS) = colnames(hap1L)

dim(fracL)
dim(fracS)
cutoff2 = 20
fracL = fracL[rowSums(!is.na(fracL))>cutoff2,]
fracS = fracS[rowSums(!is.na(fracS))>cutoff2,]
dim(fracL)
dim(fracS)

fracL = fracL - rowMeans(fracL, na.rm=TRUE)
fracS = fracS - rowMeans(fracS, na.rm=TRUE)


cutoff3 = .1
sdL = apply(fracL, 1, sd, na.rm=T)
sdS = apply(fracS, 1, sd, na.rm=T)
summary(sdL)
summary(sdS)
fracL = fracL[sdL>cutoff3,]
sdL = sdL[sdL>cutoff3]
fracS = fracS[sdS>cutoff3,]
sdS = sdS[sdS>cutoff3]
summary(sdL)
summary(sdS)
for(i in 1:nrow(fracL))fracL[i,] = fracL[i,]/sdL[i]
for(i in 1:nrow(fracS))fracS[i,] = fracS[i,]/sdS[i]
sdL0 = apply(fracL, 1, sd, na.rm=T)
sdS0 = apply(fracS, 1, sd, na.rm=T)
summary(sdL0)
summary(sdS0)


numL = effL = numS = effS = matrix(0, nrow=ncol(fracL), ncol=ncol(fracL))
for(coli in 1:ncol(fracL)){
  for(rowi in 1:ncol(fracL)){
    x = fracL[,rowi]
    y = fracL[,coli]
    kpij = !is.na(x)&!is.na(y)
    #cat(sum(kpij), " ")
    numL[rowi,coli] = sum(kpij)
    effL[rowi,coli] = cor(x[kpij],y[kpij])#cov(x[kpij],y[kpij])

    x = fracS[,rowi]
    y = fracS[,coli]
    kpij = !is.na(x)&!is.na(y)
    #cat(sum(kpij), " ")
    numS[rowi,coli] = sum(kpij)
    effS[rowi,coli] = cor(x[kpij],y[kpij])#cov(x[kpij],y[kpij])
  }
  #cat("\n")
  if(coli%%25==0)message(coli)
}
summary(c(numL))
summary(c(numS))
effL[1:4,1:4]

#prL0 = prcomp(effL)
prL  = eigen(effL)
prL$values[1:20]

prS  = eigen(effS)
prS$values[1:20]

pdf(sprintf("%s_ASE_cor_PCs_short.pdf", pref), width=6, height=6)
#short

val = prS$values[1:20]
PC1 =  prS$vectors[,1]
PC2 =  prS$vectors[,2]
PC3 =  prS$vectors[,3]
PC4 =  prS$vectors[,4]

par(mar=c(5,4,1,1), mfrow=c(2,2))
barplot(val, main="", xlab="Index", ylab="Eigen-value, short")
cols = heat.colors(n=6)
ageg = sort(names(table(phenoL$age)))
legend("topright", ageg, bty="n", text.col=cols)
legend("top", c("male", "female"), pch=1:2, bty="n")
#cols = gray.colors(n=6, start = 0, end = 0.7, gamma = 2.2, alpha=1)
coli = cols[factor(phenoL$age, levels=ageg)]
pchi = phenoL$SEX
cexi = phenoL$tp53-min(phenoL$tp53)+.5
#cexi = cexi/max(cexi)
plot(PC1, PC2, cex=cexi, col=coli, pch=pchi, bty="n", main="")
plot(PC1, PC3, cex=cexi, col=coli, pch=pchi, bty="n", main="")
plot(PC2, PC3, cex=cexi, col=coli, pch=pchi, bty="n", main="")



pchi = phenoL$DTHHRDY+1
cexi = phenoL$ctcf-min(phenoL$ctcf)+.5
barplot(val, main="", xlab="Index", ylab="Eigen-value, short")
legend("topright", legend=sort(unique(phenoL$DTHHRDY)), pch=1:5, bty="n")
legend("top", c("male", "female"), pch=1:2, bty="n")
#cexi = cexi/max(cexi)
plot(PC1, PC2, cex=cexi, col=coli, pch=pchi, bty="n", main="")
plot(PC1, PC3, cex=cexi, col=coli, pch=pchi, bty="n", main="")
plot(PC2, PC3, cex=cexi, col=coli, pch=pchi, bty="n", main="")


par(mar=c(5,4,1,1), mfrow=c(2,2))
corc = matrix(NA, nrow=nrow(covars), ncol=20)
corg = matrix(NA, nrow=nrow(fracS), ncol=20)
rownames(corg) = rownames(fracS)
rownames(corc) = rownames(covars)

for(i in 1:nrow(corc)){
y = unlist(covars[i,])
for(j in 1:20){
  corc[i,j] = cor(prS$vectors[,j], y)
}
}

for(i in 1:20){
  x = 1:nrow(corc)
  plot(x, corc[,i], bty="n", xlab="covariates", ylab="correlation", main=sprintf("PC%s", i))
}

for(pci in 1:20){
  kp = which(abs(corc[,pci])>.25)
for(covi in kp){
  plot(x=unlist(covars[covi,]), y=prS$vectors[,pci], main=sprintf("PC%sx%s", pci, rownames(covars)[covi]), xlab=rownames(covars)[covi], bty="n")
}
}

for(i in 1:nrow(corg)){
x = fracS[i,]
for(j in 1:20){
  y = prS$vectors[,j]
  kp = !is.na(x)&!is.na(y)
  corg[i,j] = cor(x[kp], y[kp])
}
}
par(mar=c(5,4,1,1), mfrow=c(2,2))
for(i in 1:20){
  x = 1:nrow(corg)
  #plot(x, corg[,i], bty="n", xlab="covariates", ylab="correlation", main=sprintf("PC%s", i))
  plot(density(corg[,i]), bty="n", xlab="correlation", main=sprintf("PC%s", i))
}
par(mar=c(5,4,1,1), mfrow=c(2,2))
for(i in 1:20){
  kp = which(abs(corg[,i])>.5)
  for(j in kp){
    plot(fracS[j,], prS$vectors[,i], xlab=rownames(fracS)[j], bty="n",ylab=sprintf("PC%s",i))
    legend("topleft", legend=round(corg[j,i],2), bty="n")
  }
}

dev.off()
write.csv(corg, sprintf("%s_ASE_cor_genexPC_short.csv", pref), quote=F)
write.csv(corc, sprintf("%s_ASE_cor_covsxPC_short.csv", pref), quote=F)




pdf(sprintf("%s_ASE_cor_PCs_long.pdf", pref), width=6, height=6)
#long
val = prL$values[1:20]
PC1 =  prL$vectors[,1]
PC2 =  prL$vectors[,2]
PC3 =  prL$vectors[,3]
PC4 =  prL$vectors[,4]

par(mar=c(5,4,1,1), mfrow=c(2,2))
barplot(val, main="", xlab="Index", ylab="Eigen-value, long")
cols = heat.colors(n=6)
ageg = sort(names(table(phenoL$age)))
legend("topright", ageg, bty="n", text.col=cols)
legend("top", c("male", "female"), pch=1:2, bty="n")
#cols = gray.colors(n=6, start = 0, end = 0.7, gamma = 2.2, alpha=1)
coli = cols[factor(phenoL$age, levels=ageg)]
pchi = phenoL$SEX
cexi = phenoL$tp53-min(phenoL$tp53)+.5
#cexi = cexi/max(cexi)
plot(PC1, PC2, cex=cexi, col=coli, pch=pchi, bty="n", main="")
plot(PC1, PC3, cex=cexi, col=coli, pch=pchi, bty="n", main="")
plot(PC2, PC3, cex=cexi, col=coli, pch=pchi, bty="n", main="")



pchi = phenoL$DTHHRDY+1
cexi = phenoL$ctcf-min(phenoL$ctcf)+.5
barplot(val, main="", xlab="Index", ylab="Eigen-value, long")
legend("topright", legend=sort(unique(phenoL$DTHHRDY)), pch=1:5, bty="n")
legend("top", c("male", "female"), pch=1:2, bty="n")
#cexi = cexi/max(cexi)
plot(PC1, PC2, cex=cexi, col=coli, pch=pchi, bty="n", main="")
plot(PC1, PC3, cex=cexi, col=coli, pch=pchi, bty="n", main="")
plot(PC2, PC3, cex=cexi, col=coli, pch=pchi, bty="n", main="")


#for(covi in 1:nrow(covars)){
#  plot(x=unlist(covars[covi,]), y=PC1, main=sprintf("PC1x%s", rownames(covars)[covi]), xlab=rownames(covars)[covi], bty="n")
#  plot(x=unlist(covars[covi,]), y=PC2, main=sprintf("PC2x%s", rownames(covars)[covi]), xlab=rownames(covars)[covi], bty="n")
#  plot(x=unlist(covars[covi,]), y=PC3, main=sprintf("PC3x%s", rownames(covars)[covi]), xlab=rownames(covars)[covi], bty="n")
#  plot(x=unlist(covars[covi,]), y=PC4, main=sprintf("PC4x%s", rownames(covars)[covi]), xlab=rownames(covars)[covi], bty="n")
#}

par(mar=c(5,4,1,1), mfrow=c(2,2))
corc = matrix(NA, nrow=nrow(covars), ncol=20)
corg = matrix(NA, nrow=nrow(fracL), ncol=20)
rownames(corg) = rownames(fracL)
rownames(corc) = rownames(covars)

for(i in 1:nrow(corc)){
y = unlist(covars[i,])
for(j in 1:20){
  corc[i,j] = cor(prL$vectors[,j], y)
}
}
for(i in 1:20){
  x = 1:nrow(corc)
  plot(x, corc[,i], bty="n", xlab="covariates", ylab="correlation", main=sprintf("PC%s", i))
}

for(pci in 1:20){
  kp = which(abs(corc[,pci])>.25)
for(covi in kp){
  plot(x=unlist(covars[covi,]), y=prS$vectors[,pci], main=sprintf("PC%sx%s", pci, rownames(covars)[covi]), xlab=rownames(covars)[covi], bty="n")
}
}

for(i in 1:nrow(corg)){
x = fracL[i,]
for(j in 1:20){
  y = prL$vectors[,j]
  kp = !is.na(x)&!is.na(y)
  corg[i,j] = cor(x[kp], y[kp])
}
}
par(mar=c(5,4,1,1), mfrow=c(2,2))
for(i in 1:20){
  x = 1:nrow(corg)
  #plot(x, corg[,i], bty="n", xlab="covariates", ylab="correlation", main=sprintf("PC%s", i))
  plot(density(corg[,i]), bty="n", ylab="corelations", main=sprintf("PC%s", i))
}
par(mar=c(5,4,1,1), mfrow=c(2,2))
for(i in 1:20){
  kp = which(abs(corg[,i])>.5)
  for(j in kp){
    plot(fracL[j,], prL$vectors[,i], xlab=rownames(fracL)[j], bty="n",ylab=sprintf("PC%s",i))
    legend("topleft", legend=round(corg[j,i],2), bty="n")
  }
}


dev.off()

write.csv(corg, sprintf("%s_ASE_cor_genexPC_long.csv", pref), quote=F)
write.csv(corc, sprintf("%s_ASE_cor_covsxPC_long.csv", pref), quote=F)



q("no")
