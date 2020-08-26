library(qvalue)
library(VGAM)

specf = "specifications.txt"
getwd()

specs = unlist(read.table(specf, as.is=T))
specs
pref = specs[1]
nsam = specs[2]
queue = specs[3]
days = specs[4]
bmem = as.numeric(specs[5])
seedval = specs[13]
wrk.dir = specs[14]
lib.dir = specs[15]
bas.dir = specs[16]
setwd(wrk.dir)

conds = c("age", "tp53", "ctcf")
v8.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8"
tis.dir = sprintf("%s/2020_01_20/%s", v8.dir, pref)
wrk.dir = sprintf("%s/qb_bb1", tis.dir)
setwd(wrk.dir)

ann.dir = sprintf("%s/Annotations", v8.dir)
cov.dir = sprintf("%s/GTEx_Analysis_v8_eQTL_covariates", ann.dir)
cov.fil = sprintf("%s/%s.v8.covariates.txt",cov.dir, pref)
covars = read.table(cov.fil, header=T, as.is=T)
rownames(covars) = covars$ID; covars = covars[,-1]
colnames(covars) = gsub("\\.", "-", colnames(covars))
covars[1:2,1:5]

phen.fil = sprintf("%s/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", ann.dir)
pheno = read.table(phen.fil, as.is=T, sep="\t", header=T)
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

lonf = sprintf("%s/short_TReCASE_%s_%s_5e+05_long.csv", tis.dir, pref, nsam)
shof = sprintf("%s/short_TReCASE_%s_%s_5e+05_short.csv", tis.dir, pref, nsam)
resL = read.csv(lonf, as.is=T)
resS = read.csv(shof, as.is=T)

splitting = function(x, split=":", block=5, convert=T){
  out = unlist(strsplit(x, split=split))[block]
  if(convert)out = as.numeric(out)
  out
}

lonf = sprintf("%s/counts_with_min_snp_%s_%s_5e+05_long.csv", tis.dir, pref, nsam)
shof = sprintf("%s/counts_with_min_snp_%s_%s_5e+05_short.csv", tis.dir, pref, nsam)
cntL = read.csv(lonf, as.is=T)
cntS = read.csv(shof, as.is=T)
colnames(cntL) = gsub("\\.", "-", colnames(cntL))
colnames(cntS) = gsub("\\.", "-", colnames(cntS))
rownames(cntL) = cntL[,1]; cntL = cntL[,-1]; cntL = as.matrix(cntL)
rownames(cntS) = cntS[,1]; cntS = cntS[,-1]; cntS = as.matrix(cntS)
cntL[1:4,1:4]
cntS[1:4,1:4]

cand = c(grep("ENSG00000141510", rownames(cntL)), grep("ENSG00000102974", rownames(cntL)))
totL = apply(cntL, 1:2, splitting, block=1)
crv = colSums(totL); crv = median(crv)/crv
nrmL = totL%*%diag(crv);summary(colSums(nrmL))

sam.ord = colnames(cntL)
m = match(sam.ord, pheno$SUBJID)
table(pheno$SUBJID[m] == sam.ord)
pheno.sub = pheno[m,]
pheno.sub$age = pheno.sub$agen
pheno.sub$tp53 = log10(nrmL[cand[1],]+1)
pheno.sub$ctcf = log10(nrmL[cand[2],]+1)
cutoff = 0.01

m = match(sam.ord, colnames(covars))
covars = covars[,m]
table(sam.ord==colnames(covars))
pheno.sub$gPC1 = unlist(covars["PC1",])
pheno.sub$gPC2 = unlist(covars["PC2",])
pheno.sub$PF1 = unlist(covars["InferredCov1",])
pheno.sub$PF2 = unlist(covars["InferredCov2",])
pheno.sub$PF3 = unlist(covars["InferredCov3",])
pheno.sub$PF4 = unlist(covars["InferredCov4",])
pheno.sub$PF5 = unlist(covars["InferredCov5",])
pheno.sub$PF6 = unlist(covars["InferredCov6",])
pheno.sub$PF7 = unlist(covars["InferredCov7",])
pheno.sub$PF8 = unlist(covars["InferredCov8",])
pheno.sub$PF9 = unlist(covars["InferredCov9",])
pheno.sub$PF10 = unlist(covars["InferredCov10",])

pheno.sub$gPC1 = pheno.sub$gPC1/sd(pheno.sub$gPC1)
pheno.sub$gPC2 = pheno.sub$gPC2/sd(pheno.sub$gPC2)
pheno.sub$PF1 = pheno.sub$PF1/sd(pheno.sub$PF1)
pheno.sub$PF2 = pheno.sub$PF2/sd(pheno.sub$PF2)
pheno.sub$PF3 = pheno.sub$PF3/sd(pheno.sub$PF3)
pheno.sub$PF4 = pheno.sub$PF4/sd(pheno.sub$PF4)
pheno.sub$PF5 = pheno.sub$PF5/sd(pheno.sub$PF5)
pheno.sub$PF6 = pheno.sub$PF5/sd(pheno.sub$PF6)
pheno.sub$PF7 = pheno.sub$PF5/sd(pheno.sub$PF7)
pheno.sub$PF8 = pheno.sub$PF5/sd(pheno.sub$PF8)
pheno.sub$PF9 = pheno.sub$PF5/sd(pheno.sub$PF9)
pheno.sub$PF10 = pheno.sub$PF5/sd(pheno.sub$PF10)

covi = 2
condi = 1

p0s = matrix(NA, nrow=12, ncol=2)
rownames(p0s) = 1:nrow(p0s)
colnames(p0s) = c("long", "short")
i = 0
pdf(sprintf("%s_ASE_vs_covariates.pdf",pref), height=6, width=6)
#for(condi in 1:3){
par(mfrow=c(2,2))
for(covi in 1:12){
i = i + 1
cnd2nm = colnames(pheno.sub)[9+covi];cnd2nm

mod = sprintf("%s_%s", pref, cnd2nm);mod
resL1.bb = read.csv(sprintf("%s_long_vglm_ohet.csv", mod), as.is=T)
resS1.bb = read.csv(sprintf("%s_short_vglm_ohet.csv", mod), as.is=T)

#resL1.qp = read.csv(sprintf("%s_long_glmqp_ohet.csv", mod), as.is=T)
#resS1.qp = read.csv(sprintf("%s_short_glmqp_ohet.csv", mod), as.is=T)

resL1.bin = read.csv(sprintf("%s_long_glmbin_ohet.csv", mod), as.is=T)
resS1.bin = read.csv(sprintf("%s_short_glmbin_ohet.csv", mod), as.is=T)

resL1.bb[1:2,]
table(resL1.bb$status)
kp = resL1.bb$status==1
summary((resL1.bb$logLik-resL1.bin$logLik)[kp])
kp = resL1.bb$status==1 & (resL1.bb$logLik-resL1.bin$logLik)>0.01
resL1.bb$status[kp] = 3
kp = resL1.bb$status==1 & (resL1.bb$logLik-resL1.bin$logLik)<=0.01
table(kp)
kp = which(kp)
resL1.bb[kp,][1:2,]
resL1.bb[kp, colnames(resL1.bin)] = resL1.bin[kp,]
resL1.bb[kp,][1:2,]
kp = resL1.bb$status %in% 0:1
table(kp)
summary(resL1.bb$p.cnd2)
summary(resL1.bb$pll.cnd2)
pi0w = qvalue(resL1.bb$p.cnd2[kp])$pi0
pi0l = qvalue(resL1.bb$pll.cnd2[kp])$pi0
pi0L = 2*mean(resL1.bb$pll.cnd2[kp] >= 0.5, na.rm=T)
pi0L[pi0L>1]=1
c(pi0w, pi0l, pi0L)

kp = resS1.bb$status==1
kp = resS1.bb$status==1 & (resS1.bb$logSik-resS1.bin$logSik)>0.01
resS1.bb$status[kp] = 3
kp = resS1.bb$status==1 & (resS1.bb$logSik-resS1.bin$logSik)<=0.01
kp = which(kp)
resS1.bb[kp, colnames(resS1.bin)] = resS1.bin[kp,]
kp = resS1.bb$status %in% 0:1
pi0w = qvalue(resS1.bb$p.cnd2[kp])$pi0
pi0l = qvalue(resS1.bb$pll.cnd2[kp])$pi0
pi0S = 2*mean(resS1.bb$pll.cnd2[kp] >= 0.5, na.rm=T)
pi0S[pi0S>1]=1
c(pi0w, pi0l, pi0S)

hist(resL1.bb$pll.cnd2[kp], main=sprintf("%s, eQTL from long", cnd2nm), xlab="p-value")
legend("topright", legend=sprintf("prop.null=%s", round(pi0L,2)), bty="n")
hist(resS1.bb$pll.cnd2[kp], main=sprintf("%s, eQTL from short", cnd2nm), xlab="p-value")
legend("topright", legend=sprintf("prop.null=%s", round(pi0S,2)), bty="n")
p0s[i,] = c(pi0L, pi0S)
rownames(p0s)[i] = sprintf("%s", cnd2nm)
#}
message(covi)
}
dev.off()
p0s

write.csv(p0s, sprintf("%s_p0s.csv", pref), quote=F)



i = 0
pdf(sprintf("%s_ASE_vs_covariates_long.pdf",pref), height=9, width=12)
#for(condi in 1:3){
par(mfrow=c(2,3), mar=c(5, 5, 4, 1))
for(covi in 1:12){
i = i + 1
cnd2nm = colnames(pheno.sub)[9+covi];cnd2nm

mod = sprintf("%s_%s", pref, cnd2nm);mod
resL1.bb = read.csv(sprintf("%s_long_vglm_ohet.csv", mod), as.is=T)
resL1.bin = read.csv(sprintf("%s_long_glmbin_ohet.csv", mod), as.is=T)

resL1.bb[1:2,]
table(resL1.bb$status)
kp = resL1.bb$status==1
summary((resL1.bb$logLik-resL1.bin$logLik)[kp])
kp = resL1.bb$status==1 & (resL1.bb$logLik-resL1.bin$logLik)>0.01
resL1.bb$status[kp] = 3
kp = resL1.bb$status==1 & (resL1.bb$logLik-resL1.bin$logLik)<=0.01
table(kp)
kp = which(kp)
resL1.bb[kp,][1:2,]
resL1.bb[kp, colnames(resL1.bin)] = resL1.bin[kp,]
resL1.bb[kp,][1:2,]
kp = resL1.bb$status %in% 0:1
pi0w = qvalue(resL1.bb$p.cnd2[kp])$pi0
pi0l = qvalue(resL1.bb$pll.cnd2[kp])$pi0
pi0L = 2*mean(resL1.bb$pll.cnd2[kp] >= 0.5, na.rm=T)
pi0L[pi0L>1]=1

hist(resL1.bb$pll.cnd2[kp], main=sprintf("%s", cnd2nm), xlab="p-value", cex=2.2, cex.main=2.2, cex.axis=2.2, cex.lab=2.2)
legend("topright", legend=sprintf("prop.null=%s", round(pi0L,2)), bty="n", cex=2.2)
#}
message(covi)
}
dev.off()


q("no")
