
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
wrk.dir = sprintf("%s/qb_bb2", tis.dir)
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

pheno.sub$gPC1 = pheno.sub$gPC1/sd(pheno.sub$gPC1)
pheno.sub$gPC2 = pheno.sub$gPC2/sd(pheno.sub$gPC2)
pheno.sub$PF1 = pheno.sub$PF1/sd(pheno.sub$PF1)
pheno.sub$PF2 = pheno.sub$PF2/sd(pheno.sub$PF2)
pheno.sub$PF3 = pheno.sub$PF3/sd(pheno.sub$PF3)
pheno.sub$PF4 = pheno.sub$PF4/sd(pheno.sub$PF4)
pheno.sub$PF5 = pheno.sub$PF5/sd(pheno.sub$PF5)

covi = 2
condi = 1

resL1 = resS1 = list(NULL, NULL, NULL)
p0s = matrix(NA, nrow=21, ncol=2)
rownames(p0s) = 1:21
colnames(p0s) = c("long", "short")
i = 0
pdf(sprintf("%s_ASE_vs_covariates.pdf",pref), height=6, width=6)
for(condi in 1:3){

mod = sprintf("%s_%s", pref, conds[condi]);mod
resL1.bb = read.csv(sprintf("%s_long_vglm_ohet.csv", mod), as.is=T)
resS1.bb = read.csv(sprintf("%s_short_vglm_ohet.csv", mod), as.is=T)

resL1.qp = read.csv(sprintf("%s_long_glmqp_ohet.csv", mod), as.is=T)
resS1.qp = read.csv(sprintf("%s_short_glmqp_ohet.csv", mod), as.is=T)

resL1.bin = read.csv(sprintf("%s_long_glmbin_ohet.csv", mod), as.is=T)
resS1.bin = read.csv(sprintf("%s_short_glmbin_ohet.csv", mod), as.is=T)

nms = c("id", "b.int", "b.od", "b.gPC1", "b.gPC2", "b.PF1", "b.PF2", "b.PF3", "b.PF4", "b.PF5", "b.cnd",
"p.int", "p.od", "p.gPC1", "p.gPC2", "p.PF1", "p.PF2", "p.PF3", "p.PF4", "p.PF5", "p.cnd", "pll.cnd", "logLik", "status")
length(nms)
colnames(resS1.bb) = colnames(resL1.bb) = nms
colnames(resS1.bin) = colnames(resL1.bin) = nms[-c(3,13,24)]
table(resL1.bb$status)
kp = resL1.bb$status==1
summary((resL1.bb$logLik-resL1.bin$logLik)[kp])
summary((resL1.bb$logLik-resL1.bin$logLik)[kp]<.01)
kp = resL1.bb$status==1 & (resL1.bb$logLik-resL1.bin$logLik)>0.01
resL1.bb[which(kp),] = NA
resL1.bb$status[kp] = 3
kp = resL1.bb$status==1 & (resL1.bb$logLik-resL1.bin$logLik)<=0.01
table(kp)
kp = which(kp)
resL1.bb[kp,][1:2,]
resL1.bb[kp, colnames(resL1.bin)] = resL1.bin[kp,]
resL1.bb[kp,][1:2,]

kp = resS1.bb$status==1 & (resS1.bb$logLik-resS1.bin$logLik)>0.01
resS1.bb[which(kp),] = NA
resS1.bb$status[kp] = 3
kp = resS1.bb$status==1 & (resS1.bb$logLik-resS1.bin$logLik)<=0.01
kp = which(kp)
resS1.bb[kp, colnames(resS1.bin)] = resS1.bin[kp,]


par(mfrow=c(2,2))
for(covi in 1:7){
i = i + 1
cnd2nm = colnames(pheno.sub)[9+covi];cnd2nm
pcnd2nm = sprintf("p.%s", cnd2nm);pcnd2nm
kpL = resL1.bb$status %in% 0:1
table(kpL)
pi0w = qvalue(resL1.bb[kpL, pcnd2nm])$pi0
pi0L = 2*mean(resL1.bb[kpL, pcnd2nm] >= 0.5, na.rm=T)
pi0L[pi0L>1]=1
c(pi0w, pi0L)

kpS = resS1.bb$status %in% 0:1
pi0w = qvalue(resS1.bb[kpS,pcnd2nm])$pi0
pi0S = 2*mean(resS1.bb[kpS,pcnd2nm] >= 0.5, na.rm=T)
pi0S[pi0S>1]=1
c(pi0w, pi0S)

hist(resL1.bb[kpL, pcnd2nm], main=sprintf("%s&%s, eQTL from long", conds[condi], cnd2nm), xlab="p-value")
legend("topright", legend=sprintf("prop.null=%s", round(pi0L,2)), bty="n")
hist(resS1.bb[kpS, pcnd2nm], main=sprintf("%s&%s, eQTL from short", conds[condi], cnd2nm), xlab="p-value")
legend("topright", legend=sprintf("prop.null=%s", round(pi0S,2)), bty="n")
p0s[i,] = c(pi0L, pi0S)
rownames(p0s)[i] = sprintf("%s.%s", conds[condi], cnd2nm)
}
message(condi)

resL1.bb = resL1.bb[kpL,]
resS1.bb = resS1.bb[kpS,]
write.csv(resL1.bb, sprintf("%s_long_final.csv", mod), row.names=F)
write.csv(resS1.bb, sprintf("%s_short_final.csv", mod), row.names=F)
resL1[[condi]] = resL1.bb
resS1[[condi]] = resS1.bb
}
dev.off()
write.csv(p0s, sprintf("%s_p0s.csv", pref), quote=F)



pdf(sprintf("%s_ASE_coi_wald.pdf",pref), height=6, width=9)
pcnd2nm = "p.cnd"
par(mfrow=c(2,3))
for(condi in 1:3){
  resL1.bb = resL1[[condi]]
  pi0L = 2*mean(resL1.bb[,pcnd2nm] >= 0.5, na.rm=T)
  qvals = qvalue(resL1.bb[,pcnd2nm])$qvalues
  legL = c(round(pi0L,2), sum(qvals<.1,na.rm=T), sum(qvals<.25, na.rm=T))
  legL = sprintf("%s=%s", c("prop.null", "q0.10", "q.25"), legL)
  hist(resL1.bb[, pcnd2nm], main=sprintf("%s, eQTL from long", conds[condi]), xlab="p-value", cex=2.2, cex.main=2.2, cex.axis=2.2, cex.lab=2.2)
  legend("topright", legend=legL, bty="n")
}  
for(condi in 1:3){
  resS1.bb = resS1[[condi]]
  pi0S = 2*mean(resS1.bb[,pcnd2nm] >= 0.5, na.rm=T)
  qvals = qvalue(resS1.bb[,pcnd2nm])$qvalues
  legS = c(round(pi0S,2), sum(qvals<.1,na.rm=T), sum(qvals<.25, na.rm=T))
  legS = sprintf("%s=%s", c("prop.null", "q0.10", "q.25"), legS)
  hist(resS1.bb[, pcnd2nm], main=sprintf("%s, eQTL from short", conds[condi]), xlab="p-value", cex=2.2, cex.main=2.2, cex.axis=2.2, cex.lab=2.2)
  legend("topright", legend=legS, bty="n")
}
dev.off()

summL = matrix(NA, nrow=3, ncol=8)
summS = matrix(NA, nrow=3, ncol=8)
rownames(summS) = rownames(summL) = conds
colnames(summS) = colnames(summL) = c("ngen", "p0", "q.01", "q.05", "q.10", "q.15", "q.20", "q.25")
pdf(sprintf("%s_ASE_coi_lrt.pdf",pref), height=6, width=9)
pcnd2nm = "pll.cnd"
par(mfrow=c(2,3), mar=c(5,5,4,1)
for(condi in 1:3){
  resL1.bb = resL1[[condi]]
  pi0L = 2*mean(resL1.bb[,pcnd2nm] >= 0.5, na.rm=T)
  qvals = qvalue(resL1.bb[,pcnd2nm])$qvalues
  legL = c(nrow(resL1.bb), pi0L, sum(qvals<.01,na.rm=T), sum(qvals<.05,na.rm=T), 
  sum(qvals<.1,na.rm=T), sum(qvals<.15,na.rm=T), sum(qvals<.20, na.rm=T), sum(qvals<.25, na.rm=T))
  summL[condi,] = legL
  legL = c(pi0L, sum(qvals<.1,na.rm=T), sum(qvals<.25, na.rm=T))
  legL = sprintf("%s=%s", c("prop.null", "q0.10", "q.25"), round(legL,2))
  hist(resL1.bb[, pcnd2nm], main=sprintf("%s, eQTL from long", conds[condi]), xlab="p-value", cex=2.2, cex.main=2.2, cex.axis=2.2, cex.lab=2.2)
  legend("topright", legend=legL, bty="n")  
}  
for(condi in 1:3){
  resS1.bb = resS1[[condi]]
  pi0S = 2*mean(resS1.bb[,pcnd2nm] >= 0.5, na.rm=T)
  qvals = qvalue(resS1.bb[,pcnd2nm])$qvalues
  legS = c(nrow(resS1.bb), pi0S, sum(qvals<.01,na.rm=T), sum(qvals<.05,na.rm=T), 
  sum(qvals<.1,na.rm=T), sum(qvals<.15,na.rm=T), sum(qvals<.20, na.rm=T), sum(qvals<.25, na.rm=T))
  summS[condi,] = legS
  legS = c(pi0S, sum(qvals<.1,na.rm=T), sum(qvals<.25, na.rm=T))
  legS = sprintf("%s=%s", c("prop.null", "q0.10", "q.25"), round(legS,2))
  hist(resS1.bb[, pcnd2nm], main=sprintf("%s, eQTL from short", conds[condi]), xlab="p-value", cex=2.2, cex.main=2.2, cex.axis=2.2, cex.lab=2.2)
  legend("topright", legend=legS, bty="n")
}
dev.off()
summL
summS
write.csv(summL, sprintf("%s_sig_cnts.csv", pref))
summL = read.csv(sprintf("%s_sig_cnts.csv", pref))
write.csv(summS, sprintf("%s_sig_cnts.csv", pref))
summS = read.csv(sprintf("%s_sig_cnts.csv", pref))
summL
summS
write.table(summL, sprintf("%s_sig_cnts.csv", pref), quote=F, sep=",", col.names=T, row.names=F)
write.table(summS, sprintf("%s_sig_cnts.csv", pref), quote=F, append=T, sep=",", col.names=T, row.names=F)



q("no")
