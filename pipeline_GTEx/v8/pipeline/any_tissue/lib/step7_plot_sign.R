pheno.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8/Annotations"
pheno = read.table(sprintf("%s/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", pheno.dir), as.is=T, sep="\t", header=T)
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

root.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8/2020_01_20"
dirs = c("Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere",
"Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Nucleus_accumbens_basal_ganglia", "Whole_Blood")
#dirs = c("Whole_Blood")
#conds = c("age", "tp53")
conds = c("age", "tp53", "ctcf")
#conds = c("age")
sams = c(194, 175, 205, 175, 202, 670)
#sams = c(670)
out.dir = sprintf("%s/%s/summ", root.dir, dirs[1])
if(!file.exists(out.dir))dir.create(out.dir)
setwd(out.dir)

splitting = function(x, split=":", block=5, convert=T){
  out = unlist(strsplit(x, split=split))[block]
  if(convert)out = as.numeric(out)
  out
}



library(qvalue)
j = 1; k = 1
for(j in 1:length(dirs)){

pref = dirs[j]
nsam = sams[j]
diri = sprintf("/pine/scr/z/h/zhabotyn/R01/GTEx/v8/2020_01_20/%s", pref)
cntL = read.csv(sprintf("%s/counts_with_min_snp_%s_%s_5e+05_long.csv", diri, pref, nsam), as.is=T)

colnames(cntL) = gsub("\\.", "-", colnames(cntL))
rownames(cntL) = cntL[,1]; cntL = cntL[,-1]; cntL = as.matrix(cntL)
cntL[1:4,1:4]
minsL = apply(cntL, 1:2, splitting)
totL = apply(cntL, 1:2, splitting, block=1)
crv = colSums(totL); crv = median(crv)/crv
nrmL = totL%*%diag(crv);summary(colSums(nrmL))

#tp53, ctcf
cand = c(grep("ENSG00000141510", rownames(cntL)), grep("ENSG00000102974", rownames(cntL)))

hap1L = apply(cntL, 1:2, splitting, block=2)
hap2L = apply(cntL, 1:2, splitting, block=3)

sam.ord = colnames(cntL)
m = match(sam.ord, pheno$SUBJID)
pheno.sub = pheno[m,]
m = match(colnames(hap1L), pheno$SUBJID);
pheni = pheno[m,]
pheni$tp53 = nrmL[cand[1],]
pheni$ctcf = nrmL[cand[2],]
pheni[1:2,]

ind = sprintf("%s/%s", root.dir, dirs[j])

for(k in 1:length(conds)){
flLi = sprintf("%s/%s_%s_long_qp_.csv", ind, dirs[j], conds[k])
flSi = sprintf("%s/%s_%s_short_qp_.csv", ind, dirs[j], conds[k])
rsLi = read.csv(flLi)
rsSi = read.csv(flSi)
kp = apply(is.na(rsLi[,2:9]), 1, any)
rsLi[kp,-1] = NA
kp = apply(is.na(rsSi[,2:9]), 1, any)
rsSi[kp,-1] = NA
colnames(rsLi) = colnames(rsSi) = c("id", "int", "sex", "cond", "sex_cond", "p_int", "p_sex", "p_cond", "p_sex_cond")
write.csv(rsLi[!is.na(rsLi$int),], sprintf("long_quasi_%s_%s.csv", dirs[j], conds[k]), quote=F, row.names=F)
write.csv(rsSi[!is.na(rsSi$int),], sprintf("short_quasi_%s_%s.csv", dirs[j], conds[k]), quote=F, row.names=F)

#eqLi = sprintf("%s/eqtl5l_longSNP.csv", ind)
#eqSi = sprintf("%s/eqtl5s_longSNP.csv", ind)
eqLi = sprintf("%s/TReCASE_%s_%s_5e+05_long.csv", ind, dirs[j], sams[j], as.is=T)
eqSi = sprintf("%s/TReCASE_%s_%s_5e+05_short.csv", ind, dirs[j], sams[j], as.is=T)

eqLTi = sprintf("%s/TReCASE_TReC_%s_%s_5e+05_long.csv", ind, dirs[j], sams[j], as.is=T)
eqLAi = sprintf("%s/TReCASE_ASE_%s_%s_5e+05_long.csv", ind, dirs[j], sams[j], as.is=T)

eqSTi = sprintf("%s/TReCASE_TReC_%s_%s_5e+05_short.csv", ind, dirs[j], sams[j], as.is=T)
eqSAi = sprintf("%s/TReCASE_ASE_%s_%s_5e+05_short.csv", ind, dirs[j], sams[j], as.is=T)

eqLi = read.csv(eqLi)
eqSi = read.csv(eqSi)

eqLTi = read.csv(eqLTi)
eqLAi = read.csv(eqLAi)
eqSTi = read.csv(eqSTi)
eqSAi = read.csv(eqSAi)

eqLTi = eqLTi[,c("MarkerRowID", "Pos", "NBod", "TReC_b", "TReC_Pvalue", "permp")]
eqLAi = eqLAi[,c("MarkerRowID", "Pos", "BBod", "ASE_b", "ASE_Pvalue", "permp")]
eqLi = eqLi[,c("MarkerRowID", "Pos", "NBod", "BBod", "TReC_b", "TReC_Pvalue", 
"ASE_b", "ASE_Pvalue", "Joint_b", "Joint_Pvalue", "trans_Pvalue", "final_Pvalue", "permp")]

eqSTi = eqSTi[,c("MarkerRowID", "Pos", "NBod", "TReC_b", "TReC_Pvalue", "permp")]
eqSAi = eqSAi[,c("MarkerRowID", "Pos", "BBod", "ASE_b", "ASE_Pvalue", "permp")]
eqSi = eqSi[,c("MarkerRowID", "Pos", "NBod", "BBod", "TReC_b", "TReC_Pvalue", 
"ASE_b", "ASE_Pvalue", "Joint_b", "Joint_Pvalue", "trans_Pvalue", "final_Pvalue", "permp")]

rownames(eqSi) = rownames(eqSTi) = rownames(eqSAi) = 
rownames(eqLi) = rownames(eqLTi) = rownames(eqLAi) = rsLi[,1]

write.csv(eqLi[!is.na(eqLi$permp),], sprintf("long_%s_TReCASE.csv", dirs[j]), quote=F)
write.csv(eqLTi[!is.na(eqLTi$permp),], sprintf("long_%s_TReC.csv", dirs[j]), quote=F)
write.csv(eqLAi[!is.na(eqLAi$permp),], sprintf("long_%s_ASE.csv", dirs[j]), quote=F)

write.csv(eqSi[!is.na(eqSi$permp),], sprintf("short_%s_TReCASE.csv", dirs[j]), quote=F)
write.csv(eqSTi[!is.na(eqSTi$permp),], sprintf("short_%s_TReC.csv", dirs[j]), quote=F)
write.csv(eqSAi[!is.na(eqSAi$permp),], sprintf("short_%s_ASE.csv", dirs[j]), quote=F)



sexi = 7
effi = 8
inti = 9
cutoff = 1

#summ = matrix(0, 

pdf(sprintf("%s_%s_count_by_snp.pdf", dirs[j], conds[k]), height=6, width=9)
par(mfrow=c(1,1))
#for(cutoff in c(1, 0.01, 1e-3, 1e-4, 1e-5, 1e-10)){
cutoff = 0.01
kp = !is.na(rsLi[,sexi]) & eqLi$permp<=cutoff
table(kp)
dim(rsLi)
dim(eqLi)
rsLj = rsLi[kp,]
eqLj = eqLi[kp,]
rsLj[1:2,]
eqLj[1:2,]

o = order(rsLj$p_cond);o[1:4]
rsLj = rsLj[o,]
eqLj = eqLj[o,]

#gnj = 1
#m = match(rsLj$id[gnj], inf$Name);m
#inf[m,]

#snpi = read.table(sprintf("%s/genotypei_%s_%s.dat",shortform, inf$chr[m], inf$count[m]), sep="\t", header=T, as.is=T)
#snps = read.table(sprintf("%s/genotypes_%s_%s.dat",shortform, inf$chr[m], inf$count[m]), sep="\t", header=T, as.is=T)
#cnts = read.table(sprintf("%s/counti_long_%s_%s.csv",shortform, inf$chr[m], inf$count[m]), sep=",", as.is=T)


gnj = 1
stp = sum(rsLj$p_cond<.05)
if(stp>100)stp=100
for(gnj in 1:stp){
#m = match(rsLj$id[gnj], inf$Name);m
rsLj[gnj,]
m = match(rsLj$id[gnj], rownames(hap1L));m
hap1L[m,1:5]
hap2L[m,1:5]

      genoL = minsL[m,]
      yL = hap1L[m,]
      nL = yL + hap2L[m,]
      yL[genoL==3] = hap2L[m, genoL==3]

min.cnt = 10
kp = which(genoL %in% c(1,3) & nL>=min.cnt)
aL = (yL+1)/(nL+2)
phenoL = pheni[kp,]
#glmL.qp = glm(cbind(yL[kp], nL[kp])~phenoL$sex+phenoL$agen+phenoL$sex*phenoL$agen, family="quasibinomial")
#summary(glmL.qp)
#par(mfrow=c(2,1))
x = pheni$agen[kp]
if(conds[k]== "tp53")x = pheni$tp53[kp]
if(conds[k]== "ctcf")x = pheni$ctcf[kp]

if(conds[k] == "age"){
  boxplot(aL[kp]~x, bty="n", main=rownames(minsL)[m], xlab=conds[k])
  xf = factor(x, levels=sort(unique(x)))
  points(xf, aL[kp], cex=log10(nL[kp]+2))
}else{
  xf = x
  plot(aL[kp]~x, bty="n", main=rownames(minsL)[m], xlab=conds[k])
}
legend("topleft", sprintf("%s=%s", c("b", "p"), format(c(rsLj$cond[gnj], rsLj$p_cond[gnj]), scientific=T, digits=1)), bty="n")
#plot(pheni$agen[-kp], aL[-kp], cex=log10(nL[-kp]+2), bty="n", main="other counts")
if(gnj%%50==0)message(gnj)
}
dev.off()

}
}
