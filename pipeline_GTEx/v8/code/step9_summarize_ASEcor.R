#rt.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8/2020_06_20"
rt.dir = "/home/groups/projects/Sun_RNA_seq/GTEx"
rs.dir = sprintf("%s/plots", rt.dir)
if(!file.exists(rs.dir))dir.create(rs.dir)
setwd(rs.dir)
library(qvalue)

dirs =
c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Artery_Aorta",
"Artery_Coronary", "Artery_Tibial", "Brain_Caudate_basal_ganglia",
"Brain_Cerebellar_Hemisphere", "Brain_Cortex", "Brain_Frontal_Cortex_BA9",
"Brain_Nucleus_accumbens_basal_ganglia", "Breast_Mammary_Tissue",
"Cells_Cultured_fibroblasts", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction",
"Esophagus_Mucosa", "Esophagus_Muscularis", 
"Heart_Atrial_Appendage", "Heart_Left_Ventricle",
"Lung", "Muscle_Skeletal", "Nerve_Tibial",
"Pancreas", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", 
"Stomach", "Testis", "Thyroid", "Whole_Blood")
sams = c(579,467,386,212,581,194,175,205,175,202,396,482,366,329,495,463,371,385,513,704,532,303,514,603,323,321,573,670)
#sams = c(579,467,386,212,581,194,175,205,175,202,396,482,366,329,495,463,513,704,532,303,514,603,323,321,573,670)
dirs[17:18]
sams[17:18]
nams = 
c("Adipose Subcutaneous", "Adipose Visc. Omentum", "Artery Aorta",
"Artery Coronary", "Artery Tibial", "Brain Caudate bas.gang.",
"Brain Cereb. Hemis.", "Brain Cortex", "Brain Frontal Cortex BA9",
"Brain Nucl.acc.bas.gang.", "Breast Mammary Tissue",
"Cells Cultured fibroblasts", "Colon Transverse", "Esophagus Gastr. Junct.",
"Esophagus Mucosa", "Esophagus Muscularis", 
"Heart Atrial App.", "Heart Left Ventr.",
"Lung", "Muscle Skeletal", "Nerve Tibial",
"Pancreas", "Skin not Sun exp. Suprap.", "Skin Sun exp. Lower leg", 
"Stomach", "Testis", "Thyroid", "Whole Blood")
length(dirs)
length(nams)

#using GTEx colors
tiscol = c(
rgb(255, 102, 0, maxColorValue = 255),#Ad. Sub  -     - #FF6600
rgb(255, 170, 0, maxColorValue = 255),#Ad. Visc -     - #FFAA00
rgb(255, 85, 85, maxColorValue = 255),#Ar. Aort -     - #FF5555
rgb(255, 170, 153, maxColorValue = 255),#Ar. Coro - 	- #FFAA99
rgb(255, 0, 0, maxColorValue = 255),#Ar. Tib -   - #FF0000
rgb(238, 238, 0, maxColorValue = 255),#Br       -     - #EEEE00
rgb(238, 238, 0, maxColorValue = 255),#Br       -     - #EEEE00
rgb(238, 238, 0, maxColorValue = 255),#Br       -     - #EEEE00
rgb(238, 238, 0, maxColorValue = 255),#Br       -     - #EEEE00
rgb(238, 238, 0, maxColorValue = 255),#Br       -     - #EEEE00
rgb(51,204,204, maxColorValue = 255),#Brst.MT  -  - #33CCCC
rgb(170, 238, 255, maxColorValue = 255),#CCf      -  - #AAEEFF
rgb(204, 153, 85, maxColorValue = 255),#CT       -  - #CC9955
rgb(136, 115, 85, maxColorValue = 255),#EGJ -  - #887355 - dark orange
rgb(85, 34, 0, maxColorValue = 255),#EMuc -  - #552200 - dark orange
rgb(136, 153, 136, maxColorValue = 255),#EMus -  - #889988 - light green
rgb(153, 0, 255, maxColorValue = 255),#HAA -  - #9900FF - dark violet
rgb(102, 0, 153, maxColorValue = 255),#HLV - #660099 - dark violet
rgb(153, 255, 0, maxColorValue = 255), #Lung - #99FF00 - light green
rgb(170, 170, 255, maxColorValue = 255),#MS , #AAAAFF - light blue
rgb(255, 215, 0, maxColorValue = 255), #NT- - #FFD700 - light yellow
rgb(153, 85, 34, maxColorValue = 255),#Panc - rgb() - #995522 - dark orange
rgb(0, 0, 255, maxColorValue = 255),#SNSES - rgb() - #0000FF - dark blue
rgb(119, 119, 255, maxColorValue = 255),#SSELl - rgb() - #7777FF - light blue
rgb(255, 221, 153, maxColorValue = 255),#Stomach - rgb() - #FFDD99 - light orange
rgb(170,170,170, maxColorValue = 255), #Testis - rgb() - #AAAAAA - light red
rgb(0,102,0, maxColorValue = 255),#Thyroid - rgb() - #006600 - dark green
rgb(255,0,187, maxColorValue = 255))#Whole Blood - rgb() - #FF00BB - dark pink
#tiscol = tiscol[-c(17,18)]
pdf(sprintf("%s/GTEx_color_scheme.pdf",rs.dir), height=6, width=6)
plot(0:1, 0:1, col="white", bty="n", xaxt="n", yaxt="n", main="",xlab="", ylab="")
legend("topleft", nams[1:13], text.col=tiscol[1:13], bty="n")
legend("topright", nams[14:26], text.col=tiscol[14:26], bty="n")
dev.off()
pchs = rep(19, length(tiscol))
ind = c(1, 2,11, 12,13,14,16,17,18,19,23,24,26)
ind = c(2,11, 12,13,16,18,19,20,25,26,28)
ind2 = c(15,21,23,24,27)
nams[ind]
pchs[ind] = 15
pchs[ind2] = 17
modL = "long"
modS = "short"
i = 1
nams[pchs==15]
nams[pchs==17]
eigS = eigL = matrix(NA, nrow=length(dirs), ncol=21)
dim(eigL)

pdf("top_PCs_long.pdf", height= 3*4, width=3*7)
par(mfrow=c(4,7), mar=c(2,2,2,0))
for(i in 1:length(dirs)){
#i = i + 1
pref = dirs[i]
nam = nams[i]
fiLi = sprintf("%s/%s/%s_ASE_cors_%s.csv", rt.dir, pref, pref, modL)
effL = read.csv(fiLi, as.is=T); effL = effL[,-1]

prL  = eigen(effL)
prL$values[1:20]
eigL[i,1:20] = val = prL$values[1:20]
eigL[i,21] = sum(prL$values)
barplot(val, main=nam, cex=2, cex.names=2, cex.main=2, xaxt="n", yaxt="n")#, xlab="Index", ylab="Eigen-value"
message(i)
}
dev.off()

i = 1
pdf("top_PCs_short.pdf", height= 3*4, width=3*7)
par(mfrow=c(4,7), mar=c(2,2,2,0))
for(i in 1:length(dirs)){
#i = i + 1
pref = dirs[i]
nam = nams[i]
fiLi = sprintf("%s/%s/%s_ASE_cors_%s.csv", rt.dir, pref, pref, modS)
effL = read.csv(fiLi, as.is=T); effL = effL[,-1]

prL  = eigen(effL)
prL$values[1:20]
eigS[i,1:20] = val = prL$values[1:20]
eigS[i,21] = sum(prL$values)
barplot(val, main=nam, cex=2, cex.names=2, cex.main=2, xaxt="n", yaxt="n")#, xlab="Index", ylab="Eigen-value"
message(i)
}
dev.off()

pdf("PCs_long.pdf", height= 3, width=6)
par(mfrow=c(1,2))
x = eigL[,1]
y = eigL[,2]
plot(x, y, xlab="eigen 1", ylab="eigen 2", bty="n", pch=19, col=tiscol)
abline(a=0, b=1, col="red", lty=3)
x = eigL[,1]
y = rowSums(eigL[,2:20])
plot(x, y, xlab="eigen 1", ylab="eigen 2-20", bty="n", pch=19, col=tiscol)
abline(a=0, b=1, col="red", lty=3)
dev.off()

pdf("PCs_short.pdf", height= 3, width=6)
par(mfrow=c(1,2))
x = eigS[,1]
y = eigS[,2]
plot(x, y, xlab="eigen 1", ylab="eigen 2", bty="n", pch=19, col=tiscol)
abline(a=0, b=1, col="red", lty=3)
x = eigS[,1]
y = rowSums(eigS[,2:20])
plot(x, y, xlab="eigen 1", ylab="eigen 2-20", bty="n", pch=19, col=tiscol)
abline(a=0, b=1, col="red", lty=3)
dev.off()

pdf("rel_PCs_long.pdf", height= 3, width=6)
par(mfrow=c(1,2))
x = eigL[,1]
y = eigL[,2]
plot(x, x/y, xlab="eigen 1", ylab="ratio", main="eigen 1/eigen 2", bty="n", pch=19, col=tiscol)
x = eigL[,1]
y = rowSums(eigL[,2:20])
plot(x, x/y, xlab="eigen 1", ylab="ratio", main="eigen 1/eigen 2-20", bty="n", pch=19, col=tiscol)
dev.off()

pdf("rel_PCs_short.pdf", height= 3, width=6)
par(mfrow=c(1,2))
x = eigS[,1]
y = eigS[,2]
plot(x, x/y, xlab="eigen 1", ylab="ratio", main="eigen 1/eigen 2", bty="n", pch=19, col=tiscol)
x = eigS[,1]
y = rowSums(eigS[,2:20])
plot(x, x/y, xlab="eigen 1", ylab="ratio", main="eigen 1/eigen 2-20", bty="n", pch=19, col=tiscol)
dev.off()

pdf("frac_PCs.pdf", height= 3, width=6)
par(mfrow=c(1,2))
x = eigL[,1]
y = rowSums(eigL[,1:20])
plot(x, x/y, xlab="eig.1", ylab="frac. first 20 eig", main="long", bty="n", pch=19, col=tiscol)

x = eigS[,1]
y = rowSums(eigS[,1:20])
plot(x, x/y, xlab="eig.1", ylab="frac. first 20 eig", main="short", bty="n", pch=19, col=tiscol)
dev.off()





#dirs =
#c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Artery_Aorta",
#"Artery_Coronary", "Artery_Tibial", "Brain_Caudate_basal_ganglia",
#"Brain_Cerebellar_Hemisphere", "Brain_Cortex", "Brain_Frontal_Cortex_BA9",
#"Brain_Nucleus_accumbens_basal_ganglia", "Breast_Mammary_Tissue",
#"Cells_Cultured_fibroblasts", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction",
#"Esophagus_Mucosa", "Esophagus_Muscularis", "Lung", "Muscle_Skeletal", "Nerve_Tibial",
#"Pancreas", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", 
#"Stomach", "Testis", "Thyroid", "Whole_Blood")
#/pine/scr/z/h/zhabotyn/R01/GTEx/v8/2020_06_20/Adipose_Visceral_Omentum/BetaBin_nocov
i = 1
kps = numeric(0)
for(i in 1:length(dirs)){
#i = 2
pref = dirs[i]
#BetaBin_nocov
#Adipose_Visceral_Omentum_age_long_glmbin_gPC_PF.csv
fiLi = sprintf("%s/%s/BetaBin_gPC_PF/%s_age_%s_glmbin_gPC_PF.csv", rt.dir, pref, pref, modL)
kp = file.exists(fiLi)
if(!kp){
  message(i)
  kps = c(kps, i)
}
}
kps
dir0 = dirs
if(length(kps)>0)dir0 = dirs[-kps]
dir0
length(dir0)
conds = c("age", "ctcf", "tp53")
condi = 1


clnms1 = c("nm", "int", "od", "b.gPC1", "b.gPC2", "b.PF1", "b.PF2", "b.PF3", "b.PF4", "b.PF5", "b.cnd",
"e.int", "e.od", "e.gPC1", "e.gPC2", "e.PF1", "e.PF2", "e.PF3", "e.PF4", "e.PF5", "e.cnd",
"pval","ll","status");length(clnms1)
ind1 = (2:11)[-2];ind1a = 2:10
ind2 = (12:21)[-2];ind2a = 11:19
ind3 = 22;ind3a = 20
ind4 = 23;ind4a = 21

clnms2 = c("nm", "int", "od", "b.gPC1", "b.gPC2", "cnd", "b.PF1", "b.PF2", "b.PF3", "b.PF4", "b.PF5",
"e.int", "e.od", "e.gPC1", "e.gPC2", "e.cnd", "e.PF1", "e.PF2", "e.PF3", "e.PF4", "e.PF5", 
"pval","ll","status");length(clnms1)
ind1 = (2:6)[-2];ind1a = 2:5
ind2 = (12:16)[-2];ind2a = 11:14
ind3 = 22;ind3a = 20
ind4 = 23;ind4a = 21

clnms3 = c("nm", "int", "od", "cnd", "b.gPC1", "b.gPC2", "b.PF1", "b.PF2", "b.PF3", "b.PF4", "b.PF5",
"e.int", "e.od", "e.cnd", "e.gPC1", "e.gPC2", "e.PF1", "e.PF2", "e.PF3", "e.PF4", "e.PF5", 
"pval","ll","status");length(clnms1)
ind1 = (2:4)[-2];ind1a = 2:3
ind2 = (12:14)[-2];ind2a = 11:12
ind3 = 22;ind3a = 20
ind4 = 23;ind4a = 21

alpha0 = 0.01
alpha0b = 0.001
alpha1 = 0.05          
alpha2 = 0.1          
alpha3 = 0.25

#for(i in 1:length(dir0)){
#  pref = dir0[i]
#  bb.cl = sprintf("%s/%s/BetaBin", rt.dir, pref)
#  system(sprintf("rm -r %s/*", bb.cl))
#}

mod = modL
for(condi in 1:3){
  for(mod in c(modL, modS)){
    sig.r = sig.r01 = sig.r001 = matrix(NA, nrow=length(dir0), ncol=12)
    for(i in 1:length(dir0)){
      #i = 2
      pref = dir0[i]
      #BetaBin_nocov
      
      mod1 = mod
#      if(mod1=="short")mod1="fix_short"
      eqL = sprintf("%s/%s/%s_%s_5e+05_%s_TReCASE_trimmed.csv", rt.dir, pref, pref, sams[i], mod1)
      #eqS = sprintf("%s/%s/%s_%s_5e+05_%s_TReCASE_trimmed.csv", rt.dir, pref, pref, sams[i], modS)
      eqL = read.csv(eqL, as.is=T)
      #eqS = read.csv(eqS, as.is=T)
      kpL = eqL[eqL$permp<alpha0, 1]
      kpLb = eqL[eqL$permp<alpha0b, 1]
      #kpS = eqS[eqS$permp<alpha0, 1]
      
      bb.cl = sprintf("%s/%s/BetaBin", rt.dir, pref)
      if(!file.exists(bb.cl))dir.create(bb.cl)
      #suffS = sprintf("%s_%s_%s", pref, conds[condi], modS)
      suff = sprintf("%s_%s_%s", pref, conds[condi], mod)
      fiL1c = sprintf("%s/%s_gPC_PF_all.csv", bb.cl, suff)
      fiL2c = sprintf("%s/%s_gPC_all.csv", bb.cl, suff)
      fiL3c = sprintf("%s/%s_nocov_all.csv", bb.cl, suff)

      fiL1c01 = sprintf("%s/%s_gPC_PF_01.csv", bb.cl, suff)
      fiL2c01 = sprintf("%s/%s_gPC_01.csv", bb.cl, suff)
      fiL3c01 = sprintf("%s/%s_nocov_01.csv", bb.cl, suff)

      fiL1c001 = sprintf("%s/%s_gPC_PF_001.csv", bb.cl, suff)
      fiL2c001 = sprintf("%s/%s_gPC_001.csv", bb.cl, suff)
      fiL3c001 = sprintf("%s/%s_nocov_001.csv", bb.cl, suff)

      fiL1 = sprintf("%s/%s/BetaBin_gPC_PF/%s_vglm_gPC_PF.csv", rt.dir, pref, suff)
      fiL2 = sprintf("%s/%s/BetaBin_gPC/%s_vglm_gPC.csv", rt.dir, pref, suff)
      fiL3 = sprintf("%s/%s/BetaBin_nocov/%s_vglm_nocov.csv", rt.dir, pref, suff)
      #fiL2 = sprintf("%s/BetaBin_nocov/%s_age_%s_vglm_nocov.csv", pref, pref, modL)
      bbr1 = read.csv(fiL1, as.is=T)
      bbr2 = read.csv(fiL2, as.is=T)
      bbr3 = read.csv(fiL3, as.is=T)
      
      fiL1 = sprintf("%s/%s/BetaBin_gPC_PF/%s_glmbin_gPC_PF.csv", rt.dir, pref, suff)
      fiL2 = sprintf("%s/%s/BetaBin_gPC/%s_glmbin_gPC.csv", rt.dir, pref, suff)
      fiL3 = sprintf("%s/%s/BetaBin_nocov/%s_glmbin_nocov.csv", rt.dir, pref, suff)
      #fiL2 = sprintf("%s/BetaBin_nocov/%s_age_%s_glmbin_nocov.csv", pref, pref, modL)
      #fiL2 = sprintf("%s/BetaBin_nocov/%s_age_%s_glmbin_nocov.csv", pref, pref, modL)
      bnr1 = read.csv(fiL1, as.is=T)
      bnr2 = read.csv(fiL2, as.is=T)
      bnr3 = read.csv(fiL3, as.is=T)
      
      #large model
      chk = which(bbr1[,24]==0)
      #summary(abs(bbr1[chk,23]-bnr1[chk,21])<.1)
      chk = which(bbr1[,24]==1 & is.na(abs(bbr1[,23]-bnr1[,21])))
      #length(chk)
      #bbr1[chk,]
      #bnr1[chk,]
      bbr1[chk,24]=4
      chk = which(bbr1[,24]==1 & abs(bbr1[,23]-bnr1[,21])>=.01)
      bbr1[chk,24]=3
      chk = which(bbr1[,24]==1 & abs(bbr1[,23]-bnr1[,21])<=.01)
      length(chk)
      bbr1[chk,3] = NA
      bbr1[chk,13] = NA
      colnames(bbr1) = clnms1
      ind0 = 1
      ind1b=(2:11);ind1 = ind1b[-2];ind1a = 2:10
      ind2b=(12:21);ind2 = ind2b[-2];ind2a = 11:19
      ind3 = 22;ind3a = 20
      ind4 = 23;ind4a = 21
      ind5 = 24
      bbr1[chk,ind1] = bnr1[chk,ind1a]
      bbr1[chk,ind2] = bnr1[chk,ind2a]
      bbr1[chk,ind3] = bnr1[chk,ind3a]
      bbr1[chk,ind4] = bnr1[chk,ind4a]
      table(bbr1[,24])
      bbr1 = bbr1[bbr1$status %in% 0:1,c(ind0, ind1b, ind2b, ind3, ind5)]
      qv = qvalue(bbr1$pval)
      bbr1$qval = qv$qvalue
      bbr1[1:5,]
      sig.r[i, 1] = sum(bbr1$qval<alpha1)
      sig.r[i, 2] = sum(bbr1$qval<alpha2)
      sig.r[i, 3] = sum(bbr1$qval<alpha3)
      sig.r[i, 4] = nrow(bbr1)
      write.table(bbr1, fiL1c, row.names=F, col.names=T, quote=F, sep=",")

      bbr1 = bbr1[bbr1$nm %in% kpL,]
      qv = qvalue(bbr1$pval)
      bbr1$qval = qv$qvalue
      sig.r01[i, 1] = sum(bbr1$qval<alpha1)
      sig.r01[i, 2] = sum(bbr1$qval<alpha2)
      sig.r01[i, 3] = sum(bbr1$qval<alpha3)
      sig.r01[i, 4] = nrow(bbr1)
      write.table(bbr1, fiL1c01, row.names=F, col.names=T, quote=F, sep=",")

      bbr1 = bbr1[bbr1$nm %in% kpLb,]
      qv = qvalue(bbr1$pval)
      bbr1$qval = qv$qvalue
      sig.r001[i, 1] = sum(bbr1$qval<alpha1)
      sig.r001[i, 2] = sum(bbr1$qval<alpha2)
      sig.r001[i, 3] = sum(bbr1$qval<alpha3)
      sig.r001[i, 4] = nrow(bbr1)
      write.table(bbr1, fiL1c001, row.names=F, col.names=T, quote=F, sep=",")

      sig.r[i, ,drop=F]
      sig.r01[i, ,drop=F]
      sig.r001[i, ,drop=F]
    
      #middle model
      colnames(bbr2) = clnms2
      chk = which(bbr2[,24]==1 & is.na(abs(bbr2[,23]-bnr2[,21])))
      length(chk)
      bbr2[chk,24]=4
      chk = which(bbr2[,24]==1 & abs(bbr2[,23]-bnr2[,21])>=.01)
      length(chk)
      bbr2[chk,24]=3
      chk = which(bbr2[,24]==1 & abs(bbr2[,23]-bnr2[,21])<=.01)
      length(chk)
      bbr2[chk,3] = NA
      bbr2[chk,13] = NA
      ind0 = 1
      ind1b = (2:6); ind1 = ind1b[-2];ind1a = 2:5
      ind2b = (12:16); ind2 = ind2b[-2];ind2a = 11:14
      ind3 = 22;ind3a = 20
      ind4 = 23;ind4a = 21
      ind5 = 24
      bbr2[chk,ind1] = bnr2[chk,ind1a]
      bbr2[chk,ind2] = bnr2[chk,ind2a]
      bbr2[chk,ind3] = bnr2[chk,ind3a]
      bbr2[chk,ind4] = bnr2[chk,ind4a]
      table(bbr2[,24])
      bbr2 = bbr2[bbr2$status %in% 0:1,c(ind0, ind1b, ind2b, ind3, ind5)]
      qv = qvalue(bbr2$pval)
      bbr2$qval = qv$qvalue
      
      sig.r[i, 5] = sum(bbr2$qval<alpha1)
      sig.r[i, 6] = sum(bbr2$qval<alpha2)
      sig.r[i, 7] = sum(bbr2$qval<alpha3)
      sig.r[i, 8] = nrow(bbr2)
      write.table(bbr2, fiL2c, row.names=F, col.names=T, quote=F, sep=",")

      bbr2 = bbr2[bbr2$nm %in% kpL,]
      qv = qvalue(bbr2$pval)
      bbr2$qval = qv$qvalue
      sig.r01[i, 5] = sum(bbr2$qval<alpha1)
      sig.r01[i, 6] = sum(bbr2$qval<alpha2)
      sig.r01[i, 7] = sum(bbr2$qval<alpha3)
      sig.r01[i, 8] = nrow(bbr2)
      write.table(bbr2, fiL2c01, row.names=F, col.names=T, quote=F, sep=",")

      bbr2 = bbr2[bbr2$nm %in% kpLb,]
      qv = qvalue(bbr2$pval)
      bbr2$qval = qv$qvalue
      sig.r001[i, 5] = sum(bbr2$qval<alpha1)
      sig.r001[i, 6] = sum(bbr2$qval<alpha2)
      sig.r001[i, 7] = sum(bbr2$qval<alpha3)
      sig.r001[i, 8] = nrow(bbr2)
      write.table(bbr2, fiL2c001, row.names=F, col.names=T, quote=F, sep=",")

      sig.r[i, ,drop=F]
      sig.r01[i, ,drop=F]
      sig.r001[i, ,drop=F]
    
      
      
      #short model
      colnames(bbr3) = clnms2
      chk = which(bbr3[,24]==1 & is.na(abs(bbr3[,23]-bnr3[,21])))
      length(chk)
      bbr3[chk,24]=4
      chk = which(bbr3[,24]==1 & abs(bbr3[,23]-bnr3[,21])>=.01)
      length(chk)
      bbr3[chk,24]=3
      chk = which(bbr3[,24]==1 & abs(bbr3[,23]-bnr3[,21])<=.01)
      length(chk)
      bbr3[chk,3] = NA
      bbr3[chk,13] = NA
      ind1b = (2:4); ind1 = ind1b[-2];ind1a = 2:3
      ind2b = (12:14); ind2 = ind2b[-2];ind2a = 11:12
      ind3 = 22;ind3a = 20
      ind4 = 23;ind4a = 21
      ind5 = 24
      bbr3[chk,ind1] = bnr3[chk,ind1a]
      bbr3[chk,ind2] = bnr3[chk,ind2a]
      bbr3[chk,ind3] = bnr3[chk,ind3a]
      bbr3[chk,ind4] = bnr3[chk,ind4a]
      table(bbr3[,24])
      bbr3 = bbr3[bbr3$status %in% 0:1,c(ind0, ind1, ind2, ind3, ind5)]
      qv = qvalue(bbr3$pval)
      bbr3$qval = qv$qvalue
      bbr3[1:5,]
      sig.r[i, 9] = sum(bbr3$qval<alpha1)
      sig.r[i, 10] = sum(bbr3$qval<alpha2)
      sig.r[i, 11] = sum(bbr3$qval<alpha3)
      sig.r[i, 12] = nrow(bbr3)
      write.table(bbr3, fiL3c, row.names=F, col.names=T, quote=F, sep=",")

      bbr3 = bbr3[bbr3$nm %in% kpL,]
      qv = qvalue(bbr3$pval)
      bbr3$qval = qv$qvalue
      sig.r01[i, 9] = sum(bbr3$qval<alpha1)
      sig.r01[i, 10] = sum(bbr3$qval<alpha2)
      sig.r01[i, 11] = sum(bbr3$qval<alpha3)
      sig.r01[i, 12] = nrow(bbr3)
      write.table(bbr3, fiL3c01, row.names=F, col.names=T, quote=F, sep=",")

      bbr3 = bbr3[bbr3$nm %in% kpLb,]
      qv = qvalue(bbr3$pval)
      bbr3$qval = qv$qvalue
      sig.r001[i, 9] = sum(bbr3$qval<alpha1)
      sig.r001[i, 10] = sum(bbr3$qval<alpha2)
      sig.r001[i, 11] = sum(bbr3$qval<alpha3)
      sig.r001[i, 12] = nrow(bbr3)
      write.table(bbr3, fiL3c001, row.names=F, col.names=T, quote=F, sep=",")

      sig.r[i, ,drop=F]
      sig.r01[i, ,drop=F]
      sig.r001[i, ,drop=F]
    
      message(i)
    }
    colnames(sig.r) = colnames(sig.r01) = colnames(sig.r001) = c("L0.05", "L0.10", "L0.25", "Lng", 
    "M0.05", "M0.10", "M0.25", "Mng", "S0.05", "S0.10", "S0.25", "Sng")
    rownames(sig.r) = rownames(sig.r01) = rownames(sig.r001) = nams
    
    write.csv(sig.r, sprintf("%s/dynamic_%s_using_%s_eQTL_all.csv", rs.dir, conds[condi], mod), quote=F)
    write.csv(sig.r01, sprintf("%s/dynamic_%s_using_%s_eQTL_01.csv", rs.dir, conds[condi], mod), quote=F)
    write.csv(sig.r001, sprintf("%s/dynamic_%s_using_%s_eQTL_001.csv", rs.dir, conds[condi], mod), quote=F)
    
    set.seed(12345)
    labs = c(0,1, 5, 10,50,100,500)
    xns = rnorm(nrow(sig.r), 0, .02)
    yns = rnorm(nrow(sig.r), 0, .02)
    pdf(sprintf("%s/dynamic_%s_allEQTL_%s.pdf", rs.dir, conds[condi], mod), height=4, width=4)
    par(mfrow=c(1,1))
    x = log10(sig.r[,"L0.05"]+1)+xns
    y = log10(sig.r[,"S0.05"]+1)+yns
    plot(x, y, main="long vs short, q<0.05", bty="n", col=tiscol, pch=pchs, 
      xlab="#sig., log10 scale", ylab="#sig., log10 scale", xaxt="n", yaxt="n")
    abline(a=0, b=1, col="red")
    axis(1, at=log10(labs+1), labels=labs)
    axis(2, at=log10(labs+1), labels=labs)
    x = log10(sig.r[,"L0.10"]+1)+xns
    y = log10(sig.r[,"S0.10"]+1)+yns
    plot(x, y, main="long vs short, q<0.10", bty="n", col=tiscol, pch=pchs, 
      xlab="#sig., log10 scale", ylab="#sig., log10 scale", xaxt="n", yaxt="n")
    abline(a=0, b=1, col="red")
    axis(1, at=log10(labs+1), labels=labs)
    axis(2, at=log10(labs+1), labels=labs)
    x = log10(sig.r[,"L0.25"]+1)+xns
    y = log10(sig.r[,"S0.25"]+1)+yns
    plot(x, y, main="long vs short, q<0.25", bty="n", col=tiscol, pch=pchs, 
      xlab="#sig., log10 scale", ylab="#sig., log10 scale", xaxt="n", yaxt="n")
    abline(a=0, b=1, col="red")
    axis(1, at=log10(labs+1), labels=labs)
    axis(2, at=log10(labs+1), labels=labs)
    dev.off()
    
    #xns = rnorm(nrow(sig.r), 0, .015)
    #yns = rnorm(nrow(sig.r), 0, .015)
    pdf(sprintf("%s/dynamic_%s_sig01EQTL_%s.pdf", rs.dir, conds[condi], mod), height=4, width=4)
    par(mfrow=c(1,1))
    x = log10(sig.r01[,"L0.05"]+1)+xns
    y = log10(sig.r01[,"S0.05"]+1)+yns
    plot(x, y, main="long vs short, q<0.05", bty="n", col=tiscol, pch=pchs, 
      xlab="#sig., log10 scale", ylab="#sig., log10 scale", xaxt="n", yaxt="n")
    axis(1, at=log10(labs+1), labels=labs)
    axis(2, at=log10(labs+1), labels=labs)
    abline(a=0, b=1, col="red")
    x = log10(sig.r01[,"L0.10"]+1)+xns
    y = log10(sig.r01[,"S0.10"]+1)+yns
    plot(x, y, main="long vs short, q<0.10", bty="n", col=tiscol, pch=pchs, 
      xlab="#sig., log10 scale", ylab="#sig., log10 scale", xaxt="n", yaxt="n")
    axis(1, at=log10(labs+1), labels=labs)
    axis(2, at=log10(labs+1), labels=labs)
    abline(a=0, b=1, col="red")
    x = log10(sig.r01[,"L0.25"]+1)#+xns
    y = log10(sig.r01[,"S0.25"]+1)#+yns
    plot(x, y, main="long vs short, q<0.25", bty="n", col=tiscol, pch=pchs, 
      xlab="#sig., log10 scale", ylab="#sig., log10 scale", xaxt="n", yaxt="n")
    axis(1, at=log10(labs+1), labels=labs)
    axis(2, at=log10(labs+1), labels=labs)
    abline(a=0, b=1, col="red")
    dev.off()
    
    pdf(sprintf("%s/dynamic_%s_sig001EQTL_%s.pdf",rs.dir,  conds[condi], mod), height=4, width=4)
    par(mfrow=c(1,1))
    x = log10(sig.r001[,"L0.05"]+1)+xns
    y = log10(sig.r001[,"S0.05"]+1)+yns
    plot(x, y, main="long vs short, q<0.05", bty="n", col=tiscol, pch=pchs, 
      xlab="#sig., log10 scale", ylab="#sig., log10 scale", xaxt="n", yaxt="n")
    axis(1, at=log10(labs+1), labels=labs)
    axis(2, at=log10(labs+1), labels=labs)
    abline(a=0, b=1, col="red")
    x = log10(sig.r001[,"L0.10"]+1)+xns
    y = log10(sig.r001[,"S0.10"]+1)+yns
    plot(x, y, main="long vs short, q<0.10", bty="n", col=tiscol, pch=pchs, 
      xlab="#sig., log10 scale", ylab="#sig., log10 scale", xaxt="n", yaxt="n")
    axis(1, at=log10(labs+1), labels=labs)
    axis(2, at=log10(labs+1), labels=labs)
    abline(a=0, b=1, col="red")                                                                        
    x = log10(sig.r001[,"L0.25"]+1)#+xns
    y = log10(sig.r001[,"S0.25"]+1)#+yns
    plot(x, y, main="long vs short, q<0.25", bty="n", col=tiscol, pch=pchs, 
      xlab="#sig., log10 scale", ylab="#sig., log10 scale", xaxt="n", yaxt="n")
    axis(1, at=log10(labs+1), labels=labs)
    axis(2, at=log10(labs+1), labels=labs)
    abline(a=0, b=1, col="red")
    dev.off()
    }
}



#for(i in 1:length(dir0)){
#  pref = dir0[i]
#  bb.cl = sprintf("%s/%s/BetaBin", rt.dir, pref)
#  system(sprintf("rm -r %s/*", bb.cl))
#}



#get lists for pathways
outdirs = c("Adipose", "Artery", "Brain", "Esophagus", "Heart", "Skin", "Blood", "allTiss")
inds = list(1:2, 3:5, 6:10, 14:16, 17:18, 23:24, c(3:5,28), 1:28)
for(tisi in 1:length(outdirs)){
  #outdir = "allTiss"
  outdir = outdirs[tisi]
  ind = inds[[tisi]]
  
  if(!file.exists(outdir))dir.create(outdir)
  mod = modL
  for(asc in c("all", "001", "01")){
    for(condi in 1:3){
      for(mod in c(modL, modS)){
        lstS10 = lstS25 = 
        lstM10 = lstM25 = 
        lstL10 = lstL25 = 
        lstSB = lstMB = lstLB = character(0)
        for(i in ind){
          #i = 2
          pref = dir0[i]
          mod1 = mod
          bb.cl = sprintf("%s/%s/BetaBin", rt.dir, pref)
          suff = sprintf("%s_%s_%s", pref, conds[condi], mod)
          suf0 = sprintf("%s_%s", conds[condi], mod)
          fiL1c = sprintf("%s/%s_gPC_PF_%s.csv", bb.cl, suff, asc)
          fiL2c = sprintf("%s/%s_gPC_%s.csv", bb.cl, suff, asc)
          fiL3c = sprintf("%s/%s_nocov_%s.csv", bb.cl, suff, asc)
    
    
          #large model
          bbr1 = read.csv(fiL1c, as.is=T)
          lstL10 = c(lstL10, bbr1$nm[bbr1$qval<alpha2])
          lstL25 = c(lstL25, bbr1$nm[bbr1$qval<alpha3])
          lstLB = c(lstLB, bbr1$nm)
          
          #middle model
          bbr1 = read.csv(fiL2c, as.is=T)
          lstM10 = c(lstM10, bbr1$nm[bbr1$qval<alpha2])
          lstM25 = c(lstM25, bbr1$nm[bbr1$qval<alpha3])
          lstMB = c(lstMB, bbr1$nm)
    
          bbr1 = read.csv(fiL3c, as.is=T)
          lstS10 = c(lstS10, bbr1$nm[bbr1$qval<alpha2])
          lstS25 = c(lstS25, bbr1$nm[bbr1$qval<alpha3])
          lstSB = c(lstSB, bbr1$nm)
          #message(i)
        }
      
        out = sprintf("%s/L10_%s_%s.txt", outdir, suf0, asc)
        write.table(lstL10, out, row.names=F, col.names=F, quote=F)
        out = sprintf("%s/L25_%s_%s.txt", outdir, suf0, asc)
        write.table(lstL25, out, row.names=F, col.names=F, quote=F)
        out = sprintf("%s/LB_%s_%s.txt", outdir, suf0, asc)
        write.table(lstLB, out, row.names=F, col.names=F, quote=F)
        
        out = sprintf("%s/M10_%s_%s.txt", outdir, suf0, asc)
        write.table(lstM10, out, row.names=F, col.names=F, quote=F)
        out = sprintf("%s/M25_%s_%s.txt", outdir, suf0, asc)
        write.table(lstM25, out, row.names=F, col.names=F, quote=F)
        out = sprintf("%s/MB_%s_%s.txt", outdir, suf0, asc)
        write.table(lstMB, out, row.names=F, col.names=F, quote=F)
        
        out = sprintf("%s/S10_%s_%s.txt", outdir, suf0, asc)
        write.table(lstS10, out, row.names=F, col.names=F, quote=F)
        out = sprintf("%s/S25_%s_%s.txt", outdir, suf0, asc)
        write.table(lstS25, out, row.names=F, col.names=F, quote=F)
        out = sprintf("%s/SB_%s_%s.txt", outdir, suf0, asc)
        write.table(lstSB, out, row.names=F, col.names=F, quote=F)
      }
    }
  }
  message(tisi)
}

c(length(lstL10), length(lstL25), length(lstLB))
c(length(lstM10), length(lstM25), length(lstMB))
c(length(lstS10), length(lstS25), length(lstSB))

table(table(lstL10))
table(table(lstM10))
table(table(lstS10))

table(table(lstL25))
table(table(lstM25))
table(table(lstS25))

table(table(lstLB))
table(table(lstMB))
table(table(lstSB))


q("no")

#