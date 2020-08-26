#--begin=21:00:00 -t 03-00:00:00 
#sbatch -p general -o step4_collectU.out --mem=2g --wrap="R CMD BATCH step4_collectU.R step4_collectU.Rout"
args=(commandArgs(TRUE))
#args = c("5e5")

cis_window = as.numeric(args[1])
numpoints = 100

model = "long"
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
out.dir = sprintf("res%s", cis_window)
if(!file.exists(out.dir))dir.create(out.dir)

nsub = as.numeric(nsam)
#root.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8"
preprdir = sprintf("%s_prepr", pref)
#geno_dir = sprintf("%s/WGS_VCF", root.dir)
int.dir = sprintf("%s_%s_%s", pref, nsub, cis_window)

boot.dir = sprintf("boot_%s_%s_%s_%s_%s", pref, nsub, cis_window, model, numpoints)

#hg38 = read.csv("/pine/scr/z/h/zhabotyn/R01/info/hg38/hg38.csv")

info = read.table(sprintf("%s/geneInfo_prepr_%s.txt", preprdir, model), header=T)
dim(info)

bfls = list.files(boot.dir, pattern="upd_eigenMT");length(bfls)

#Xl = read.csv(sprintf("%s/Xmat_long.csv", int.dir), header=F)
#Xs = read.csv(sprintf("%s/Xmat_short.csv", int.dir), header=F)
samples = read.table(sprintf("%s/samples.dat", int.dir))
cntsM = matrix(NA, nrow=length(bfls), ncol=nsub)# = cntsA = cntsT
rownames(cntsM) = 1:nrow(cntsM)
colnames(cntsM) = samples[,1]
cntsRJ = cntsRT = cntsRA = cntsM = list(cntsM, cntsM)

assJ = assT = assA = data.frame(matrix(NA, nrow=length(bfls), ncol=24))
tokp = rep(FALSE, nrow(assJ));tokp = list(tokp, tokp)

assJm = assTm = assAm = list(assJ, assJ)
emtres = data.frame(matrix(NA, nrow=length(bfls), ncol=21))
bootres = data.frame(matrix(NA, nrow=length(bfls), ncol=10))
colnames(bootres) = c("id", "snpid", "minp", sprintf("lm%s", 1:3), sprintf("glm%s", 1:4))
emtres = list(emtres, emtres)
bootres = list(bootres, bootres)

logiti = function(x){1/(1+exp(-x))}
model = c("short", "long")

for(i in 1:length(bfls)){
#for(i in c(1,13100:13200)){
#for(i in 1:100){
#message(i)
for(j in 1:length(model)){
#i = i + 1
  outdir0 = sprintf("run_%s_%s_%s_%s", pref, nsub, cis_window, model[[j]])
  modelj = model[j]
  boot.dir = sprintf("boot_%s_%s_%s_%s_%s", pref, nsub, cis_window, modelj, numpoints)

  fli = sprintf("%s/%s", boot.dir, bfls[i])
  if(file.size(fli)>0){
    kp = T
  }
  if(kp){
    resi = read.csv(fli, as.is=T)
    bootres[[j]][i,1:3] = resi[1,c("gene", "SNP", "p.value")]
    if(i == 1)colnames(emtres[[j]]) = colnames(resi)
    emtres[[j]][i,] = resi[1,]

    xmin = -log10(resi$p.value)
  }
  
  suffi = gsub(".csv", "", gsub("upd_eigenMT_","", bfls[i]))
  genfi = read.table(sprintf("%s/genotypes_%s.dat", int.dir, suffi), as.is=T, header=T)
  infi = read.table(sprintf("%s/genotypei_%s.dat", int.dir, suffi), as.is=T, header=T)
  #tot, hap1, hap2
  cntT = cntM = read.csv(sprintf("%s/counti_%s_%s.csv", int.dir, model[[j]], suffi), as.is=T, header=F)

  #timi = read.table(sprintf("%s/%s_time.txt", outdir0, resi$gene), sep="\t", header=T, as.is=T)
  #assi0 = read.table(sprintf("%s/%s_eqtl.txt", outdir0a, resi$gene), sep="\t", header=T, as.is=T)
  #timi0 = read.table(sprintf("%s/%s_time.txt", outdir0a, resi$gene), sep="\t", header=T, as.is=T)
  assi = read.table(sprintf("%s/%s_eqtl.txt", outdir0, resi$gene), sep="\t", header=T, as.is=T)
  if(nrow(assi)>0){
  tokp[[j]][i] = TRUE
  if(i == 1)colnames(assJm[[j]])=colnames(assTm[[j]])=colnames(assAm[[j]]) = colnames(assi)
  oJ = order(assi$final_Pvalue)
  assJm[[j]][i,] = assi[oJ[1],]
  oT = order(assi$TReC_Pvalue)
  assTm[[j]][i,] = assi[oT[1],]
  oA = order(assi$ASE_Pvalue)
  assAm[[j]][i,] = assi[oA[1],] 
  #assi$ASE_Pvalue[oA[1]]
  
  m = match(resi$SNP, infi$snpid)
  cntM$snpM = unlist(genfi[m,])
  m = match((assJm[[j]])$MarkerRowID[i], infi$snpid);m
  cntM$snpJ = unlist(genfi[m,])
  m = match((assTm[[j]])$MarkerRowID[i], infi$snpid);m
  cntM$snpT = unlist(genfi[m,])
  m = match((assAm[[j]])$MarkerRowID[i], infi$snpid);m
  cntM$snpA = unlist(genfi[m,])
  #cntM
  #c(assJm[[j]]$final_Pvalue[i], assTm[[j]]$TReC_Pvalue[i], assAm[[j]]$ASE_Pvalue[i])
  cntsM[[j]][i,] = apply(cntM, 1, paste, collapse="|")
  rownames(cntsM[[j]])[i] = resi$gene
  #cntM order:
  #total, hap1, hap2, MatrixEQTL, Final, Trec, ASE
  #recode cntsM to format used for fitting beta-binomial submodel for dynamic eQTL

  #taking the most significant SNP classified into 4 groups according to TReCASE model AA, AB, BA and BB
  #if SNP was from the group 0, 1 or 3 - keep the values as they are
  #for the group 2 - swap hap1 and hap2 values.
  #Call these processed values as A1 and A2.
  #Also we classified whether a SNP for a given individual and gene is heterozygous: 0 - homozygous, 1 - heterozygous. Call these SNPInd
  #The output file is written in a format: A1|A2|SNPInd.
  cntM0 = cntM
  h1 = cntM0[cntM0$snpJ==2, 2]
  cntM0[cntM0$snpJ==2, 2] = cntM0[cntM0$snpJ==2, 3]
  cntM0[cntM0$snpJ==2, 3] = h1
  cntM0[cntM$snpJ %in% c(0,3), 4] = 0
  cntM0[cntM$snpJ %in% c(1,2), 4] = 1
  cntsRJ[[j]][i,] = apply(cntM0[,c(2:4)], 1, paste, collapse="|")
  rownames(cntsRJ[[j]])[i] = resi$gene

  cntM0 = cntM
  h1 = cntM0[cntM0$snpT==2, 2]
  cntM0[cntM0$snpT==2, 2] = cntM0[cntM0$snpT==2, 3]
  cntM0[cntM0$snpT==2, 3] = h1
  cntM0[cntM$snpT %in% c(0,3), 4] = 0
  cntM0[cntM$snpT %in% c(1,2), 4] = 1
  cntsRT[[j]][i,] = apply(cntM0[,c(2:4)], 1, paste, collapse="|")
  rownames(cntsRT[[j]])[i] = resi$gene

  cntM0 = cntM
  h1 = cntM0[cntM0$snpA==2, 2]
  cntM0[cntM0$snpA==2, 2] = cntM0[cntM0$snpA==2, 3]
  cntM0[cntM0$snpA==2, 3] = h1
  cntM0[cntM$snpA %in% c(0,3), 4] = 0
  cntM0[cntM$snpA %in% c(1,2), 4] = 1
  cntsRA[[j]][i,] = apply(cntM0[,c(2:4)], 1, paste, collapse="|")
  rownames(cntsRA[[j]])[i] = resi$gene
#}}


  fli = sprintf("%s/%s", boot.dir, gsub("upd_eigenMT", "short_boot_pval", bfls[i]))
  if(kp)kp = file.exists(fli)
  if(kp)kp = file.size(fli)>0
  kp
  if(kp){
    bresi = read.csv(fli, as.is=T)
    x = -log10(bresi[,2])
    y = bresi[,3]
    kp1 = y>=0.002 & y<=0.2
    kp2 = y>=0.002 & y<=0.3
    kp3 = y>=0     & y<=0.3
    kp3a = y>0     & y<=0.3
    y2 = y*1000
    n2 = 1000-y2
    y = -log10(y)
    lm1 = lm(y[kp1]~x[kp1])
    lm2 = lm(y[kp2]~x[kp2])
    lm3 = lm(y[kp3a]~x[kp3a])
    10^-c(lm1$coef[1]+lm1$coef[2]*xmin, lm2$coef[1]+lm2$coef[2]*xmin, lm3$coef[1]+lm3$coef[2]*xmin)
    glm1 = glm(cbind(y2, n2)[kp1,]~x[kp1], family="binomial")
    glm2 = glm(cbind(y2, n2)[kp2,]~x[kp2], family="binomial")
    glm3 = glm(cbind(y2, n2)[kp3,]~x[kp3], family="binomial")
    glm4 = glm(cbind(y2, n2)~x, family="binomial")
    lmr = 10^-c(lm1$coef[1]+lm1$coef[2]*xmin, lm2$coef[1]+lm2$coef[2]*xmin, lm3$coef[1]+lm3$coef[2]*xmin)
    glmr = logiti(c(glm1$coef[1]+glm1$coef[2]*xmin, glm2$coef[1]+glm2$coef[2]*xmin, glm3$coef[1]+glm3$coef[2]*xmin, glm4$coef[1]+glm4$coef[2]*xmin))
    bootres[[j]][i, 4:10] = c(lmr, glmr)
  }
  }else{
    tokp[[j]][i] = FALSE
    message(suffi)
  }
  rownames(assJm[[j]])[i]=rownames(assTm[[j]])[i]=rownames(assAm[[j]])[i] = resi$gene
  }
  if(i %% 100 == 0) message(i, " out of ", length(bfls))
}

#based on our past checks the settings used for glm, version 3 (kp3) will be used
#y>=0 & y<=0.3
setwd(out.dir)
j = 1
for(j in 1:length(model)){
  bootres[[j]] = bootres[[j]][tokp[[j]],]
  assJm[[j]] = assJm[[j]][tokp[[j]],]
  assTm[[j]] = assTm[[j]][tokp[[j]],]
  assAm[[j]] = assAm[[j]][tokp[[j]],]
  emtres[[j]] = emtres[[j]][tokp[[j]],]

  cntsM[[j]] = cntsM[[j]][tokp[[j]],]
  cntsRJ[[j]] = cntsRJ[[j]][tokp[[j]],]
  cntsRT[[j]] = cntsRT[[j]][tokp[[j]],]
  cntsRA[[j]] = cntsRA[[j]][tokp[[j]],] 
  
  bootres[[j]]$numtest = bootres[[j]]$glm3/bootres[[j]]$minp
  #dont get less then 1 effective SNP
  bootres[[j]]$numtest[which(bootres[[j]]$numtest<1)] = 1
  kp = which(bootres[[j]]$numtest>emtres[[j]]$TESTS)
  #don't exceed number of SNPs
  bootres[[j]]$numtest[kp] = emtres[[j]]$TESTS[kp]
  
  #calculate joint estimated permutation p-value
  assJm[[j]]$permp = assJm[[j]]$final_Pvalue*bootres[[j]]$numtest
  kp = which(is.na(assJm[[j]]$permp))
  assJm[[j]]$permp[kp] = assJm[[j]]$ASE_Pvalue[kp]*bootres[[j]]$numtest[kp]
  
  #calculate TReC estimated permutation p-value
  assTm[[j]]$permp = assTm[[j]]$TReC_Pvalue*bootres[[j]]$numtest
  #calculate ASE estimated permutation p-value
  assAm[[j]]$permp = assAm[[j]]$ASE_Pvalue*bootres[[j]]$numtest
  
  #don't exceed 1
  assJm[[j]]$permp[assJm[[j]]$permp>1] = 1
  assAm[[j]]$permp[assAm[[j]]$permp>1] = 1
  assTm[[j]]$permp[assTm[[j]]$permp>1] = 1
  
  modelj = model[[j]]
  suff = sprintf("%s_%s_%s_%s", pref, nsub, cis_window, modelj)
  write.csv(assJm[[j]],  sprintf("TReCASE_%s.csv", suff))
  write.csv(assTm[[j]],  sprintf("TReCASE_TReC_%s.csv", suff))
  write.csv(assAm[[j]],  sprintf("TReCASE_ASE_%s.csv", suff))
  write.csv(bootres[[j]],  sprintf("permp_est_%s.csv", suff))
  write.csv(emtres[[j]],  sprintf("eigenMT_est_%s.csv", suff))
  
  #counts along with minimum SNP for each of the models
  write.csv(cntsM[[j]], sprintf("counts_with_min_snp_%s.csv", suff))
  #add converted counts for joint model minimum SNP so that 
  write.csv(cntsRJ[[j]], sprintf("rec_counts_for_TReCASE_%s.csv", suff))
  write.csv(cntsRT[[j]], sprintf("rec_counts_for_TReC_%s.csv", suff))
  write.csv(cntsRA[[j]], sprintf("rec_counts_for_ASE_%s.csv", suff))

  #keep shorter version
  #available columns
  #Pos	MarkerRowID	NBod	BBod	NBod0	BBod0	
  #TReC_b	TReC_Chisq	TReC_df	TReC_Pvalue	
  #ASE_b	ASE_Chisq	ASE_df	ASE_Pvalue	
  #Joint_b	Joint_Chisq	Joint_df	Joint_Pvalue	
  #n_TReC	n_ASE	n_ASE_Het	trans_Chisq	trans_Pvalue	final_Pvalue
  outT = assTm[[j]][,c("MarkerRowID", "Pos", "NBod", "TReC_b", "TReC_Pvalue", "permp")]
  outA = assAm[[j]][,c("MarkerRowID", "Pos", "BBod", "ASE_b", "ASE_Pvalue", "permp")]
  outJ = assJm[[j]][,c("MarkerRowID", "Pos", "NBod", "BBod", "TReC_b", "TReC_Pvalue", 
  "ASE_b", "ASE_Pvalue", "Joint_b", "Joint_Pvalue", "trans_Pvalue", "final_Pvalue", "permp")]
  
  #removing some of the technical information
  write.csv(outJ[!is.na(outJ$permp),], sprintf("%s_TReCASE_trimmed.csv", suff), quote=F)
  write.csv(outT[!is.na(outT$permp),], sprintf("%s_TReC_trimmed.csv", suff), quote=F)
  write.csv(outA[!is.na(outA$permp),], sprintf("%s_ASE_trimmed.csv", suff), quote=F)


  alpha = 0.05
  message(sum(tokp[[j]]), " ",
  mean(emtres[[j]]$BF<alpha, na.rm=T), " ",
  mean(bootres[[j]]$glm3<alpha, na.rm=T), " ",
  mean(assJm[[j]]$permp<alpha, na.rm=T), " ",
  mean(assTm[[j]]$permp<alpha, na.rm=T), " ",
  mean(assAm[[j]]$permp<alpha, na.rm=T))
}

q("no")
