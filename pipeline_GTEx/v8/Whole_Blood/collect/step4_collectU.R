#--begin=21:00:00 -t 03-00:00:00 
#sbatch -p general -o step4_collectU.out -t 07-00:00:00 --mem=8g --wrap="R CMD BATCH step4_collectU.R step4_collectU.Rout"

nsub = 670
cis_window = 5e5
numpoints = 100

model = "long"
pref = "Whole_Blood"

root.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8"
inputdir = sprintf("%s/TReC_ASReC/%s", root.dir, pref)
preprdir = sprintf("%s/%s_prepr", root.dir, pref)
outputdir = sprintf("%s/%s_res", root.dir, pref)
cov.dir = sprintf("%s/Annotations/GTEx_Analysis_v8_eQTL_covariates", root.dir)
geno_dir = sprintf("%s/WGS_VCF", root.dir)
cnt.dir = sprintf("%s", inputdir)
outdir = sprintf("%s/trecase_res_example", inputdir)

rvcf.dir = sprintf("%s", geno_dir)
rvcf.sub = sprintf("%s/vcf_%s", geno_dir, pref)

int.dir = sprintf("%s_%s_%s", pref, nsub, cis_window)


boot.dir = sprintf("boot_%s_%s_%s_%s_%s", pref, nsub, cis_window, model, numpoints)

#hg38 = read.csv("/pine/scr/z/h/zhabotyn/R01/info/hg38/hg38.csv")

info = read.table(sprintf("%s/geneInfo_prepr_%s.txt", preprdir, model), header=T)
dim(info)

bfls = list.files(boot.dir, pattern="upd_eigenMT");length(bfls)

#Xl = read.csv(sprintf("%s/Xmat_long.csv", int.dir), header=F)
#Xs = read.csv(sprintf("%s/Xmat_short.csv", int.dir), header=F)
samples = read.table(sprintf("%s/samples.dat", int.dir))
cntsM = cntsA = cntsT = matrix(NA, nrow=length(bfls), ncol=nsub)
rownames(cntsM) = 1:nrow(cntsM)
colnames(cntsM) = samples[,1]
cntsM = list(cntsM, cntsM)

assJ = assT = assA = data.frame(matrix(NA, nrow=length(bfls), ncol=24))
assJm = assTm = assAm = list(assJ, assJ)
emtres = data.frame(matrix(NA, nrow=length(bfls), ncol=21))
bootres = data.frame(matrix(NA, nrow=length(bfls), ncol=10))
colnames(bootres) = c("id", "snpid", "minp", sprintf("lm%s", 1:3), sprintf("glm%s", 1:4))
emtres = list(emtres, emtres)
bootres = list(bootres, bootres)

logiti = function(x){1/(1+exp(-x))}
model = c("short", "long")

logiti = function(x){1/(1+exp(-x))}
for(i in 1:length(bfls)){
#for(i in 1:100){
  for(j in 1:length(model)){
#i = i + 1
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
  
  outdir0 = sprintf("run_%s_%s_%s_%s", pref, nsub, cis_window, model[[j]])
  suffi = gsub(".csv", "", gsub("upd_eigenMT_","", bfls[i]))
  genfi = read.table(sprintf("%s/genotypes_%s.dat", int.dir, suffi), as.is=T, header=T)
  infi = read.table(sprintf("%s/genotypei_%s.dat", int.dir, suffi), as.is=T, header=T)
  cntT = cntM = read.csv(sprintf("%s/counti_%s_%s.csv", int.dir, model[[j]], suffi), as.is=T, header=F)

  assi = read.table(sprintf("%s/%s_eqtl.txt", outdir0, resi$gene), sep="\t", header=T, as.is=T)
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

  m = match(resi$SNP, infi$snpid)
  cntM$snpM = unlist(genfi[m,])
  m = match((assJm[[j]])$MarkerRowID[i], infi$snpid);m
  cntM$snpJ = unlist(genfi[m,])
  m = match((assTm[[j]])$MarkerRowID[i], infi$snpid);m
  cntM$snpT = unlist(genfi[m,])
  m = match((assAm[[j]])$MarkerRowID[i], infi$snpid);m
  cntM$snpA = unlist(genfi[m,])
  #cntM
  c(assJm[[j]]$final_Pvalue[i], assTm[[j]]$TReC_Pvalue[i], assAm[[j]]$ASE_Pvalue[i])
  cntsM[[j]][i,] = apply(cntM, 1, paste, collapse=":")
  rownames(cntsM[[j]])[i] = resi$gene
  
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
  
  }
  if(i %% 1000 == 0) message(i, " out of ", length(bfls))
}

j = 1
for(j in 1:2){
  bootres[[j]]$numtest = bootres[[j]]$glm3/bootres[[j]]$minp
  bootres[[j]]$numtest[which(bootres[[j]]$numtest<1)] = 1
  kp = which(bootres[[j]]$numtest>emtres[[j]]$TESTS)
  bootres[[j]]$numtest[kp] = emtres[[j]]$TESTS[kp]
  
  assJm[[j]]$permp = assJm[[j]]$final_Pvalue*bootres[[j]]$numtest
  kp = which(is.na(assJm[[j]]$permp))
  assJm[[j]]$permp[kp] = assJm[[j]]$ASE_Pvalue[kp]*bootres[[j]]$numtest[kp]
  
  assTm[[j]]$permp = assTm[[j]]$TReC_Pvalue*bootres[[j]]$numtest
  assAm[[j]]$permp = assAm[[j]]$ASE_Pvalue*bootres[[j]]$numtest
  
  modelj = model[[j]]
  write.csv(assJm[[j]],  sprintf("TReCASE_%s_%s_%s_%s.csv", pref, nsub, cis_window, modelj))
  write.csv(assTm[[j]],  sprintf("TReCASE_TReC_%s_%s_%s_%s.csv", pref, nsub, cis_window, modelj))
  write.csv(assAm[[j]],  sprintf("TReCASE_ASE_%s_%s_%s_%s.csv", pref, nsub, cis_window, modelj))
  write.csv(bootres[[j]],  sprintf("permp_est_%s_%s_%s_%s.csv", pref, nsub, cis_window, modelj))
  write.csv(emtres[[j]],  sprintf("eigenMT_est_%s_%s_%s_%s.csv", pref, nsub, cis_window, modelj))
  
  write.csv(cntsM[[j]], sprintf("counts_with_min_snp_%s_%s_%s_%s.csv", pref, nsub, cis_window, modelj))
  alpha = 0.05
  message(mean(emtres[[j]]$BF<alpha, na.rm=T), " ",
  mean(bootres[[j]]$glm3<alpha, na.rm=T), " ",
  mean(assJm[[j]]$permp<alpha, na.rm=T), " ",
  mean(assTm[[j]]$permp<alpha, na.rm=T), " ",
  mean(assAm[[j]]$permp<alpha, na.rm=T))
}
