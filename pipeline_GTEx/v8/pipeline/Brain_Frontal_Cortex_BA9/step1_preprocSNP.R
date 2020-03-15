args=(commandArgs(TRUE))
# args = c("21", "2e5")
# args = c("22", "5e5")
chri = as.numeric(args[1])
cis_window = as.numeric(args[2])
if(is.na(cis_window))cis_window=2e5
#cis_window = 2e5
chri 
cis_window

#module add samtools; module add tabix

permute = F
set.seed(1565691)
pref = "Brain_Frontal_Cortex_BA9"
modelS = "short"
modelL = "long"
library(Matrix)
source("helpers.R")

pars4gr = "/nas/longleaf/home/zhabotyn/progs/parser/parser4gr"
pars4grNA = "/nas/longleaf/home/zhabotyn/progs/parser/parser4grNA"
pars4grNAunf = "/nas/longleaf/home/zhabotyn/progs/parser/parser4grNAunf"
# ------------------------------------------------------------
# read in data
# ------------------------------------------------------------
root.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8"

inputdir = sprintf("%s/TReC_ASReC/%s", root.dir, pref)
preprdir = sprintf("%s/%s_prepr", root.dir, pref)
outputdir = sprintf("%s/%s_res", root.dir, pref)
if(!file.exists(outputdir))dir.create(outputdir)
if(!file.exists(preprdir))dir.create(preprdir)
cov.dir = sprintf("%s/Annotations/GTEx_Analysis_v8_eQTL_covariates", root.dir)
geno_dir = sprintf("%s/WGS_VCF", root.dir)
cnt.dir = sprintf("%s", inputdir)
outdir = sprintf("%s/trecase_res_example", inputdir)

rvcf.dir = sprintf("%s", geno_dir)
rvcf.sub = sprintf("%s/vcf_%s", geno_dir, pref)
if(!file.exists(rvcf.sub))dir.create(rvcf.sub)


  cvrtsh = read.csv(sprintf("%s/covariates_%s.csv", preprdir, modelS), as.is=T)
  cvrtlo = read.csv(sprintf("%s/covariates_%s.csv", preprdir, modelL), as.is=T)

  covariates_file_name = sprintf("%s/%s.v8.covariates.txt", cov.dir, pref)
  cvrt = read.table(covariates_file_name, as.is=T)
  rownames(cvrt) = cvrt[,1]; cvrt = cvrt[,-1]
  colnames(cvrt) = gsub("\\.", "-", colnames(cvrt))
  cvrt = data.frame(t(cvrt))

  #to.rm = c("GTEX-111YS", "GTEX-14C5O", "GTEX-148VJ", "GTEX-1497J", "GTEX-1AYCT", "GTEX-1S82P")
  to.rm = c("")
  dim(cvrt)
  rownames(cvrt) = cvrt[,1]
  cvrt = cvrt[,-1]
  kp.ind = !(rownames(cvrt) %in% to.rm)
  which(!kp.ind)
  cvrt = cvrt[kp.ind,]
  samcov = rownames(cvrt)
  dim(cvrt)
  
#chri = 22
  snpsfl = sprintf("%s/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_chr%s.vcf.gz", geno_dir, chri)

  tmpfl = gzfile(snpsfl)
  lns = readLines(tmpfl, n=3e4)
  close(tmpfl)
  lni = grep("CHROM", lns)
  samvcf = gsub("#", "", unlist(strsplit(lns[lni], split="\t")))

    m = match(samcov, samvcf)
    lns[lni] = paste(samvcf[m], collapse="\t")
    headervcf = sprintf("%s/header.vcf", rvcf.sub)
    if(!file.exists(headervcf))write.table(lns[c(1:3,4)], headervcf, row.names=F, col.names=F, quote=F)
    

vcfsam = read.table(sprintf("%s/sample.txt", rvcf.sub),as.is=T)
vcfout = sprintf("%s/SNP_chr%s.vcf", rvcf.sub, chri)
nms = vcfsam[,1]
table(nms[!(nms %in% to.rm)]==vcfsam[,1])



MAFcut = 0.05
#vcfout = sprintf("%s/GTEx_v8_WGS_chr%s.vcf.gz", rvcf.dir, chri)
comi = sprintf("zcat %s  | sed 's/ /\t/g'  | sed /^#/d  | cut  -f '10-' | %s | cut -f '1-3'", snpsfl, pars4grNA)
p = pipe(comi)
x4u = read.table(p, colClasses=c("integer"), fill=TRUE, row.names=NULL)
x4u[1:4,]
table(x4u[,3])
apply(x4u, 2, max)

m = match(samcov, samvcf[-(1:9)]);m
kp = x4u[,2] %in% m
table(kp)
x4u = x4u[kp,]
x4u[,2] = as.numeric(factor(x4u[,2], levels=m))
apply(x4u,2,max)
#table(x4u[,2])

geno = t(sparseMatrix(i=x4u[,2], j=x4u[,1], x=x4u[,3]))
dim(geno)


#comi = sprintf("zcat %s  | sed 's/ /\t/g' | sed /^#/d  | cut  -f '1-3'", "/pine/scr/z/h/zhabotyn/R01/GTEx/v8/WGS_VCF/GTEx_v8_WGS_chr22.vcf.gz")

#comi = sprintf("cat %s  | sed 's/ /\t/g' | sed /^#/d  | cut  -f '1-5'", SNP_file_name)
comi = sprintf("zcat %s  | sed 's/ /\t/g' | sed /^#/d  | cut  -f '1-3'", snpsfl)
p = pipe(comi)
pos = read.table(p, fill=TRUE, row.names=NULL)
pos[1:2,]
pos[,1] = gsub("chr", "", pos[,1])

SNPInfo = pos[1:nrow(geno), c(3,1,2)]
dim(SNPInfo)
dim(geno)
no.nas = apply(geno==-1, 1, sum)==0; table(no.nas)
geno = geno[no.nas,]
SNPInfo = SNPInfo[no.nas,]

#MAFu = apply(geno, 1, function(x) (sum(x %in% c(1,2)) + 2*sum(x==3)) /(2*length(x)))
#summary(MAFu)
#table(MAFu>=MAFcut & MAFu<=(1-MAFcut))


#load back counts
  ginfo = sprintf("%s/geneInfo_prepr_%s.txt", preprdir, modelS); ginfo
  geneInfo = read.table(ginfo, as.is=T)
  geneInfo[,2] = gsub("chr", "", geneInfo[,2])
  
  ase1f = sprintf("%s/ase1_prepr_%s.txt", preprdir, modelS); ase1f
  ase1 = read.table(ase1f, as.is=T)

#  ase1f = sprintf("%s/ase1_prepr_%s.txt", preprdir, modelL); ase1f
#  ase1l = read.table(ase1f, as.is=T)

  ase2f = sprintf("%s/ase2_prepr_%s.txt", preprdir, modelS); ase2f
  ase2 = read.table(ase2f, as.is=T)

  totp = sprintf("%s/Tcnt_prepr_%s.txt", preprdir, modelS); totp
  trecDs = read.table(totp, as.is=T)

  totp = sprintf("%s/Tcnt_prepr_%s.txt", preprdir, modelL); totp
  trecDl = read.table(totp, as.is=T)
  table(unlist(trecDl==trecDs))
  colnames(trecDl) = gsub("\\.", "-", colnames(trecDl))
  colnames(trecDs) = gsub("\\.", "-", colnames(trecDs))
  #add normalized total expr
  depthL = 10^cvrtlo[,ncol(cvrtlo)]
  depthL = depthL/median(depthL)
  exprL = trecDl
  exprL = log1p(as.matrix(exprL)%*%diag(depthL))
  exprL = t(apply(exprL, 1, normscore))

  depthS = 10^cvrtsh[,ncol(cvrtsh)]
  depthS = depthL/median(depthS)
  exprS = trecDs
  exprS = log1p(as.matrix(exprS)%*%diag(depthS))
  exprS = t(apply(exprS, 1, normscore))

  normS = sprintf("%s/Tcnt_norm_%s.txt", preprdir, modelS)
  write.table(exprS, normS, row.names=F, col.names=F, quote=F)
  normL = sprintf("%s/Tcnt_norm_%s.txt", preprdir, modelL)
  write.table(exprL, normL, row.names=F, col.names=F, quote=F)


  kpchr = geneInfo[,2]==chri;table(kpchr)
  trecDl = trecDl[kpchr,]
  trecDs = trecDs[kpchr,]
  ase1 = ase1[kpchr,]
  ase2 = ase2[kpchr,]
  geneInfo = geneInfo[kpchr,]

  exprL = exprL[kpchr,]
  exprS = exprS[kpchr,]

  summary(ase1[,1:5])
  summary(ase2[,1:5])
  summary(colSums(ase1))
  summary(colSums(ase2))
  summary(colSums(ase1)/colSums(ase1+ase2))

numsam = ncol(geno)  
nsub = ncol(geno)
subsA = sample(1:nsub)


for(nsub in numsam){

  subs = subsA[1:nsub]
  subse = c(0, subs)+1
  subsf = c(1:9, subs+9)
  subs
  subsf
  int.dir = sprintf("%s_%s_%s", pref, nsub, cis_window)
  if(!file.exists(int.dir))dir.create(int.dir)

  sam.ord = cvrtlo[subs,1]
  sam.ord2 = colnames(trecDs)[subs]
  table(sam.ord==sam.ord2)
  
  sam.ordf = sprintf("%s/samples.dat", int.dir)
  sam.ordf
  if(!file.exists(sam.ordf))write.table(sam.ord, sam.ordf, row.names=F, col.names=F, quote=F)
  #MAF = apply(geno, 1, function(x) (sum(x %in% c(1,2)) + 2*sum(x==3)) /(2*length(x)))
  #summary(MAF)
  
  # filter out maf < 0.05
  genos = geno[,subs] #matrix(, ncol=nsub)
  MAF = apply(genos, 1, function(x) (sum(x %in% c(1,2)) + 2*sum(x==3)) /(2*length(x)))
  summary(MAF)
  
  geno2kp = which(MAF >= MAFcut & MAF  <= (1-MAFcut))
  length(geno2kp)
  genos = genos[geno2kp, ] #matrix(, ncol=nsub)
  dim(genos)
  
  SNPInfo = SNPInfo[geno2kp, ]
  #colnames(SNPInfo) = c("snp","chr",'pos')
  colnames(SNPInfo) = c("snpid", "chr", "pos")
  dim(SNPInfo)
  head(SNPInfo)
  


  # ------------------------------------------------------------
  # asSeq2
  # ------------------------------------------------------------
  #setwd(outdir0)
  exprL = exprL[,subs]
  exprS = exprS[,subs]
  
  trecDsi = t(trecDs[,subs])
  trecDsi[1:5,1:5]
  
  trecDli = t(trecDl[,subs])
  trecDsi[1:5,1:5]
  
  ase1i = t(ase1[,subs])
  ase1i[1:5,1:5]
  
  #ase1si = 
  
  ase2i = t(ase2[,subs])
  ase2i[1:5,1:5]
  
  
  
  #genos = t(genos)
  
  #X2 = model.matrix(~., data = cvrt[subs,-c(5,10)])
  tmp = model.frame(cvrtsh[subs,-1])
  #f = paste("~", paste(colnames(cvrtsh)[-1],collapse="+"), collapse="")
  #X2sh = model.matrix(f, data = tmp)
  X2sh = model.matrix(~., data = tmp)
  dim(X2sh)
  Xsh = X2sh[,-1]
  dim(Xsh)
  
  tmp = model.frame(cvrtlo[subs,-1])
  X2lo = model.matrix(~., data = tmp)
  dim(X2lo)
  Xlo = X2lo[,-1]
  dim(Xlo)
  
  nr = nrow(genos);nr
  time1 = proc.time()
  genos[genos == 3] = 4
  genos[genos == 2] = 3 
  genos[1:4,1:14]
  genos = matrix(genos, nrow=nr) 
  genos[1:4,1:14]
  
  eChr = rep(chri, nrow(geneInfo))
  ePos = as.numeric((geneInfo[,3] + geneInfo[,4])/2)
  
  mChr = rep(chri, nr)
  geni = 1
  numsnps = rep(0, nrow(geneInfo));nrow(geneInfo)
  for(geni in 1:nrow(geneInfo)){
  #geni = 1
  #message("geno: ", nrow(geno), " ", ncol(geno), " trecD: ", nrow(trecD), " ", ncol(trecD), " ", nrow(geneInfo))
  #geni = 389
  
  #geni = indi
  
    local.distances = as.numeric((geneInfo[,4] - geneInfo[,3]))/2+cis_window
    kp = which(SNPInfo[,3]>=(ePos[geni]-local.distances[geni]) & SNPInfo[,3]<(ePos[geni]+local.distances[geni]));length(kp)
    numsnps[geni] = length(kp)
    if(length(kp)>0){
      genfi = genos[kp,,drop=F]
      SNPInfoi = SNPInfo[kp,,drop=F]
      dim(genfi)
      dim(SNPInfoi)
      suffi = sprintf("%s_%s", chri, geni)
      write.table(genfi, sprintf("%s/genotypes_%s.dat", int.dir, suffi), row.names=F, col.names=T, quote=F, sep="\t")
      write.table(SNPInfoi, sprintf("%s/genotypei_%s.dat", int.dir, suffi), row.names=F, col.names=T, quote=F, sep="\t")

      expression_file_nameS = sprintf("%s/GE_norm_%s_%s.dat", int.dir, modelS, suffi)
      exprSj = t(matrix(exprS[geni,]))
      write.table(exprSj, expression_file_nameS, row.names=F, col.names=F, quote=F, sep="\t")

      expression_file_nameL = sprintf("%s/GE_norm_%s_%s.dat", int.dir, modelL, suffi)
      exprLj = t(matrix(exprL[geni,]))
      write.table(exprLj, expression_file_nameL, row.names=F, col.names=F, quote=F, sep="\t")

      write.table(cbind(t(t(trecDsi[,geni])), t(t(ase1i[,geni])), t(t(ase2i[,geni]))),
        sprintf("%s/counti_%s_%s.csv", int.dir, modelS, suffi), row.names=F, col.names=F, quote=F, sep=",")
  
      write.table(cbind(t(t(trecDli[,geni])), t(t(ase1i[,geni])), t(t(ase2i[,geni]))),
        sprintf("%s/counti_%s_%s.csv", int.dir, modelL, suffi), row.names=F, col.names=F, quote=F, sep=",")
  
      Xmatfils = sprintf("%s/Xmat_%s.csv", int.dir, modelS)
      if(!file.exists(Xmatfils))write.table(Xsh, Xmatfils, row.names=F, col.names=F, quote=F, sep=",")
  
      Xmatfillo = sprintf("%s/Xmat_%s.csv", int.dir, modelL)
      if(!file.exists(Xmatfillo))write.table(Xlo, Xmatfillo, row.names=F, col.names=F, quote=F, sep=",")
    }
    message(geni, " out of ", nrow(geneInfo), " #snps: ", length(kp))
  }
}
summary(numsnps)

q(save = 'no')


