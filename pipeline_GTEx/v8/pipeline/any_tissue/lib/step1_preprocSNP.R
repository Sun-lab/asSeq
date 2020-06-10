args=(commandArgs(TRUE))
# args = c("22", "5e5")
chri = as.numeric(args[1])
cis_window = as.numeric(args[2])
if(is.na(cis_window))cis_window=2e5
chri 
cis_window
cis_window_char = as.character(cis_window)

options(scipen=999)
#module add samtools; module add tabix

permute = F
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
#"GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_chr"
vcf.pre = specs[17]
setwd(wrk.dir)

set.seed(seedval)

pars4grNAunf = specs[8]
numsam = nsub = nsam

modelS = "short"
modelL = "long"
library(Matrix)
source("../helpers.R")

# ------------------------------------------------------------
# read in data
# ------------------------------------------------------------

inputdir = sprintf("%s/TReC_ASReC/%s", bas.dir, pref)
preprdir = sprintf("%s_prepr", pref)
outputdir = sprintf("%s_res", pref)
if(!file.exists(outputdir))dir.create(outputdir)
if(!file.exists(preprdir))dir.create(preprdir)
cov.dir = sprintf("%s/Annotations/GTEx_Analysis_v8_eQTL_covariates", bas.dir)
geno_dir = sprintf("%s/WGS_VCF", bas.dir)
cnt.dir = sprintf("%s", inputdir)

rvcf.dir = sprintf("%s", geno_dir)
rvcf.sub = sprintf("%s/vcf_%s", geno_dir, pref)
if(!file.exists(rvcf.sub))dir.create(rvcf.sub)


cvrtsh = read.csv(sprintf("%s/covariates_%s.csv", preprdir, modelS), as.is=T)
cvrtlo = read.csv(sprintf("%s/covariates_%s.csv", preprdir, modelL), as.is=T)

covariates_file_name = sprintf("%s/%s.v8.covariates.txt", cov.dir, pref)
cvrt = read.table(covariates_file_name, as.is=T, header=T)
rownames(cvrt) = cvrt[,1]; cvrt = cvrt[,-1]
colnames(cvrt) = gsub("\\.", "-", colnames(cvrt))
cvrt = data.frame(t(cvrt))

#to.rm = c("GTEX-111YS", "GTEX-14C5O", "GTEX-148VJ", "GTEX-1497J", "GTEX-1AYCT", "GTEX-1S82P")
to.rm = rownames(cvrt)[is.na(match(rownames(cvrt),cvrtlo[,1]))]
to.rm
dim(cvrt)

kp.ind = !(rownames(cvrt) %in% to.rm)
which(!kp.ind)
cvrt = cvrt[kp.ind,]
samcov = rownames(cvrt)
dim(cvrt)
#chri = 22
snpsfl = sprintf("%s/%s%s.vcf.gz", geno_dir, vcf.pre, chri)

tmpfl = gzfile(snpsfl)
lns = readLines(tmpfl, n=3e4)
close(tmpfl)
lni = grep("CHROM", lns)
samvcf = gsub("#", "", unlist(strsplit(lns[lni], split="\t")))

m = match(samcov, samvcf)
lns[lni] = paste(samvcf[m], collapse="\t")
headervcf = sprintf("%s/header.vcf", rvcf.sub)
if(!file.exists(headervcf))write.table(lns[c(1:3,4)], headervcf, 
                               row.names=F, col.names=F, quote=F)
    

vcfsam = read.table(sprintf("%s/sample.txt", rvcf.sub),as.is=T)
vcfout = sprintf("%s/SNP_chr%s.vcf", rvcf.sub, chri)
nms = vcfsam[,1]
table(nms[!(nms %in% to.rm)]==vcfsam[,1])



MAFcut = 0.05
#row, column value
#vcfout = sprintf("%s/GTEx_v8_WGS_chr%s.vcf.gz", rvcf.dir, chri)
comi = sprintf("zcat %s  | sed 's/ /\t/g'  | sed /^#/d  | cut  -f '10-' | %s | cut -f '1-3'", snpsfl, pars4grNA)
p = pipe(comi)
x4u = read.table(p, colClasses=c("integer"), fill=TRUE, row.names=NULL)
x4u[1:4,]

m = match(samcov, samvcf[-(1:9)]);m
kp = x4u[,2] %in% m
x4u = x4u[kp,]
x4u[,2] = as.numeric(factor(x4u[,2], levels=m))

comi = sprintf("zcat %s  | sed 's/ /\t/g' | sed /^#/d  | cut  -f '1-3'", snpsfl)
p = pipe(comi)
pos = read.table(p, fill=TRUE, row.names=NULL)
pos[1:2,]
pos[,1] = gsub("chr", "", pos[,1])
SNPInfo = pos[, c(3,1,2)]

no.nas = aggregate(x4u[,3]==-1, by=list(x4u[,1]), FUN=any)
no.nas[1:4,]

to.kp = no.nas[no.nas[,2]==FALSE,1]
to.kp1 = x4u[,1] %in% to.kp
table(to.kp1)
x4u = x4u[to.kp1,]
x4u[,1] = as.numeric(factor(x4u[,1], levels=to.kp))
SNPInfo = SNPInfo[to.kp,]


#load back counts
ginfo = sprintf("%s/geneInfo_prepr_%s.txt", preprdir, modelL); ginfo
geneInfo = read.table(ginfo, as.is=T)
geneInfo[,2] = gsub("chr", "", geneInfo[,2])
  
ase1f = sprintf("%s/ase1_prepr_%s.txt", preprdir, modelL); ase1f
ase1 = read.table(ase1f, as.is=T)

ase2f = sprintf("%s/ase2_prepr_%s.txt", preprdir, modelL); ase2f
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

subsA = sample(1:nsub)


x4u[which(x4u[,3]==3),3] = 4
x4u[which(x4u[,3]==2),3] = 3
colnames(SNPInfo) = c("snpid", "chr", "pos")
for(nsub in numsam){
  subs = subsA[1:nsub]
  subse = c(0, subs)+1
  subsf = c(1:9, subs+9)
  subs
  subsf
  int.dir = sprintf("%s_%s_%s", pref, nsub, cis_window_char)
  if(!file.exists(int.dir))dir.create(int.dir)

  sam.ord = cvrtlo[subs,1]
  sam.ord2 = colnames(trecDs)[subs]
  table(sam.ord==sam.ord2)
  
  sam.ordf = sprintf("%s/samples.dat", int.dir)
  sam.ordf
  if(!file.exists(sam.ordf))write.table(sam.ord, sam.ordf, row.names=F, col.names=F, quote=F)

  exprL = exprL[,subs]
  exprS = exprS[,subs]
  
  trecDsi = t(trecDs[,subs])
  trecDsi[1:5,1:5]
  
  trecDli = t(trecDl[,subs])
  trecDsi[1:5,1:5]
  
  ase1i = t(ase1[,subs])
  ase1i[1:5,1:5]
  
  ase2i = t(ase2[,subs])
  ase2i[1:5,1:5]
  
  tmp = model.frame(cvrtsh[subs,-1])
  X2sh = model.matrix(~., data = tmp)
  dim(X2sh)
  Xsh = X2sh[,-1]
  dim(Xsh)
  
  tmp = model.frame(cvrtlo[subs,-1])
  X2lo = model.matrix(~., data = tmp)
  dim(X2lo)
  Xlo = X2lo[,-1]
  dim(Xlo)
  
  nr = nrow(SNPInfo)
  time1 = proc.time()
  
  eChr = rep(chri, nrow(SNPInfo))
  ePos = as.numeric((geneInfo[,3] + geneInfo[,4])/2)
  
  mChr = rep(chri, nr)
  geni = 1
  numsnps = rep(0, nrow(geneInfo));nrow(geneInfo)
  for(geni in 1:nrow(geneInfo)){
    local.distances = as.numeric((geneInfo[,4] - geneInfo[,3]))/2+cis_window
    to.kp = which(SNPInfo[,3]>=(ePos[geni]-local.distances[geni]) & SNPInfo[,3]<(ePos[geni]+local.distances[geni]));length(to.kp)
    proceed=FALSE
    if(length(to.kp)>0)proceed=TRUE
    if(proceed){
      to.kp1 = x4u[,1] %in% to.kp
      x4u0 = x4u[to.kp1,]
      x4u0[,1] = tmp = as.numeric(factor(x4u0[,1], levels=to.kp))
      SNPInfo0 = SNPInfo[to.kp,]
      geno0 = sparseMatrix(i=x4u0[,1], j=x4u0[,2], x=x4u0[,3], dims=c(length(to.kp), numsam))
      c(nrow(geno0), nrow(SNPInfo0))
      genos = geno0[,subs,,drop=F] #matrix(, ncol=nsub)
      MAF = apply(genos, 1, function(x) (sum(x %in% c(1,3)) + 2*sum(x==4)) /(2*length(x)))
      summary(MAF)
      summary(MAF>=MAFcut & MAF<=(1-MAFcut))
      geno2kp = which(MAF >= MAFcut & MAF  <= (1-MAFcut))
      numsnps[geni] = length(geno2kp)
      if(numsnps[geni]==0)proceed=FALSE
    }
    if(proceed){
      genfi = genos[geno2kp, ,drop=F]
      SNPInfoi = SNPInfo0[geno2kp, ]
      dim(genfi)
      dim(SNPInfoi)
      suffi = sprintf("%s_%s", chri, geni)
      write.table(as.matrix(genfi), sprintf("%s/genotypes_%s.dat", int.dir, suffi), row.names=F, col.names=T, quote=F, sep="\t")
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
    message(geni, " out of ", nrow(geneInfo), " #snps: ", length(geno2kp))
  }
}
summary(numsnps)

q(save = 'no')


