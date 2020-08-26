#module add r;module add samtools; module add tabix

#sbatch -p general -N 1 -t 07-00:00:00 -o step0_trim_long.out --mem=8g --wrap="R CMD BATCH step0_trim_long.R step0_trim_long.Rout"
#srun --mem=32g --time=8:00:00 --pty tcsh

#args=(commandArgs(TRUE))
#specf = args[1]
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
rcmd = specs[19]
setwd(wrk.dir)

set.seed(seedval)
model = "short"

library(Matrix)
library(DESeq2)
source(sprintf("%s/helpers.R", lib.dir))

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

hg38 = sprintf("%s/Reference/hg38.csv", bas.dir)
gns = read.csv(hg38, as.is=T)
gns = gns[order(gns$gene_id),]
geneInfo = gns[,c(6, 1, 2, 3, 4, 5)]; 
colnames(geneInfo) = c("Name", "chr", "start", "end", "width", "strand")
geneInfo$leno = geneInfo$end-geneInfo$start+1

# -------------
# covariates
# -------------
covariates_file_name = sprintf("%s/%s.v8.covariates.txt", cov.dir, pref)
covariates_file_upd = sprintf("%s/%s.v8.covariates_upd.txt", cov.dir, pref)
cvrt = read.table(covariates_file_name, header =T, as.is=T)
rownames(cvrt) = cvrt[,1]; cvrt = cvrt[,-1]
colnames(cvrt) = gsub("\\.", "-", colnames(cvrt))
cvrt = data.frame(t(cvrt))
rm.sex = all(cvrt$sex==1)|all(cvrt$sex==0);rm.sex
dim(cvrt)
if(rm.sex){cvrt = cvrt[,colnames(cvrt)!="sex"]}
dim(cvrt)

totp = sprintf("%s/Tcnt_ori.txt", preprdir); totp
trecD = read.table(totp, as.is=T)
colnames(trecD) = gsub("\\.","-", colnames(trecD))
colnames(trecD)[1:4]
m = match(colnames(trecD), rownames(cvrt))
cvrt = cvrt[m,]
table(colnames(trecD)==rownames(cvrt))
mi = match(rownames(trecD), geneInfo$Name)
geneInfo = geneInfo[mi,]
table(rownames(trecD)==geneInfo$Name)

cookCutOff = 4/nrow(cvrt) #1 #qf(0.99, 7, 348)
cookCutOff

trecDc = matrix(as.numeric(unlist(trecD)), nrow=nrow(trecD))

cov.to.rm = grep("InferredCov", colnames(cvrt))
covnms = colnames(cvrt)[-cov.to.rm]
design = formula(paste("~", paste(covnms, collapse=" + ")), collapse=" ")
 
cvrt$pcr = as.factor(cvrt$pcr)
cvrt$platform = as.factor(cvrt$platform)
if(!rm.sex)cvrt$sex = as.factor(cvrt$sex)
dds <- DESeqDataSetFromMatrix(countData = (trecDc+1),
                                colData = cvrt,
                                design = design)
dds = DESeq(dds)
names(assays(dds))
res <- results(dds)
summary(res)
ind = apply(assays(dds)[["cooks"]],1,function(x) sum(x>cookCutOff))
table(ind)
  
dds2 = replaceOutliers(dds, trim = 0.2, cooksCutoff = cookCutOff)
table(counts(dds2) == assays(dds2)[["originalCounts"]])
names(assays(dds2))

diag.dir = sprintf("%s_diag", pref)
if(!file.exists(diag.dir))dir.create(diag.dir)
filt = table(counts(dds2) == assays(dds2)[["originalCounts"]])
pref1 = sprintf("%s_%s", pref, model)
write.csv(filt, sprintf("%s/%s_updtot.csv", diag.dir, pref1), quote=F, row.names=F)
  
trecD = counts(dds2)-1
trecD[trecD<0] = 0 
trecD[1:2,]


o = order(factor(geneInfo$chr, levels=sprintf("chr%s",c(1:22, "X", "Y", "M"))), geneInfo$start)
geneInfo = geneInfo[o,]
trecD = trecD[o,]
ginfoL = sprintf("%s/geneInfo_prepr_%s.txt", preprdir, "long"); ginfoL
geneInfoL = read.table(ginfoL, as.is=T)
table(geneInfoL$Name==geneInfo$Name)
 
cvrtf = sprintf("%s/covariates_%s.csv", preprdir, model); cvrtf
write.csv(cvrt[,covnms], cvrtf)

ginfo = sprintf("%s/geneInfo_prepr_%s.txt", preprdir, model); ginfo
write.table(geneInfo, ginfo, quote=F, row.names=T, col.names=T)

totp = sprintf("%s/Tcnt_prepr_%s.txt", preprdir, model); totp
write.table(trecD, totp, quote=F, row.names=T, col.names=T)

totb = sprintf("%s/Tcnt_prepr_%s.bin", preprdir, model)
con = file(totb, "wb")
writeBin(as.double(t(trecD)), con=con, double())
close(con)
  
  
kbin = sprintf("%s/K_ex_%s.bin", preprdir, model)
xbin = sprintf("%s/X_ex_%s.bin", preprdir, model)
  
if(!file.exists(xbin)){
  con = file(xbin, "wb")
  writeBin(as.double(c(as.matrix(cvrt[,covnms]))), con=con, double())
  close(con)
}
  
if(!file.exists(kbin)){
  Ks = 10^cvrt$rd
  Ks = matrix(rep(Ks/mean(Ks),each=nrow(trecD)), nrow=nrow(trecD))
  con = file(kbin,"wb")
  writeBin(as.double(c(t(Ks))), con=con, double())
  close(con)
}
  
kin = readBin(kbin, n=200, double())
kin
tin = readBin(totb, n=200, double())
tin
xin = readBin(xbin, n=200, double())
xin

q(save = 'no')


