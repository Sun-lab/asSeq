#sbatch -p general -N 1 -t 07-00:00:00 -o step0_trim_long.out --mem=8g --wrap="R CMD BATCH step0_trim_long.R step0_trim_long.Rout"
#srun --mem=32g --time=8:00:00 --pty tcsh

#args=(commandArgs(TRUE))
# args = c("21", "1")
# args = c("22", "1")
#chri = as.numeric(args[1])
#chri 

#module add samtools; module add tabix

permute = F
set.seed(1565691)
pref = "Brain_Cerebellar_Hemisphere"
model = "long"

library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
gene_anno_file = "gencode.v26.GRCh38.genes.gtf"
gene_anno_file = sprintf("/pine/scr/z/h/zhabotyn/R01/GTEx/v8/Reference/%s", gene_anno_file)
gene_anno_file2 = "gencode.v26.GRCh38.genes.sqlite"
gene_anno_file2 = sprintf("/pine/scr/z/h/zhabotyn/R01/GTEx/v8/Reference/%s", gene_anno_file2)


hg38 = sprintf("/pine/scr/z/h/zhabotyn/R01/GTEx/v8/Reference/hg38.csv")
if(!file.exists(hg38)){
  gtfFile <- system.file("extdata","GTF_files",gene_anno_file,
                         package="GenomicFeatures")
  txdb2 <- makeTxDbFromGFF(file=gene_anno_file,
               dataSource="ensemblgenomes",
               organism="Homo sapiens")
  library(AnnotationDbi)
  saveDb(txdb2, gene_anno_file2)
  names(txdb2)
  #library(GenomicTools)
  #importGtf(filePath=gene_anno_file, saveObjectAsRds = TRUE)

  gns = genes(txdb2)
  write.csv(gns, hg38, row.names=F)
}
#hg38 = read.csv(hg38)
#hg38[1:2,]

#gene_anno_file = "exon_by_genes_gencode.v26.GRCh38.rds"
#gene_anno_file = sprintf("/pine/scr/z/h/zhabotyn/R01/GTEx/v8/Reference/%s", gene_anno_file)
#genes   = readRDS(gene_anno_file)

library(Matrix)
library(DESeq2)
source("helpers.R")

pars4gr = "/nas/longleaf/home/zhabotyn/progs/parser/parser4gr"
pars4grNA = "/nas/longleaf/home/zhabotyn/progs/parser/parser4grNA"
pars4grNAunf = "/nas/longleaf/home/zhabotyn/progs/parser/parser4grNAunf"
# ------------------------------------------------------------
# read in data
# ------------------------------------------------------------
root.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8"

#vcf
#/pine/scr/z/h/zhabotyn/R01/GTEx/v8/WGS_VCF
#GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz

##suppose counts are already split into subtissue
#in new format counts are provided by each sample within a subtissue
#for now - whole blood is in
#/pine/scr/w/e/weisun/_GTEx/v8/TReC_ASReC/Whole_Blood
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

#library("GenomicFeatures")
#db =loadDb(sprintf("%s/Homo_sapiens_gencode_v24lift37_nchr.sqlite", root.dir))
#names(db)
#genes = genes(db)
#exons = exonsBy(db, by=c("gene"))
##exons = unlist(exons)
#length(exons)

# -------------
# covariates
# -------------

  
  covariates_file_name = sprintf("%s/%s.v8.covariates.txt", cov.dir, pref)
  covariates_file_upd = sprintf("%s/%s.v8.covariates_upd.txt", cov.dir, pref)
  cvrt = read.table(covariates_file_name, header =T, as.is=T)
  rownames(cvrt) = cvrt[,1]; cvrt = cvrt[,-1]
  colnames(cvrt) = gsub("\\.", "-", colnames(cvrt))
  cvrt = data.frame(t(cvrt))
  dim(cvrt)
  
#  to.rm = c("GTEX-111YS", "GTEX-14C5O", "GTEX-148VJ", "GTEX-1497J", "GTEX-1AYCT", "GTEX-1S82P")
  to.rm = c("")
  dim(cvrt)
  kp.ind = !(rownames(cvrt) %in% to.rm)
  which(!kp.ind)
  cvrt = cvrt[kp.ind,]
  samcov = rownames(cvrt)
  dim(cvrt)
  #cvrt$rd = log(cvrt$rd)
#  cvrt$sex = factor(cvrt$sex)
  head(cvrt)
  dim(cvrt)

  gns = read.csv(hg38, as.is=T)
  gns = gns[order(gns$gene_id),]
  geneInfo = gns[,c(6, 1, 2, 3, 4, 5)]; 
  colnames(geneInfo) = c("Name", "chr", "start", "end", "width", "strand")
#  geneInfo = trecD[,c(1,3,4,5,2)]
  geneInfo$leno = geneInfo$end-geneInfo$start+1
  geneInfo[1:2,]  

  
  trecFls = list.files(inputdir, pattern="trecase")
  trecNms = gsub(".trecase.txt", "", trecFls)
  trecNms = gsub(sprintf("%s_", pref), "", trecNms)
#  trecNms = sapply(trecNms, get_blocks, split="_", blocks=3)  
  trecNms[1:4]
  m = match(rownames(cvrt), trecNms);
  trecFls = trecFls[m]
  trecNms = trecNms[m]
  i = 1
  for(i in 1:length(trecFls)){
    cnti = read.table(sprintf("%s/%s", inputdir, trecFls[i]), as.is=T)
    if(i == 1){
      trecD = ase1 = ase2 = aseN = matrix(0, nrow=nrow(cnti), ncol=length(trecFls))
      rownames(trecD) = rownames(ase1) = rownames(ase2) = rownames(aseN) = rownames(cnti)
    }
    trecD[,i] = cnti[,1]
    ase1[,i] = cnti[,2]
    ase2[,i] = cnti[,3]
    aseN[,i] = cnti[,4]
    if(i%%25==0)message(i, " out of ", length(trecFls))
  }
  colnames(trecD) = colnames(ase1) = colnames(ase2) = colnames(aseN) = trecNms
  table(rownames(trecD)==gns$gene_id)
  #m = match(gns$gene_id, rownames(trecD));
  
  cvrt$rd = log10(colSums(trecD))
  cutoff = 10; prop = .2
  kp = apply(trecD>=cutoff, 1, mean)>prop;table(kp)
  trecD = trecD[kp,]
  ase1 = ase1[kp,]
  ase2 = ase2[kp,]
  aseN = aseN[kp,]
  geneInfo = geneInfo[kp,]

  diag.dir = sprintf("diag_%s", model)
  if(!file.exists(diag.dir))dir.create(diag.dir)

#
#
#
fracN = aseN/(ase1+ase2+aseN)>.1 & (ase1+ase2+aseN)>10
frac1 = ase1/(ase1+ase2)
kp = (ase1+ase2)>10

frac1ag = t(sapply(1:nrow(frac1), classify, frac1=frac1, fracN=fracN, kp=kp))
frac1ag[1:4,]
colnames(frac1ag) = c("extr.h1", "hapN", "n.sam")
cut.hN = .05; min.samples = 20

for(min.samples in c(5, 10, 20)){

pdf(sprintf("%s/%s_fracWrong_by_fracExtr_%s.pdf", diag.dir, pref, min.samples), height=6, width=6)
par(mfrow=c(2,2))
kp = frac1ag[,3]>min.samples
kp1 = kp & frac1ag[,3]>min.samples & frac1ag[,1]<=0
kp1a = kp1 & frac1ag[,2]>0
#table(kp1)
plot(density(frac1ag[kp1a,2]), bty="n", main="%samp.w/hapN, no extr.hap1")#, ylim=c(0,1)
legend("topright", legend=c(sprintf("%%withHapN=%s", round(mean(frac1ag[kp1,2]>0),2)), sprintf("based on %s genes", sum(kp1))), bty="n")
  for(cut.hN in c(0.05, 0.1, 0.25)){
    kp2 = kp & frac1ag[,3]>min.samples & frac1ag[,1]>cut.hN
    kp2a = kp2 & frac1ag[,2]>0
    plot(density(frac1ag[kp2a,2]), bty="n", xlim=0:1, main=sprintf("%%samp.w/hapN, %%extr.hap1>%s", cut.hN))#, ylim=c(0,1)
    legend("topright", legend=c(sprintf("%%withHapN=%s", round(mean(frac1ag[kp2,2]>0),2)), sprintf("based on %s genes", sum(kp2))), bty="n")#for genes with extreme fraction>XXX what fraction is with hapN
    message(sum(kp), " ", sum(kp1), " ", sum(kp2))
  }
dev.off()
}
#many samples with extreme frac among samples with samples with notable frac.hN
head(frac1ag[kp2,])


for(min.samples in c(5, 10, 20)){

pdf(sprintf("%s/%s_fracExtr_byfracWrong_%s.pdf", diag.dir, pref, min.samples), height=6, width=6)
par(mfrow=c(2,2))
kp = frac1ag[,3]>min.samples
kp1 = kp & frac1ag[,3]>min.samples & frac1ag[,2]<=0
kp1a = kp1 & frac1ag[,1]>0
#table(kp1)
plot(density(frac1ag[kp1a,1]), bty="n", main="%extr.hap1, no samp.w/hapN")#, ylim=c(0,1)
legend("topright", legend=c(sprintf("%%withExtr=%s", round(mean(frac1ag[kp1,1]>0),2)), sprintf("based on %s genes", sum(kp1))), bty="n")
  for(cut.hN in c(0.01, 0.05, 0.1)){
    kp2 = kp & frac1ag[,3]>min.samples & frac1ag[,2]>cut.hN
    kp2a = kp2 & frac1ag[,1]>0
    plot(density(frac1ag[kp2a,1]), bty="n", xlim=0:1, main=sprintf("%%extr.hap1, %%samp.w/hapN>%s", cut.hN))#, ylim=c(0,1)
    legend("topright", legend=c(sprintf("%%withExtr=%s", round(mean(frac1ag[kp2,1]>0),2)), sprintf("based on %s genes", sum(kp2))), bty="n")
    message(sum(kp), " ", sum(kp1), " ", sum(kp2))
  }
dev.off()
}




frac9 = rep(0, ncol(ase1))
fracN = aseN/(ase1+ase2+aseN)>.1 & (ase1+ase2)>10#(ase1+ase2-aseN)>5
  i = 1  
  pdf(sprintf("%s/%s_diagn.pdf", diag.dir, pref), height=4, width=8)
  for(i in 1:ncol(ase1)){
  par(mfrow=c(1,2))
    frac1 = ase1[,i]/(ase1+ase2)[,i]
    kp = (ase1+ase2)[,i]>10
    plot(density(na.omit(frac1[kp])), bty="n", main=trecNms[i])
    plot(density(na.omit(frac1[which(kp&!fracN[,i])])), bty="n", main=trecNms[i])
    frac9[i] = sum((na.omit(frac1[which(kp&!fracN[,i])])-.5)>.4)
    if(i%%25==0)message(i, " out of ", length(trecNms))
  }
  dev.off()

pdf(sprintf("%s/%s_extreme.pdf", diag.dir, pref), height=4, width=8)
par(mfrow=c(1,2))
hist(frac9, main=pref)
legend("topright", legend=ncol(ase1), bty="n")
plot(log10(colSums(ase1+ase2)+1), frac9, bty="n", main=pref, xlab="#ase", ylab="num extr")
dev.off()


  pdf(sprintf("%s/%s_ase_reads.pdf", diag.dir, pref), height=4, width=4)
    x = colSums(trecD)
    y = colSums(ase1+ase2)
    plot(log10(x+1), log10(y+1), bty="n", xlab="total", ylab="allele specific", main="")
  dev.off()
  kp = log10(y+10)>4; table(kp)
  trecNms[!kp]
  #no ASE counts
  #"GTEX-1AYD5"
  #"GTEX-1S831"
  trecD = trecD[,kp]
  ase1 = ase1[,kp]
  ase2 = ase2[,kp]
  aseN = aseN[,kp]
  cvrt = cvrt[kp,]
  

  fracN = aseN/(ase1+ase2+aseN)>.1 & (ase1+ase2)>10#(ase1+ase2-aseN)>5
  fracs = apply(fracN, 1, sum, na.rm=T)
  o = order(fracs, decreasing=T)
  fracs[1:5]
  table(fracs)
  c(sum(fracs == 0), sum(fracs == 1), sum(fracs == 2), sum(fracs %in% 3:5), sum(fracs %in% 6:10), sum(fracs %in% 11:1000), ncol(ase1))
  ex1 = fracs[o[1]]
  ex2 = fracs[o[2]]
  ex3 = fracs[o[3]]
  ex4 = fracs[o[4]]
   
  pdf(sprintf("%s/%s_top_wrong.pdf", diag.dir, pref), height=4, width=8)
  par(mfrow=c(1,2))
  for(exi in c(ex1, ex2, ex3, ex4)){
  frac1 = ase1[fracs==exi,]/(ase1+ase2+aseN)[fracs==exi,]
  frac2 = ase2[fracs==exi,]/(ase1+ase2+aseN)[fracs==exi,]
  fracN = aseN[fracs==exi,]/(ase1+ase2+aseN)[fracs==exi,]
  frac1s = ase1[fracs==exi,]/(ase1+ase2)[fracs==exi,]
  frac2s = ase2[fracs==exi,]/(ase1+ase2)[fracs==exi,]
  plot(frac1s, fracN, main=rownames(ase1)[exi], bty="n")
#  plot(frac2s, fracN, main=rownames(ase1)[exi], bty="n")
  kp = frac1s<=.5
  plot(density(na.omit(frac1s[kp]),from=0, to=.5), main=rownames(ase1)[exi], bty="n", xlim=0:1, ylim=c(0,4))
  lines(density(na.omit(frac1s[!kp]),from=.5,to=1))
  }
  dev.off()


  fracN = aseN/(ase1+ase2+aseN)>.1 & (ase1+ase2+aseN)>10#(ase1+ase2-aseN)>5
  fracNp = apply(fracN, 1, mean)
  summary(fracNp>0.05)
  summary(c(fracN))
  
  #The whole gene if fraction of samples with non-zero hN is >0.05
  #The whole gene if fraction of samples with non-zero hN is >0.01 and percent of samples with h1E is bigger then 20%
  #Any particular sample-gene count if hN>0.1
  cut1 = 0.05
  aseo = ase1+ase2
  #replace all the cases with too high fraction of bad ASE samples to zero ASE
  rm1 = fracNp>cut1; table(rm1)
  ase1[rm1,] = 0
  ase2[rm1,] = 0
  
  #replace all the cases with too many extreme hap1 if fractoin of bad ASE samples is at least 1%
  cut2 = 0.01
  cut3 = 0.20
  fracE = apply(abs(ase1/(ase1+ase2)-0.5)>.4 & (ase1+ase2)>10, 1, mean)
  rm2 = fracNp>cut2 & fracE>cut3; table(rm2)
  table(rm2)
  ase1[rm2,] = 0
  ase2[rm2,] = 0
  
  #replace all the cases  
  fracN = aseN/(ase1+ase2+aseN)>.1 & (ase1+ase2+aseN)>10#(ase1+ase2-aseN)>5
  ase1[fracN] = 0
  ase2[fracN] = 0
  table(c(aseo!=(ase1+ase1)))
  ase1[1:4,1:4]
  ase2[1:4,1:4]
  
  filt = c(cut1=sum(rm1), cut2=sum(rm2), allrm=sum(aseo!=(ase1+ase2)))
  write.csv(filt, sprintf("%s/%s_filt.csv", diag.dir, pref), quote=F)
  filt

#
#
#



  samcov = colnames(trecD)
  samout = sprintf("%s/sample.txt", rvcf.sub)
  if(!file.exists(samout))write.table(samcov, samout, row.names=F, col.names=F, quote=F)

  # ------------------------------------------------------------
  # DESeq2 remove outlier (default 99% f-distribution (p, m-p))
  # ------------------------------------------------------------
  library(DESeq2)
  
  cookCutOff = 4/nrow(cvrt) #1 #qf(0.99, 7, 348)
  cookCutOff
  #design = model.matrix(~., cvrt)
  trecDc = matrix(as.numeric(unlist(trecD)), nrow=nrow(trecD))

  #design = formula(paste("~", paste(c(colnames(cvrt)), collapse=" + ")), collapse=" ")

#  cov.to.rm = grep("InferredCov", colnames(cvrt))
#  covnms = colnames(cvrt)[-cov.to.rm]
  covnms = colnames(cvrt)

  design = formula(paste("~", paste(covnms, collapse=" + ")), collapse=" ")
 
  cvrt$pcr = as.factor(cvrt$pcr)
  cvrt$platform = as.factor(cvrt$platform)
  cvrt$sex = as.factor(cvrt$sex)
  dds <- DESeqDataSetFromMatrix(countData = (trecDc+1),
                                colData = cvrt,
                                design = design)
  dds = DESeq(dds)
  names(assays(dds))
  res <- results(dds)
  summary(res)
  ind = apply(assays(dds)[["cooks"]],1,function(x) sum(x>cookCutOff))
  table(ind)
  
  #save.image("tmp.Rdata");q("no")
  #load("/pine/scr/z/h/zhabotyn/R01/tmp.Rdata");library(DESeq2)
  
  
  dds2 = replaceOutliers(dds, trim = 0.2, cooksCutoff = cookCutOff)
  table(counts(dds2) == assays(dds2)[["originalCounts"]])
  names(assays(dds2))
  
  trecD = counts(dds2)-1
  trecD[trecD<0] = 0 
  trecD[1:2,]


  o = order(factor(geneInfo$chr, levels=sprintf("chr%s",c(1:22, "X", "Y", "M"))), geneInfo$start)
  geneInfo = geneInfo[o,]
  trecD = trecD[o,]
  ase1 = ase1[o,]
  ase2 = ase2[o,]
  aseN = aseN[o,]
 
  cvrtf = sprintf("%s/covariates_%s.txt", preprdir, model); cvrtf
  write.csv(cvrt[,covnms], cvrtf)

  ginfo = sprintf("%s/geneInfo_prepr_%s.txt", preprdir, model); ginfo
  write.table(geneInfo, ginfo, quote=F, row.names=T, col.names=T)

  ase1f = sprintf("%s/ase1_prepr_%s.txt", preprdir, model); ase1f
  write.table(ase1, ase1f, quote=F, row.names=T, col.names=T)

  ase2f = sprintf("%s/ase2_prepr_%s.txt", preprdir, model); ase2f
  write.table(ase2, ase2f, quote=F, row.names=T, col.names=T)

  aseNf = sprintf("%s/aseN_prepr_%s.txt", preprdir, model); aseNf
  write.table(aseN, aseNf, quote=F, row.names=T, col.names=T)

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
      Ks[1:5,1:5]
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


nrow(trecD)



q(save = 'no')


