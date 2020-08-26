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
pref = "Brain_Caudate_basal_ganglia"
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





q(save = 'no')


