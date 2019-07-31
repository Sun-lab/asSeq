args=(commandArgs(TRUE))
# args = c("22")
chri = as.numeric(args[1])
chri 

permute = T
set.seed(1565691)


# ------------------------------------------------------------
# read in data
# ------------------------------------------------------------
inputdir = "/fh/fast/sun_w/licai/_GTEx/R_batch4_whole_blood/stepC/output/"
step4a   = "step4a_ASReC_per_gene_hap1_filterIt"
step4b   = "step4b_ASReC_per_gene_hap2_filterIt"
step5    =  "step5_TReC_per_gene_filterIt"
geno_dir = "/fh/fast/sun_w/licai/_GTEx/data_genotype_all/genotype"
outdir = '/fh/fast/sun_w/licai/_GTEx/wb_eQTL_output_R_batch7'
cnt.dir = sprintf("%s/../cnt", geno_dir)

# setwd(inputdir)
if(chri == 1){
  if(permute){
    system(sprintf("echo 'chr\ttime_hr' >  %s/asSeq_permute_time.txt", outdir))
  }else{
    system(sprintf("echo 'chr\ttime_hr' >  %s/asSeq_time.txt", outdir))
  }
}

# -------------
# covariates
# -------------
geneInfo = read.table(sprintf("%s/Info_ex_chr%s.dat", cnt.dir, chri), 
                    header = T, as.is = T, sep ='\t')


covariates_file_name = sprintf("%s/../whole_blood_covs_eQTL.txt",geno_dir)
cvrt = read.table(covariates_file_name, sep ="\t", header =T)
cvrt$rd = log(cvrt$rd)
cvrt$SEX = factor(cvrt$SEX)
head(cvrt)
sam2kp = rownames(cvrt)


# -------------
# Read Counts
# -------------

samples = read.table(sprintf("%s/../phased/phasing_chr22.sample",geno_dir),
                     header=T,as.is=T)
samples = samples[-1,2]
length(samples)

SNP_file_name = sprintf("%s/chr%s.txt",geno_dir, chri)
expression_file_name = sprintf("%s/%s/gene_level_counts_filter_out_low_expressed_chr%s.txt",
                               inputdir, step5, chri)
step4a_file = sprintf("gene_level_counts_filter_out_low_expressed_hap1_chr%s.txt",
                      chri)
step4b_file = sprintf("gene_level_counts_filter_out_low_expressed_hap2_chr%s.txt",
                      chri)


trecD = read.table(expression_file_name, header = T, 
                   sep = '\t', as.is = T, check.names = F)
trecD = data.matrix(trecD)

anyNA(trecD)
trecD   = trecD[, sam2kp]
dim(trecD)
trecD[1:5,1:5]

ase1 = read.table(sprintf("%s/%s/%s",inputdir , step4a, step4a_file), as.is = T,
                  sep = "\t", check.names = F)
ase1 = data.matrix(ase1)
ase1    = ase1[, sam2kp]
dim(ase1)
ase1[1:5,1:5]
anyNA(ase1)

ase2 = read.table(sprintf("%s/%s/%s",inputdir, step4b, step4b_file), as.is = T,
                  sep = "\t", check.names = F)
ase2 = data.matrix(ase2)
ase2    = ase2[, sam2kp]
dim(ase2)
ase2[1:5,1:5]
anyNA(ase2)

# ------------------------------------------------------------
# DESeq2 remove outlier (default 99% f-distribution (p, m-p))
# ------------------------------------------------------------
library(DESeq2)

cookCutOff = 4/355 #1 #qf(0.99, 7, 348)
cookCutOff

dds <- DESeqDataSetFromMatrix(countData = (trecD+1),
                              colData = cvrt,
                              design = ~ rd + C1 + C2 + C3  + AGE + SEX)
dds = DESeq(dds)
names(assays(dds))
res <- results(dds)
summary(res)
ind = apply(assays(dds)[["cooks"]],1,function(x) sum(x>cookCutOff))
table(ind)

replaceOutliers <- function(object, trim=.2, cooksCutoff, minReplicates=1, whichSamples) {
  if (is.null(attr(object,"modelMatrix")) | !("cooks" %in% assayNames(object))) {
    stop("first run DESeq, nbinomWaldTest, or nbinomLRT to identify outliers")
  }
  if (minReplicates < 1) {
    stop("at least 3 replicates are necessary in order to indentify a sample as a count outlier")
  }
  stopifnot(is.numeric(minReplicates) & length(minReplicates) == 1)
  p <- ncol(attr(object,"modelMatrix"))
  m <- ncol(object)
  if (m <= p) {
    assays(object)[["originalCounts"]] <- counts(object)
    return(object)
  }
  if (missing(cooksCutoff)) {
    cooksCutoff <- qf(.99, p, m - p)
  }
  idx <- which(assays(object)[["cooks"]] > cooksCutoff)
  mcols(object)$replace <- apply(assays(object)[["cooks"]], 1, function(row) any(row > cooksCutoff))
  mcols(mcols(object),use.names=TRUE)["replace",] <- DataFrame(type="intermediate",description="had counts replaced")
  trimBaseMean <- apply(counts(object,normalized=TRUE),1,mean,trim=trim)
  # build a matrix of counts based on the trimmed mean and the size factors
  replacementCounts <- if (!is.null(normalizationFactors(object))) {
    as.integer(matrix(rep(trimBaseMean,ncol(object)),ncol=ncol(object)) * 
                 normalizationFactors(object))
  } else {
    as.integer(outer(trimBaseMean, sizeFactors(object), "*"))
  }
  # replace only those values which fall above the cutoff on Cook's distance
  newCounts <- counts(object)
  newCounts[idx] <- replacementCounts[idx]
  
  if (missing(whichSamples)) {
    whichSamples <- nOrMoreInCell(attr(object,"modelMatrix"), n = minReplicates)
  }
  stopifnot(is.logical(whichSamples))
  object$replaceable <- whichSamples
  mcols(colData(object),use.names=TRUE)["replaceable",] <- DataFrame(type="intermediate",
                                                                     description="outliers can be replaced")
  assays(object)[["originalCounts"]] <- counts(object)
  if (sum(whichSamples) == 0) {
    return(object)
  }
  counts(object)[,whichSamples] <- newCounts[,whichSamples,drop=FALSE]
  object
}

nOrMoreInCell <- function(modelMatrix, n) {
  numEqual <- sapply(seq_len(nrow(modelMatrix)), function(i) {
    modelMatrixDiff <- t(t(modelMatrix) - modelMatrix[i,])
    sum(apply(modelMatrixDiff, 1, function(row) all(row == 0)))
  })
  numEqual >= n
}

dds2 = replaceOutliers(dds, trim = 0.2, cooksCutoff = cookCutOff)
table(counts(dds2) == assays(dds2)[["originalCounts"]])
names(assays(dds2))

trecD = counts(dds2)-1
trecD[trecD<0] = 0 

totb = sprintf("%s/Tcnt_ex_after_replace_chr%s.bin", cnt.dir, chri)
con = file(totb, "wb")
writeBin(as.double(t(trecD)), con=con, double())
close(con)



# ------------------------------------------------------------
# create geneInfo and snpInfo file 
# ------------------------------------------------------------

geno = read.table(SNP_file_name, as.is = T)
dim(geno)
geno[, 1] = gsub('chr', "", geno[,1])
geno[1:4, 1:9]
anyNA(geno)

SNPInfo = geno[, c(5, 1, 2)]

genoSam2kp = match(colnames(trecD), samples)
geno = geno[,-c(1:5)]
geno = geno[, genoSam2kp]
geno = data.matrix(geno)
dim(geno)

if(permute){
  for(genei in 1:nrow(geno)){
    geno[genei,] = sample(geno[genei, ])
  }
}

# filter out maf < 0.05
MAF = apply(geno, 1, function(x) (sum(x %in% c(1,2)) + 2*sum(x==3)) /(2*length(x)))
geno2kp = which(MAF >= 0.05 & MAF  <= 1- 0.05)
length(geno2kp)
geno = geno[geno2kp, ]
dim(geno)

SNPInfo = SNPInfo[geno2kp, ]
colnames(SNPInfo) = c("snp","chr",'pos')
dim(SNPInfo)
head(SNPInfo)


# ------------------------------------------------------------
# asSeq2
# ------------------------------------------------------------
setwd(outdir)
trecD = t(trecD)
trecD[1:5,1:5]

ase1 = t(ase1)
ase1[1:5,1:5]

ase2 = t(ase2)
ase2[1:5,1:5]

geno = t(geno)
geno[1:5,1:5]


X2 = model.matrix(~., data = cvrt)
X2[1:2,]
# -------------
# check sample match
# -------------
library(asSeq2)
table(samples[genoSam2kp] == rownames(trecD))
table(samples[genoSam2kp] == rownames(ase1))
table(samples[genoSam2kp] == rownames(ase2))
table(samples[genoSam2kp] == rownames(X2))

if(permute){
  file_trecase = sprintf("%s/wb_asSeq2_eps5e-5_maf0.05_permute_replace4toN_chr%s.txt",
                         outdir, chri)
}else{
  file_trecase = sprintf("%s/wb_asSeq2_eps5e-5_maf0.05_replace4toN_chr%s.txt",
                         outdir, chri)
}

time1 = Sys.time()
trecase(trecD, ase1, ase2, geno, X2, 
        SNPInfo, geneInfo, 
        file_trecase = file_trecase, useASE = 1, min_ASE_total = 5L, 
        min_nASE = 5L, min_nASE_het = 5L, cis_window = 1e5,
        eps = 5e-5, max_iter = 4000)
time2 = Sys.time()

if(permute){
  system(sprintf("echo '%s\t%s' >>  %s/asSeq_permute_time.txt", 
                 chri, time2 - time1, outdir))
}else{
  system(sprintf("echo '%s\t%s' >>  %s/asSeq_time.txt", 
                 chri, time2 - time1, outdir))
}

time2 - time1 
# ------------------------------------------------------------
# trec_fast
# ------------------------------------------------------------

library(asSeq)

output.tag     = sprintf("%s/wb_trecase_fast_5e-5_chr%s",
                         outdir, chri)
if(permute){
  output.tag     = sprintf("%s/wb_asSeq1_eps5e-5_permute_replace4toN_chr%s",
                           outdir, chri)
}else{
  output.tag     = sprintf("%s/wb_asSeq1_eps5e-5_replace4toN_chr%s",
                           outdir, chri)
}
time1 = Sys.time()
geno[geno == 3] = 4
geno[geno == 2] = 3 
time1 = Sys.time()
asSeq:::trecase(trecD, ase1, ase2, X, geno, output.tag = output.tag, 
                p.cut=1,  local.distance = 1.5e5, 
                eChr = as.numeric(geneInfo[,2]), 
                ePos = as.numeric((geneInfo[,3] + geneInfo[,4])/2), 
                mChr = as.numeric(SNPInfo[,2]), mPos = as.numeric(SNPInfo[,3]),
                maxit = 4000, 
                min.AS.reads = 5, min.AS.sample = 5, min.n.het = 5)

time2 = Sys.time()
time2 - time1 
q(save = 'no')


asSeq2_replaceO = read.table("/fh/fast/sun_w/licai/_GTEx/wb_eQTL_output/wb_asSeq2_eps5e-5_maf0.05_permute_chr1.txt", header = T, as.is = T)
mean(asSeq2_replaceO$TReC_Pvalue < 0.05)
mean(asSeq2_replaceO$Joint_Pvalue<0.05, na.rm = T)
mean(asSeq2_replaceO$ASE_Pvalue<0.05, na.rm = T)
mean(asSeq2_replaceO$p.final<0.05, na.rm = T)


pdf("../figures6/permuted/asSeq2_before_and_after_replace_outlier.pdf",
    height=4, width = 8)
par(mfrow = c(1,2))
hist(asSeq2$TReC_Pvalue, main = "")
hist(asSeq2_replaceO$TReC_Pvalue, main = "")
dev.off()


for(chri in 1:22){ # step7_trecase_permute # 
  com = sprintf("sbatch -t 3-0 --exclusive R CMD BATCH '--args %s' \\
                step1_trecase.R step1_log/step1_trecase_replace4toN_chr%s.Rout\n",
                chri, chri)
  message(com)
}

for(chri in 1:22){ # step7_trecase_permute # 
  com = sprintf("sbatch -t 3-0 --exclusive R CMD BATCH '--args %s' step1_trecase.R \\
                 step1_log/step1_trecase_permute_replace4toN_chr%s.Rout\n",
                chri, chri)
  message(com)
}
