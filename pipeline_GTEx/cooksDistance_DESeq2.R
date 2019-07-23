
chri = 22

proj_dir = "/fh/fast/sun_w/licai/_GTEx/wb_eQTL_output/"

# ------------------------------------------------------------
# read in data
# ------------------------------------------------------------
inputdir = "/fh/fast/sun_w/licai/_GTEx/R_batch4_whole_blood/stepC/output/"
step5    =  "step5_TReC_per_gene_filterIt"
geno_dir = "/fh/fast/sun_w/licai/_GTEx/data_genotype_all/genotype"
outdir = '/fh/fast/sun_w/licai/_GTEx/wb_eQTL_output'


# covariates 
covariates_file_name = sprintf("%s/../whole_blood_covs_eQTL.txt",geno_dir)
cvrt = read.table(covariates_file_name, sep ="\t", header =T)
cvrt$rd = log(cvrt$rd)
cvrt$SEX = factor(cvrt$SEX)
head(cvrt)
sam2kp = rownames(cvrt)

# total read counts
expression_file_name = sprintf("%s/%s/gene_level_counts_filter_out_low_expressed_chr%s.txt",
                               inputdir, step5, chri)

trecD = read.table(expression_file_name, header = T, 
                   sep = '\t', as.is = T, check.names = F)
trecD = data.matrix(trecD)

anyNA(trecD)
trecD   = trecD[, sam2kp]
dim(trecD)
trecD[1:5,1:5]


library(DESeq2)

cookCutOff = 4/355 #1 #qf(0.95, 7, 348)
cookCutOff

dds <- DESeqDataSetFromMatrix(countData = (trecD+1), # add 1 to avoid zero counts
                              colData = cvrt,
                              design = ~ rd + C1 + C2 + C3  + AGE + SEX)
dds = DESeq(dds)
names(assays(dds))
# res <- results(dds)
# summary(res)
ind = apply(assays(dds)[["cooks"]],1,function(x) sum(x>cookCutOff))
table(ind)

replaceOutliers <- function(object, trim=.2, cooksCutoff, minReplicates=1, whichSamples) {
  if (is.null(attr(object,"modelMatrix")) | !("cooks" %in% assayNames(object))) {
    stop("first run DESeq, nbinomWaldTest, or nbinomLRT to identify outliers")
  }
  if (minReplicates < 1) { ## changed minReplicates to 1
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

q("no")
