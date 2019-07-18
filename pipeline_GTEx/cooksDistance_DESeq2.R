
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

q("no")
