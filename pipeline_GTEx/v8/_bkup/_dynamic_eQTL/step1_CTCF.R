
library(data.table)
library(GenomicRanges)

# --------------------------------------------------------------------
# read in resutls
# --------------------------------------------------------------------

ctcf_file = "long_quasi_Whole_Blood_ctcf.csv"
ctcf = fread(file.path("../results/GTEx8_Whole_Blood_summary/", ctcf_file))
dim(ctcf)
ctcf[1:2,]

summary(ctcf$p_cond)
pi0 = 2*mean(ctcf$p_cond > 0.5)
pi0

qval = nrow(ctcf)*pi0*ctcf$p_cond/rank(ctcf$p_cond)
table(ctcf$p_cond < 0.01)
table(qval < 0.01)
table(qval < 0.05)
table(qval < 0.10)
table(qval < 0.15)
table(qval < 0.20)

ctcf$qval = qval

# --------------------------------------------------------------------
# read in gene annoation information
# --------------------------------------------------------------------

gene_file = "gencode.v26.GRCh38.genes_gene_level_anno.txt"
gene_file = file.path("../Reference", gene_file)
genes = fread(gene_file)
dim(genes)
genes[1:2,]

table(ctcf$id %in% genes$geneId)
mat1 = match(ctcf$id, genes$geneId)

ctcf = cbind(ctcf, genes[mat1,])
dim(ctcf)
ctcf[1:2,]

ctcf$promoter_start = rep(NA, nrow(ctcf))
ctcf$promoter_end   = rep(NA, nrow(ctcf))

wn = which(ctcf$strand == "-")
wp = which(ctcf$strand == "+")

ctcf$promoter_start[wn] = ctcf$end[wn]
ctcf$promoter_end[wn]   = ctcf$end[wn] + 199

ctcf$promoter_start[wp] = ctcf$start[wp] - 199
ctcf$promoter_end[wp]   = ctcf$start[wp]

dim(ctcf)
ctcf[1:5,]

gr1 = makeGRangesFromDataFrame(ctcf, ignore.strand=TRUE, 
                               seqnames.field="chr",
                               start.field="promoter_start", 
                               end.field="promoter_end")

# --------------------------------------------------------------------
# read in CTCF binding site information
# --------------------------------------------------------------------

ff1  = "CTCFBSDB_all_exp_sites_Sept12_2012_hg38_loci.bed.gz"
ff1  = file.path("../Reference/CTCF/", ff1)

ctcf.bs = fread(ff1)
dim(ctcf.bs)
ctcf.bs[1:5,]

names(ctcf.bs) = c("chr", "start", "end")
table(ctcf.bs$chr)

lens = ctcf.bs$end - ctcf.bs$start + 1
summary(lens)

pdf("../Reference/CTCF/CTCFBS_len_hist.pdf", width=6, height=4)
par(mar=c(5,4,1,1), bty="n")
hist(log10(lens), xlab="log10(CTCF BS length)", main="", breaks=100)
abline(v=log10(200))
dev.off()

table(lens < 400)/length(lens)
table(lens < 300)/length(lens)
table(lens < 200)/length(lens)

gr2 = makeGRangesFromDataFrame(ctcf.bs, ignore.strand=TRUE, 
                               seqnames.field="chr",
                               start.field="start", 
                               end.field="end")

gr3 = makeGRangesFromDataFrame(ctcf.bs[which(lens < 200),], 
                               ignore.strand=TRUE, 
                               seqnames.field="chr",
                               start.field="start", 
                               end.field="end")

fun1 <- function(x) sum(width(reduce(x, ignore.strand=T)))
fun1(gr1)

width2 = fun1(gr2)
width3 = fun1(gr3)
width2
width3

prp2 = width2/(3234.83*10^6)
prp3 = width3/(3234.83*10^6)
prp2
prp3

mtch2 = findOverlaps(gr1, gr2, select="first")
table(!is.na(mtch2))

mtch3 = findOverlaps(gr1, gr3, select="first")
table(!is.na(mtch3))

table(qval < 0.1, !is.na(mtch2))
table(qval < 0.1, !is.na(mtch3))

chisq.test(qval < 0.05, !is.na(mtch2))
chisq.test(qval < 0.05, !is.na(mtch3))

chisq.test(qval < 0.1, !is.na(mtch2))
chisq.test(qval < 0.1, !is.na(mtch3))

gc()
sessionInfo()
q(save = "no")


