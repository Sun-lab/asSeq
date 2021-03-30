
library(parallel)

# -----------------------------------------------------------------
# select all the transcripts per gene
# -----------------------------------------------------------------

gene_anno = readRDS("exon_by_genes_gencode.v26.GRCh38.rds")
length(gene_anno)
names(gene_anno)[1]
gene_anno[[1]]

date()
gene_exonic_length = mclapply(gene_anno,
                              function(x){sum(width(reduce(x)))},
                              mc.cores=8)
date()

gene_exonic_length[1:5]

saveRDS(gene_exonic_length, file = "gene_exonic_length.v26.GRCh38.rds")

q(save="no")

