
library("GenomicFeatures")

gtfFile = "gencode.v26.GRCh38.genes.gtf"

path = "https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf/"

txdb = makeTxDbFromGFF(file=gtfFile, format="gtf",
  dataSource=paste(path, gtfFile, sep=""),
  organism="Homo sapiens")

saveDb(txdb, file="gencode.v26.GRCh38.genes.sqlite")

seqlevels(txdb)
columns(txdb)
keytypes(txdb)

genes = exonsBy(txdb, by="gene")
saveRDS(genes, file = "exon_by_genes_gencode.v26.GRCh38.rds")

sessionInfo()

q(save="no")
