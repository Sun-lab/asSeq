
library(biomaRt)
library(AnnotationDbi)

# -----------------------------------------------------------------
# select all the transcripts per gene
# -----------------------------------------------------------------

txdb = loadDb("gencode.v26.GRCh38.genes.sqlite")

seqlevels(txdb)
columns(txdb)
keytypes(txdb)

cols   = c("TXCHROM", "TXSTRAND", "TXSTART", "TXEND")
enids  = keys(txdb, keytype="GENEID")
length(enids)
enids[1:5]

ptm    = proc.time()
annots = select(txdb, keys=enids, columns=cols, keytype="GENEID")
proc.time() - ptm

dim(annots)
head(annots)

table(annots$TXEND - annots$TXSTART > 0)

length(unique(annots$GENEID))

# -----------------------------------------------------------------
# generate gene-level annotation
# -----------------------------------------------------------------

chrs    = tapply(annots$TXCHROM, annots$GENEID, unique)
strands = tapply(annots$TXSTRAND, annots$GENEID, unique)

is.list(chrs)
is.list(strands)

head(chrs)
head(strands)

starts = tapply(annots$TXSTART, annots$GENEID, min)
ends   = tapply(annots$TXEND, annots$GENEID, max)

table(names(chrs) == names(strands))
table(names(chrs) == names(starts))
table(names(chrs) == names(ends))

anno1 = data.frame(geneId = names(chrs), chr=chrs, strand=strands,
  start=starts, end=ends, row.names=NULL, stringsAsFactors=FALSE)

dim(anno1)
head(anno1)

pdf("hist_gene_length.pdf", width=4, height=3)
par(mar=c(5,4,1,1))
hist(log10(anno1$end - anno1$start), xlab="log10(gene length)", main="")
dev.off()

# -----------------------------------------------------------------
# get more information from biomart
# -----------------------------------------------------------------

enids = strsplit(anno1$geneId, split=".", fixed=TRUE)
table(sapply(enids, length))
enids = matrix(unlist(enids), byrow=TRUE, ncol=2)

dim(enids)
enids[1:5,]

length(unique(enids[,1]))

enids = enids[,1]
enids[1:5]

anno1$ensembl_gene_id = enids

mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")

ptm   = proc.time()
anno2 = getBM(c("hgnc_symbol","description","chromosome_name","band", 
                "strand","start_position","end_position","ensembl_gene_id"),
  filters="ensembl_gene_id", values=enids, mart)
proc.time() - ptm

dim(anno2)
anno2[1:2,]

dim(anno1)
anno1[1:2,]
length(unique(anno1$ensembl_gene_id))

table(anno2$chromosome_name)

tbls = table(anno2$ensembl_gene_id)
tbls = sort(tbls, decreasing=TRUE)

tbls[1:15]

anno2[anno2$ensembl_gene_id=="ENSG00000187510",]
anno2[anno2$ensembl_gene_id=="ENSG00000230417",]
anno2[anno2$ensembl_gene_id=="ENSG00000276085",]

# -----------------------------------------------------------------
# obtain merged annotation
# -----------------------------------------------------------------

anno3 = merge(anno1, anno2, by="ensembl_gene_id", all=TRUE)

dim(anno3)
anno3[1:2,]

table(anno3$strand.x, anno3$strand.y, useNA="ifany")
table(anno3$chr == paste("chr", anno3$chromosome_name, sep=""))

ww1 = which(anno3$chr != paste("chr", anno3$chromosome_name, sep=""))
table(anno3[ww1,"chr"], anno3[ww1,"chromosome_name"], useNA="ifany")

summary(anno3$start_position - anno3$start)
summary(anno3$end_position - anno3$end)

fun1 <- function(v){paste(unique(v), collapse=";")}

hgnc_symbol = tapply(anno3$hgnc_symbol, anno3$geneId, fun1)
description = tapply(anno3$description, anno3$geneId, fun1)

table(names(hgnc_symbol) == anno1$geneId)
table(names(description) == anno1$geneId)

anno1$hgnc_symbol = hgnc_symbol
anno1$description = description

dim(anno1)
anno1[1:2,]

write.table(anno1, file = "gencode.v26.GRCh38.genes_gene_level_anno.txt", 
            append = FALSE, quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)

gc()

sessionInfo()

q(save="no")

