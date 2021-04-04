rm(list = ls())
setwd("/fh/fast/sun_w/licai/_tumor_eQTL/GitHub/asSeq/pipeline_GTEx/v8/cell_type_composition/")

calculate_TPM = function(count, gene_length) {
  if (nrow(count) != length(gene_length)) {
    stop("Number of rows of the count matrix does not match gene lengths!")
  }
  TPM = count / gene_length
  t(t(TPM)*1e6/colSums(TPM))
}

# -----------------------------------------------------------------
# Read in total read count
# -----------------------------------------------------------------
TReC = read.table("data/Tcnt_prepr_long.txt", as.is = T, header = T, 
                  check.names = F)
TReC[1:5, 1:5]
dim(TReC)

Info = read.table("data/geneInfo_prepr_long.txt", as.is = T, header = T, 
                  check.names = F )
Info[1:5,]
dim(Info)

rownames(TReC) = Info$Name

# -----------------------------------------------------------------
# gene name and length information 
# -----------------------------------------------------------------

geneInfo = read.table("../Reference/gencode.v26.GRCh38.genes_gene_level_anno.txt",
                      header = T, as.is = T, sep="\t")
head(geneInfo)

geneLength = readRDS("../Reference/gene_exonic_length.v26.GRCh38.rds")
geneLength[1:5]
length(geneLength)
geneLength = unlist(geneLength)

table(rownames(TReC) %in% names(geneLength))

gene_length = geneLength[rownames(TReC)]
gene_length[1:5]

# -----------------------------------------------------------------
# Transform count data to TPM data 
# -----------------------------------------------------------------
table(rownames(TReC) == names(gene_length))

TPM = calculate_TPM(TReC, gene_length)
TPM[1:5,1:5]
dim(TPM)


table(rownames(TPM) %in%  geneInfo$geneId) 

# -----------------------------------------------------------------
# hgnc_symbol
# -----------------------------------------------------------------

hgnc = geneInfo$hgnc_symbol[match(rownames(TPM), geneInfo$geneId)]
sum(duplicated(hgnc))
sum(is.na(hgnc) | hgnc == "")

gene2rm = which(is.na(hgnc) | hgnc == "")
length(gene2rm)
sum(duplicated(hgnc[-gene2rm]))

TPM = TPM[-gene2rm, ]
rownames(TPM) = hgnc[-gene2rm]
dim(TPM)
TPM[1:5,1:5]



# -----------------------------------------------------------------
# reference data
# -----------------------------------------------------------------

file = "/fh/fast/sun_w/licai/cell_type_association_summary/act/TCGA_COAD/data/LM22.txt"

LM22 = read.table(file, header = T, sep = "\t", as.is = T)
LM22[1:5,1:5]

table(LM22$Gene.symbol %in% rownames(TPM))
table(LM22$Gene.symbol %in% geneInfo$hgnc_symbol)

TPM = cbind(Gene.symbol = rownames(TPM), TPM)
TPM = TPM[which(rownames(TPM) %in% LM22$Gene.symbol), ]
write.table(TPM, file = "data/TPM.txt", quote = FALSE,
            row.names = F, sep="\t")


q(save = "no")