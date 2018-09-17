# Count total reads using function summarizeOverlaps from R package GenomicAlignments.
args = commandArgs(TRUE)
#i = 1
i = as.numeric(args[1])
i

meta1 = read.table('/fh/scratch/delete30/sun_w/COAD_tumor_eQTL/200/gdc_sample_sheet.2018-06-29.tsv',
                   sep = '\t', header = T, as.is = T)
meta1[1:2,]

barcode = meta1$Sample.ID

library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")

setwd("/fh/fast/sun_w/licai/_tumor_eQTL/R_batch1/stepC_collect_read_count/pipeline")
pipeline_dir = getwd()

db_input = "/fh/fast/sun_w/licai/_tumor_eQTL"
txdb = loadDb(sprintf("%s/Homo_sapiens_gencode_v28_GRCh38.sqlite", db_input))
seqlevels(txdb)
columns(txdb)
keytypes(txdb)

genes = exonsBy(txdb, by="gene")

proj_output = "/fh/fast/sun_w/licai/_tumor_eQTL/R_batch1/stepC_collect_read_count/output/step5_TReC_per_gene_filterIt"

if (!dir.exists(proj_output)) dir.create(proj_output)
proj_input = "/fh/scratch/delete30/sun_w/COAD_tumor_eQTL/sort_bam/"
#setwd(proj_input)
#filenames = list.files(pattern="_prepBAM_uniq_filtered.bam$")
setwd(proj_input)
filenames = list.files(pattern="_filtered_sorted_byQname.bam$")[i]

ind = grep(substr(filenames, 1,38), meta1$File.Name)
sams = barcode[ind]
sams

# This seems not to be working...
library(BiocParallel)
#register(MulticoreParam(workers = 1), default=TRUE)
register(SerialParam(), default=TRUE)
registered()
bpparam()
bamfiles  = BamFileList(filenames, yieldSize=1000000)
# "Smaller chunks consume less memory,  but are a little less efficient to process." From the vignette
# bamfiles  = BamFileList(filenames, yieldSize=1000000)
bamfiles

cts = rep(0, length(bamfiles))

for(i in 1:length(bamfiles)){
  
  cat(i, date(), "\n")
  
  bami = bamfiles[i]
  ct1    = countBam(bami)
  cts[i] = ct1$records
}

se = summarizeOverlaps(features=genes, reads=bamfiles, mode="Union",
                       singleEnd=TRUE, ignore.strand=FALSE, fragments=FALSE )

se
colData(se)
#rowData(se)  # old version of Bioconductor
rowRanges(se)
str(metadata(rowRanges(se)))

as1 = assay(se)

head(as1)
cts = rbind(cts, colSums(as1))

setwd(pipeline_dir)
colnames(as1) = sams
write.table(as1, file = sprintf("%s/%s_gene_level_counts_filterIt_total.txt", proj_output,sams), append = FALSE,
            quote = FALSE, sep = "\t", eol = "\n", row.names = TRUE,
            col.names = TRUE)

row.names(cts) = c("nMappedReads", "nFragMappedToGene")
write.table(cts, file = sprintf("%s/%s_TReC_liberal_filterIt_total.txt", proj_output,sams), append = FALSE,
            quote = FALSE, sep = "\t", eol = "\n", row.names = TRUE,
            col.names = TRUE)

sessionInfo()
q(save="no")

