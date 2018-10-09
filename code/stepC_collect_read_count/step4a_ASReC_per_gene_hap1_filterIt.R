#Count allele-specific reads using function summarizeOverlaps from R package GenomicAlignments.

# sbatch -t 5-0 --mem=32000 R CMD BATCH --no-save --no-restore step4a_ASReC_per_gene_hap1_filterIt.R

library("GenomicAlignments")
library("GenomicFeatures")

pipeline_dir = getwd()

db_input = "/fh/fast/sun_w/licai/_tumor_eQTL"
txdb = loadDb(sprintf("%s/Homo_sapiens_gencode_v28_GRCh38.sqlite", db_input))
seqlevels(txdb)
columns(txdb)
keytypes(txdb)

genes = exonsBy(txdb, by="gene")

proj_input = "../../../data_cloud/step3_extract_asReads_filterIt"
proj_output = "../output/step4a_ASReC_per_gene_hap1_filterIt"

if (!dir.exists(proj_output)) dir.create(proj_output)
setwd(proj_input)

filenames = list.files(pattern="_filterIt_hap1.bam$")
bamfiles  = BamFileList(filenames, yieldSize=1000000)
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
rowRanges(se)
str(metadata(rowRanges(se)))

as1 = assay(se)
head(as1)
cts = rbind(cts, colSums(as1))

setwd(pipeline_dir)
write.table(as1, file = sprintf("%s/gene_level_counts_filterIt_hap1.txt", proj_output), append = FALSE,
            quote = FALSE, sep = "\t", eol = "\n", row.names = TRUE,
            col.names = TRUE)

row.names(cts) = c("nMappedReads", "nFragMappedToGene")
write.table(cts, file = sprintf("%s/ASReC_liberal_filterIt_hap1.txt", proj_output), append = FALSE,
            quote = FALSE, sep = "\t", eol = "\n", row.names = TRUE,
            col.names = TRUE)

sessionInfo()
q(save="no")

