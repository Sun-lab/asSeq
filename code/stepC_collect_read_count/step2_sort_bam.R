args = commandArgs(TRUE)
#i = 1
i = as.numeric(args[1])
i

# bsub -q day -M 32 R CMD BATCH --no-save --no-restore step2_sort_bam.R

# Sort the bam file by qname. This is because we want the paired-end reads to be next to each other,
# so that we can extract them together when we extract allele-specific reads. 
# source("https://bioconductor.org/biocLite.R")
# biocLite("Rsamtools")
library(Rsamtools)

# -----------------------------------------------------------------
# read in sample information
# -----------------------------------------------------------------

proj_output = "/fh/scratch/delete30/sun_w/COAD_tumor_eQTL/sort_bam/521"
proj_input = "/fh/scratch/delete30/sun_w/COAD_tumor_eQTL/521"

if (!dir.exists(proj_output)) dir.create(proj_output, recursive=TRUE)
setwd(proj_input)
ffs = list.files(pattern="rehead.bam$", recursive = T)
ffs


ffi = ffs[i]
cat(i, ffi, date(), "\n")
# The following lines use Rsamtools::sortBam().
fsi = strsplit(ffi, '/', fixed = T)[[1]][2]
fsi = gsub(".bam", "_filtered_sorted_byQname", fsi)
fsi = sprintf("%s/%s", proj_output, fsi)
sortBam(ffi, fsi, byQname=TRUE, maxMemory=16384)
# Alternatively, the following lines directly call samtools.
maxMemory = "16384M"
# fsi = gsub("_hg38_sorted_RG_dedup_realigned.bam", "_filtered_sorted_byQname.bam", ffi)
# fsi = sprintf("%s/%s", proj_output, fsi)
# system(sprintf("samtools sort -m %s -n %s -o %s", maxMemory, fsi, ffi))

cat('OK Sort \n')

sessionInfo()

q('no')


