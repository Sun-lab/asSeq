args = commandArgs(TRUE)
#i = 1
i = as.numeric(args[1])
i
# ----------------------------------------------------------------------
# match filename with TCGA barcode
# ----------------------------------------------------------------------

meta1 = read.table('/fh/scratch/delete30/sun_w/COAD_tumor_eQTL/200/gdc_sample_sheet.2018-06-29.tsv',
                   sep = '\t', header = T, as.is = T)
meta1[1:2,]

barcode = meta1$Sample.ID

library(asSeq)
setwd('/fh/fast/sun_w/licai/_tumor_eQTL/R_batch1/stepC_collect_read_count/pipeline')
proj_pipeline =getwd()
proj_input = "/fh/scratch/delete30/sun_w/COAD_tumor_eQTL/sort_bam/521"
proj_output = "../output/step3_extract_asReads_filterIt"
proj_prepareBAM = "/fh/scratch/delete30/sun_w/COAD_tumor_eQTL/prepare_bam"
if (!dir.exists(proj_output)) dir.create(proj_output)
if (!dir.exists(proj_prepareBAM)) dir.create(proj_prepareBAM)
setwd(proj_input)
files = list.files(pattern="_filtered_sorted_byQname.bam")
#files = list.files(pattern="_prepBAM_uniq_filtered.bam")
length(files)
files

sams  = gsub("_filtered_sorted_byQname.bam", "", files)
#sams  = gsub("_prepBAM_uniq_filtered.bam", "", files)
sams

ind = sapply(sams, function(x) grep(x, meta1$File.Name))
sams = barcode[ind]
sams


sam1 = sams[i]
cat("\n", sam1, date(), "\n")
  #snpList = sprintf("../%s_hetSNP.txt", sam1)
  #snpFolderPrefix = "/lustre/scr/c/h/chongjin/single_cell/stepB_phasing_and_imputation/output/step15_hg38/snps/ind1_SRR1163508"
  #snpFolder = dir(snpFolderPrefix, pattern=sam1)
  #snpFolder = paste0(snpFolderPrefix, snpFolder)
snpFolder = "/fh/fast/sun_w/licai/_tumor_eQTL/data_snp"
snpFolder
tmpSnpFile = tempfile()
system(sprintf("cd %s/%s ;  cat combined_hg38.txt | cut -f1,3,4,5  >> %s", snpFolder, substr(sams[i], 1, 12), tmpSnpFile)) # hg38
#system(sprintf("cd %s/%s ; rm combined.txt; for i in `seq 22` ; do cut -d' ' -f1-4 chr${i}.txt >> %s ; done", snpFolder, substr(sams[i], 1, 12), tmpSnpFile)) # hg19

if(file.info(tmpSnpFile)$size == 0){
  stop(sprintf("skip %s since there is no heterzygous SNP file", sam1))
}

input      = files[i]
outputTag  = sprintf("%s_asCounts_hetSNP_filterIt", sam1)
# outputTag  = sprintf("ind%s_%s_asCounts_hetSNP_filterIt",i, sam1)
input = paste(proj_input, input, sep="/")
output_prepareBAM = paste(proj_prepareBAM, outputTag, sep="/")
outputTag = paste(proj_pipeline, proj_output, outputTag, sep="/")


# There is no need to sort and filter since Paul has already done this step.
system.time(prepareBAM(input, output_prepareBAM, sortIt=FALSE, filterIt=TRUE, min.avgQ=20, min.mapQ=20, getUniqMapping=TRUE))
cat("\n", sam1, date(), "\n")

system.time(extractAsReads(paste0(output_prepareBAM, "_uniq_filtered.bam"), tmpSnpFile, outputTag))
#system.time(extractAsReads(input, tmpSnpFile, outputTag))
# Error in extractAsReads(input, tmpSnpFile, outputTag) : 
# expect a paired-read, but see something else
file.remove(tmpSnpFile)
cat("\n", sam1, date(), "\n")

sessionInfo()

q(save="no")




