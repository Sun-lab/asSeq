# ----------------------------------------------------------------------
# submit jobs
# ----------------------------------------------------------------------

setwd("/fh/fast/sun_w/licai/_tumor_eQTL/R_batch1/stepC_collect_read_count/pipeline")
for(i in 1:121){
  cat(sprintf("sbatch --mem=32000 R CMD BATCH --no-save --no-restore '--args %s' step3_extract_asReads_filterIt.R ./step3_log/step3_extract_asReads_filterIt_%s.Rout", i,i), '\n')
}
q('no')