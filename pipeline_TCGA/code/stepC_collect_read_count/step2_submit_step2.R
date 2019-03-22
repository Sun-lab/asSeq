# -----------------------------------------------------------------
# submit jobs
# -----------------------------------------------------------------

setwd("/fh/fast/sun_w/licai/_tumor_eQTL/R_batch1/stepC_collect_read_count/pipeline")
for(i in 1:121){
  cat(sprintf("sbatch --mem=32000 R CMD BATCH --no-save --no-restore '--args %02d' step2_sort_bam.R ./step2_log/step2_sort_bam_%02d.Rout", i,i), '\n')
}
q('no')