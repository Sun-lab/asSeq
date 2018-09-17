# ----------------------------------------------------------------------
# submit jobs
# ----------------------------------------------------------------------

for(i in 1:521){
  cat(sprintf("sbatch -t 0-1 R CMD BATCH --no-save --no-restore '--args %s'  step5_TReC_per_gene_filterIt.R ./step5_log/step5_TReC_per_gene_filterIt_%s.Rout", i,i), '\n')
}
q('no')