# sbatch R CMD BATCH --no-save --no-restore step11_submit_step11.R
for(i in (1:22)){
  com<-sprintf("sbatch --exclusive R CMD BATCH --no-save --no-restore '--args %s' step11_get_snp.R step11_get_snp_%s.Rout",i,i)
  message(com)
  system(com)
}

q("no")
