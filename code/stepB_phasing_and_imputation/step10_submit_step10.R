# sbatch R CMD BATCH --no-save --no-restore step10_submit_step10.R
for(i in (1:22)){
  com<-sprintf("sbatch R CMD BATCH --no-save --no-restore '--args %s' step10_get_most_likely_genotype.R step10_get_most_likely_genotype_%s.Rout",i,i)
  message(com)
  system(com)
}

q("no")
