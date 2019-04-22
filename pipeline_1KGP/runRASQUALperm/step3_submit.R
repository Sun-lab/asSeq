#module add tabix;module add gcc;module add gsl
#export PATH=/nas/longleaf/home/zhabotyn/progs/rasqual-master/bin/:$PATH

queue="general"
days = 7
mem = 8
for(gni in 1:2)){
  com = sprintf("R CMD BATCH '--args %s ' step3_RASQUAL.R routrb/step3_RASQUAL_%s.Rout",
                                     gni,              gni)
  qout = sprintf("boutrb/out_%s.out",
                             gni)           
  com2 = sprintf("sbatch -t 0%s-00:00:00 %s --mem=%sg --wrap=\"%s\"", days, qout, mem, com)        
  message(com2)
  system(com2)     
}