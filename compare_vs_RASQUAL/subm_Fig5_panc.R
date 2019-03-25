#if you want to run multiple b0's:
#about 8.5 seconds per iteration for all profiles
#expect 71 minutes per 500 iterations
#if you run only b0=0 - under 15 minutes
workdir = "/pine/scr/z/h/zhabotyn/R01/2019_03_20"
setwd(workdir)

niter = 5e2
queue = "general"
mem = 8000
ss = 2
#b0s = c(0, .125,.25,.5, 1)
b0s = c(.5, 1)
#div = c(.5, 1, 2, 4, 8, Inf)
div = 1#a subset for figure (a)
ods = c(0.01, 0.10, 0.5, 2)
#if(!file.exists("rasq"))dir.create("rasq")
if(!file.exists("checkrout"))dir.create("checkrout")

days = 7
mem = 8
b1 = 0
nsnp = c(4)
st = proc.time()
for(nsnpi in nsnp){
  for(ssi in ss){
  for(b0i in b0s){
      for(dvi in div){
        for(odi in ods){
          rout = sprintf("checkrout/simush_%s_%s_%s_%s_%s_%s.Rout",
                          b0i, odi, dvi, niter, ssi, nsnpi) 
          com = sprintf("R CMD BATCH '--args %s %s %s %s %s %s' simush.R %s",
                                        b0i, odi, dvi, niter, ssi, nsnpi, rout)              
          qout = sprintf("--output=checkrout/outsh_%s_%s_%s_%s_%s_%s.out",
                                        b0i, odi, dvi, niter, ssi, nsnpi)           
  #        com2 = sprintf("bsub -q %s -M 8 \"%s\"", queue, com)        
          com2 = sprintf("sbatch -t 0%s-00:00:00 %s --mem=%sg --wrap=\"%s\"", 
                        days, qout, mem, com)        
          message(com)
          system(com)
          message(ss, " ", b0, " ", th)      
          }
        }
      }
  }
}
en = proc.time()
en[3]-st[3]
