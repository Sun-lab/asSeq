#or you can run all combinations in a few hours
niter = 5e2
#niter = 2e0
queue = "general"
ss = 2
b0s = c(0, .125,.25,.5, 1)
#b0s = 0
#div = 1#a subset for figure (b)
div = c(.5, 1, 2, 4, 8, Inf)
ods = c(0.01, 0.10, 0.5, 2)
if(!file.exists("checkrout"))dir.create("checkrout")

days = 1
mem = 4
b1 = 0
st = proc.time()
for(ssi in ss){
  for(b0i in b0s){
      for(dvi in div){
        for(odi in ods){
          rout = sprintf("checkrout/simush_%s_%s_%s_%s_%s.Rout",
                          b0i, odi, dvi, niter, ssi) 
          com = sprintf("R CMD BATCH '--args %s %s %s %s %s' simush.R %s",
                                        b0i, odi, dvi, niter, ssi, rout)              
          qout = sprintf("--output=checkrout/outsh_%s_%s_%s_%s_%s.out",
                                        b0i, odi, dvi, niter, ssi)           
  #        com2 = sprintf("bsub -q %s -M 8 \"%s\"", queue, com)        
          com2 = sprintf("sbatch -t 0%s-00:00:00 %s --mem=%sg --wrap=\"%s\"", 
                        days, qout, mem, com)        
          message(com)
          #system(com)
          system(com2)
          message(ssi, " ", b0i, " ", dvi, " ", odi, " ")      
          }
        }
      }
}
en = proc.time()
en[3]-st[3]
