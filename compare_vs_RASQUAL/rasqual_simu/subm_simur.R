#niter = 1e1
niter=1e4
ssi = c(3:1)
b0i = c(0, .125,.25,.5, 1)
div = c(.5, 1, 2, 4, 8, Inf)
od = c(0.01, 0.10, 0.5, 2)
#b0i = c(0, .5, 1);div=c(1, Inf); od=c(0.01, 0.5)
if(!file.exists("rasq"))dir.create("rasq")
if(!file.exists("rout"))dir.create("rout")
queue = "general"

days = 7
mem = 8
for(ss in ssi){
for(b0 in b0i){
    for(dv in div){
      for(th in od){
        com = sprintf("R CMD BATCH '--args %s %s %s %s %s' simur.R rout/simur_%s_%s_%s_%s_%s.Rout",
                                              b0,th,dv,niter,ss,              b0,th,dv,niter,ss)
        qout = sprintf("--output=rout/outsh_%s_%s_%s_%s_%s.out",
                                          b0,th,dv,niter,ss)           
        com2 = sprintf("sbatch -t 0%s-00:00:00 %s --mem=%sg --wrap=\"%s\"", days, qout, mem, com)        
        message(com2)
        system(com2)
             
        }
      }
    }
}
