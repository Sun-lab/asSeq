#module add tabix;module add gcc;module add gsl


pipe.dir = sprintf("%s/pipe", root.dir)
wrk.dir = sprintf("%s/2018_10_30", pipe.dir)
setwd(wrk.dir)
if(!file.exists("routrb"))dir.create("routrb")
if(!file.exists("boutrb"))dir.create("boutrb")
res = read.csv(sprintf("%s/number_of_snps_per_gene.csv", wrk.dir),as.is=T)
tbl = table(res[,5])
nms = gsub("chr", "", names(tbl))
o = order(as.numeric(nms))
tbl = tbl[o]
nms[o]
tbl
len = tbl


queue="general"
block = 1
chrmax = 22
days = 7
mem = 8
for(subs in c(280, 200, 100, 50, 25)){
  for(chri in chrmax:1){
    for(blockj in 1:ceiling(len[chri]/block)){
    stj = (blockj-1)*block+1
    enj = blockj*block
    if(enj>len[chri])enj=len[chri]    
        com = sprintf("R CMD BATCH '--args %s %s %s %s' step3_RASQUAL.R routrb/step3_RASQUAL_%s_%s_%s_%s.Rout",
                                              chri,stj,enj,subs,              chri,stj,enj,subs)
        qout = sprintf("boutrb/out_%s_%s_%s_%s.out",
                                 chri,stj,enj,subs)           
        com2 = sprintf("sbatch -t 0%s-00:00:00 %s --mem=%sg --wrap=\"%s\"", days, qout, mem, com)        
        message(com2)
        system(com2)     
      }
  }
}