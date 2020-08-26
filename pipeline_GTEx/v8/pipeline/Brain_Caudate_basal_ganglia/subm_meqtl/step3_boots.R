#sbatch -p general -N 1 -t 01-00:00:00 -o batch.out --mem=32g --wragenepos_p="R CMD BATCH step1_MatrixEQTL.R step1_MatrixEQTL.Rout"
args = commandArgs(trailingOnly = TRUE)
#args = c("22", "1", "20")
#args = c("10", "1", "20")
#args = c("8", "1", "20")
#args = c("6", "43")
#args = c("2", "43")
#args = c("3", "43")
#args = c("21", "20")
#args = c("21", "50")
#need to fix chromosomes 7 and 8
args
chri = as.numeric(args[1])
numpoints = as.numeric(args[2])
        
set.seed(12345)
nsam = 280
nsub = 280
cisDist = 2e5
subs = sample(1:nsam)[1:nsub]
subse = c(0, subs)+1
subsf = c(1:9, subs+9)
subs
subsf

eigenMTdir = "/nas/longleaf/home/zhabotyn/pythlab/eigenMT"

#norm = F
#norm = T
#data.dir = "../data"
#cnt.dir = sprintf("%s/cnt", data.dir)
#geno.dir = "../datagen" 


root.dir = "/pine/scr/z/h/zhabotyn/R01"
data.dir = sprintf("%s/omnidata", root.dir)
cnt.dir = sprintf("%s/parplus_280", data.dir)
expr = read.table(sprintf("%s/gene_Tcount_chr%s.dat", cnt.dir, chri), as.is=T)

#work.dir = sprintf("%s/2019_09_23", root.dir)
#if(!file.exists(work.dir))dir.create(work.dir)
#setwd(work.dir)
library(MatrixEQTL)
useModel = modelLINEAR; 
source("helpers.R")
perm.dir = sprintf("boot%s_seed1", numpoints)
if(!file.exists(perm.dir))dir.create(perm.dir)


  blocki = 1
  countjobs = 0
  for(blocki in 1:nrow(expr)){
   results = tryCatch({
      suff0 = sprintf("%s_%s", chri, blocki)
      queue="general"
      block = 1
      days = 7
      mem = "1g"

        com = sprintf("R CMD BATCH '--args %s %s %s %s' step2_runboot.R routM/step2_runboot_%s_%s_%s_%s.Rout",
                                           chri,blocki,nsub,numpoints,                 chri,blocki,nsub,numpoints)
        qout = sprintf("boutB/step2_runboot_%s_%s_%s_%s.out",
                                       chri,blocki,nsub,numpoints) 
       filout = sprintf("%s/time_%s.csv", perm.dir, suff0)
       if(!file.exists(filout)){
          com2 = sprintf("sbatch -p %s -N 1 -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", queue, days, qout, mem, com)        
          message(com2)
          paral = 1e6
          if(blocki%%paral==0){
            system(com)
          }else{
            system(com2)
          }
          countjobs = countjobs + 1
       }
 0
 }, error = function(e){
 1
 })
 if(results==1)message(blocki)
}

#R CMD BATCH '--args 22' step1_MatrixEQTL.R step1_MatrixEQTL_22.Rout
message(countjobs)

q("no")


