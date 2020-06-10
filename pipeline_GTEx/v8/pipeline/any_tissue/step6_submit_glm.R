#R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt' step6_submit_glm.R step6_submit_glm.Rout
args=(commandArgs(TRUE))
#args = c("specifications_Adipose_Visceral_Omentum.txt")
args

rt.dir = getwd()
specf = args[1]

specs = unlist(read.table(specf, as.is=T))
pref = specs[1]
wrk.dir = sprintf("%s/%s", rt.dir, pref)
setwd(wrk.dir)
specs = unlist(read.table("specifications.txt", as.is=T))
specs
nsam = specs[2]
queue = specs[3]
days = specs[4]
lib.dir = sprintf("%s/lib", rt.dir)
mem = "2g"

source(sprintf("%s/helpers.R", lib.dir))

rinpdir = lib.dir
routdir = sprintf("%s/rout_%s", wrk.dir, pref)
boutdir = sprintf("%s/bout_%s", wrk.dir, pref)
if(!file.exists(routdir))dir.create(routdir)
if(!file.exists(boutdir))dir.create(boutdir)

#submit only marginal model with one covariate: genotyping PC1, genotyping PC2,
#PEER factor 1, ..., PEER factor 10
for(cnd2i in 1:12){
  rinp = sprintf("%s/step6_glm_marginal.R", rinpdir)
  rout = sprintf("%s/step6_glm_marginal_%s.Rout", boutdir, cnd2i)
  bout = sprintf("%s/step6_glm_marginal_%s.out", routdir, cnd2i)
  com = sprintf("R CMD BATCH '--args %s' %s %s",
                                   cnd2i, rinp, rout)
  com2 = sprintf("sbatch -p %s -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                         queue, days, bout, mem, com)        
  if(!file.exists(bout)){
     message(com2)
     #system(com2)
  }
}

#submit condition of interest (age, CTCF or TP53) along with 2 genotyping PCs
for(cndi in 1:3){
  rinp = sprintf("%s/step6_glm_gPC.R", rinpdir)
  bout = sprintf("%s/step6_glm_gPC_%s.out", boutdir, cndi)
  rout = sprintf("%s/step6_glm_gPC_%s.Rout", routdir, cndi)
  com = sprintf("R CMD BATCH '--args %s' %s %s",
                                   cnd2i, rinp, rout)
  com2 = sprintf("sbatch -p %s -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                        queue, days, bout, mem, com)        
  if(!file.exists(bout)){
    message(com2)
    #system(com2)
  }
}

#submit condition of interest (age, CTCF or TP53) along with 2 genotyping PCs and 5 PEER factors
for(cndi in 1:3){
  rinp = sprintf("%s/step6_glm_gPC_PF.R", rinpdir)
  bout = sprintf("%s/step6_glm_gPC_PF_%s.out", boutdir, cndi)
  rout = sprintf("%s/step6_glm_gPC_PF_%s.Rout", routdir, cndi)
  com = sprintf("R CMD BATCH '--args %s' %s %s",
                                   cnd2i, rinp, rout)
  com2 = sprintf("sbatch -p %s -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                     queue, days, bout, mem, com)        
  if(!file.exists(bout)){
    message(com2)
    #system(com2)
  }
}


q("no")
