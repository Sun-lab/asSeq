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
#PEER factor 1, ..., PEER factor 10 - overview for all together
rinp = sprintf("%s/step7_plotindp_marginal.R", rinpdir)
rout = sprintf("%s/step7_plotindp_marginal.Rout", boutdir)
bout = sprintf("%s/step7_plotindp_marginal.out", routdir)
com = sprintf("R CMD BATCH %s %s", rinp, rout)
com2 = sprintf("sbatch -p %s -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                         queue, days, bout, mem, com)        
message(com2)
#system(com2)


#submit condition of interest (age, CTCF or TP53) along with 2 genotyping PCs
for(cndi in 1:3){
  rinp = sprintf("%s/step7_plotindp_gPC.R", rinpdir)
  bout = sprintf("%s/step7_plotindp_gPC.out", boutdir)
  rout = sprintf("%s/step7_plotindp_gPC.Rout", routdir)
  com = sprintf("R CMD BATCH '--args %s' %s %s",
                                   cndi, rinp, rout)
  com2 = sprintf("sbatch -p %s -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                        queue, days, bout, mem, com)        
  if(!file.exists(bout)){
    message(com2)
    #system(com2)
  }
}

#submit condition of interest (age, CTCF or TP53) along with 2 genotyping PCs and 5 PEER factors
for(cndi in 1:3){
  rinp = sprintf("%s/step7_plotindp_gPC_PF.R", rinpdir)
  bout = sprintf("%s/step7_plotindp_gPC_PF.out", boutdir, cndi)
  rout = sprintf("%s/step7_plotindp_gPC_PF.Rout", routdir, cndi)
  com = sprintf("R CMD BATCH '--args %s' %s %s",
                                   cndi, rinp, rout)
  com2 = sprintf("sbatch -p %s -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                     queue, days, bout, mem, com)        
  if(!file.exists(bout)){
    message(com2)
    #system(com2)
  }
}


q("no")
