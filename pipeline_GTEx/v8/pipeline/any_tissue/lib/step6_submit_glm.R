queue = "general"
mem = "4g"
days = 7
routdir = "rout"
boutdir = "bout"
if(!file.exitst(routdir))dir.create(routdir)
if(!file.exitst(boutdir))dir.create(boutdir)
#submit only marginal nmodel with one covariate: genotyping PC1, genotyping PC2,
#PEER factor 1, ..., PEER factor 10
for(cnd2i in 1:12){ # step7_trecase_permute # 
  qout = sprintf("%s/step6_glm_marginal_%s.out", routdir, cnd2i)
  rout = sprintf("%s/step6_glm_marginal_%s.Rout", boutdir, cnd2i)
  com = sprintf("R CMD BATCH '--args %s' step6_glm_marginal.R %s",
                                   cnd2i,        rout)
  com2 = sprintf("sbatch -p %s -N 1 -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                         queue, days, qout, mem, com)        
  if(!file.exists(qout)){
     message(com2)
     system(com2)
  }
}

#submit condition of interest (age, CTCF or TP53) along with 2 genotyping PCs
for(cndi in 1:3){
  qout = sprintf("%s/step6_glm_gPC_%s.out", boutdir, cndi)
  rout = sprintf("%s/step6_glm_gPC_%s.Rout", routdir, cndi)
  com = sprintf("R CMD BATCH '--args %s' step6_glm_gPC.R %s",
                                   cndi,        rout)
  com2 = sprintf("sbatch -p %s -N 1 -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                        queue, days, qout, mem, com)        
  if(!file.exists(qout)){
    message(com2)
    system(com2)
  }
}

#submit condition of interest (age, CTCF or TP53) along with 2 genotyping PCs and 5 PEER factors
for(cndi in 1:3){
  qout = sprintf("%s/step6_glm_gPC_PF_%s.out", boutdir, cndi)
  rout = sprintf("%s/step6_glm_gPC_PF_%s.Rout", routdir, cndi)
  com = sprintf("R CMD BATCH '--args %s' step6_glm_gPC_PF.R %s",
                                   cndi,        rout)
  com2 = sprintf("sbatch -p %s -N 1 -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                     queue, days, qout, mem, com)        
  if(!file.exists(qout)){
    message(com2)
    system(com2)
  }
}


q("no")
