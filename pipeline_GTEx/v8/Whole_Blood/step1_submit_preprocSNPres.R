#cnt.dir = sprintf("/pine/scr/z/h/zhabotyn/R01/GTEx/_GTEx/data_genotype_all/cnt")
queue = "general"
days = 1
mem = "64g"
chri = 6
dir.create("rout")
dir.create("bout")
for(cis_window in c(5e5, 2e5)){
for(chri in c(1:3,5:6)){
#geneInfo = read.table(sprintf("%s/Info_ex_chr%s.dat", cnt.dir, chri), 
#                    header = T, as.is = T, sep ='\t')
  if(chri %in% c(1:12)){
    mem = "64g"
  }else{
    mem = "48g"
  }
    qout = sprintf("bout/step1_preprocSNP_%s_%s.out", chri, cis_window)
        com = sprintf("R CMD BATCH '--args %s %s' step1_preprocSNP.R rout/step1_preprocSNP_%s_%s.Rout",
                                            chri,cis_window,              chri,cis_window)
        com2 = sprintf("sbatch -p %s -N 1 -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", queue, days, qout, mem, com)        
        message(com2)
        system(com2)   
#  message(chri, " ", nrow(geneInfo))
}
}

q("no")
