#R CMD BATCH '--args 21' step0_get_het_snps.R rout/step0_get_het_snps_21.Rout"

gtex.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/_GTEx"

queue = "general"
days = 7
chri = 6
mem = "1g"
nsub = 352
cnt = 0

routdir = sprintf("rout")
boutdir = sprintf("bout")
if(!file.exists(routdir))dir.create(routdir)
if(!file.exists(boutdir))dir.create(boutdir)
for(chri in 1:22){
  if(chri %in% 1:10){
    mem = "24g"
  }else{
    mem = "16g"
  }
  qout = sprintf("%s/step0_get_het_snps_%s.out", boutdir, chri)
  com = sprintf("R CMD BATCH '--args %s' step0_get_het_snps.R %s/step0_get_het_snps_%s.Rout",
                                     chri,                   routdir,         chri)
  com2 = sprintf("sbatch -p %s -N 1 -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", queue, days, qout, mem, com)        

  message(com2)
  system(com2)
}

q("no")
