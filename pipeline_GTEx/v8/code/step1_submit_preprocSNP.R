#R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt' step1_submit_preprocSNP.R step1_submit_preprocSNP.Rout
args=(commandArgs(TRUE))
#args = "specifications_Adipose_Visceral_Omentum.txt"
args

rt.dir = getwd()
specf = args[1]

specs = unlist(read.table(specf, as.is=T))
pref = specs[1]
wrk.dir = sprintf("%s/%s", rt.dir, pref)
setwd(wrk.dir)
specs = unlist(read.table("specifications.txt", as.is=T))
nsam = specs[2]
queue = specs[3]
days = specs[4]
bmem = as.numeric(specs[5])
rcmd = specs[19]
wrk.dir = sprintf("%s/%s", rt.dir, pref)
lib.dir = sprintf("%s/lib", rt.dir)
setwd(wrk.dir)

mem0 = sprintf("%sg", round(bmem/3))
mem1 = sprintf("%sg", bmem)
mem2 = sprintf("%sg", round(bmem*1.5))
c(mem1, mem2)

rinpdir = lib.dir
routdir = sprintf("%s/rout_%s", wrk.dir, pref)
boutdir = sprintf("%s/bout_%s", wrk.dir, pref)
if(!file.exists(routdir))dir.create(routdir)
if(!file.exists(boutdir))dir.create(boutdir)

for(cis_window in c(5e5)){
  for(chri in 1:22){
    mem = mem1
    if(chri %in% c(1:12))mem = mem2
  
    rinp = sprintf("%s/step1_preprocSNP.R", rinpdir)
    bout = sprintf("%s/step1_preprocSNP_%s_%s.out", boutdir, chri, cis_window)
    rout = sprintf("%s/step1_preprocSNP_%s_%s.Rout", routdir, chri, cis_window)
   
    com = sprintf("%s CMD BATCH '--args %s %s' %s %s",
                                       rcmd, chri,cis_window, rinp, rout)
    com2 = sprintf("sbatch -p %s -N 1 -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                            queue, days, bout, mem, com)        
    message(com2)
    system(com2)   
  }
}

q("no")
