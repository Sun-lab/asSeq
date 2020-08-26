#R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt' step8_submit_dynamic.R step6_submit_dynamic.Rout
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
mem = "8g"

source(sprintf("%s/helpers.R", lib.dir))

rinpdir = lib.dir
routdir = sprintf("%s/rout_%s", wrk.dir, pref)
boutdir = sprintf("%s/bout_%s", wrk.dir, pref)
if(!file.exists(routdir))dir.create(routdir)
if(!file.exists(boutdir))dir.create(boutdir)

#submit only marginal model with one covariate: genotyping PC1, genotyping PC2,
#PEER factor 1, ..., PEER factor 10
rinp = sprintf("%s/step8_dynamic.R", rinpdir)
rout = sprintf("%s/step8_dynamic.Rout", routdir)
bout = sprintf("%s/step8_dynamic.out", boutdir)
com = sprintf("R CMD BATCH %s %s",
                          rinp, rout)
com2 = sprintf("sbatch -p %s -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                       queue, days, bout, mem, com)        
message(com2)
system(com2)


q("no")
