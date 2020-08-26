args=(commandArgs(TRUE))
#args = c("specifications_Brain_Caudate_basal_ganglia.txt", "5e5")
args

rt.dir = getwd()
specf = args[1]
cis_window = as.numeric(args[2])

specs = unlist(read.table(specf, as.is=T))
pref = specs[1]
nsam = specs[2]
queue = specs[3]
days = specs[4]
#wrk.dir = sprintf("%s/%s", rt.dir, pref)
bas.dir = specs[16]
wrk.dir = specs[14]
setwd(wrk.dir)
lib.dir = specs[15]
mem = sprintf("%sg", specs[5])
mem

rinpdir = lib.dir
routdir = sprintf("%s/rout_%s", wrk.dir, pref)
boutdir = sprintf("%s/bout_%s", wrk.dir, pref)
if(!file.exists(routdir))dir.create(routdir)
if(!file.exists(boutdir))dir.create(boutdir)

rinp = sprintf("%s/step6_PC.R", rinpdir)
bout = sprintf("%s/step6_PC.out", boutdir)
rout = sprintf("%s/step6_PC.Rout", routdir)
com = sprintf("R360 CMD BATCH %s %s",
                        rinp, rout)
com2 = sprintf("sbatch -p %s -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                         queue, days, bout, mem, com)        

message(com2)
system(com2)





q("no")
