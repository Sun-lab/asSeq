#R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt 5e5' step5_submit_collect.R step5_submit_collect.Rout
args=(commandArgs(TRUE))
#args = c("specifications_Adipose_Visceral_Omentum.txt", "5e5")
args

rt.dir = getwd()
specf = args[1]
cis_window = as.numeric(args[2])

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

nsub = nsam

rinpdir = lib.dir
routdir = sprintf("%s/rout_%s", wrk.dir, pref)
boutdir = sprintf("%s/bout_%s", wrk.dir, pref)
if(!file.exists(routdir))dir.create(routdir)
if(!file.exists(boutdir))dir.create(boutdir)

rinp = sprintf("%s/step5_collectU.R", rinpdir)
bout = sprintf("%s/step5_collectU_%s.R", 
                boutdir, cis_window)
rout = sprintf("%s/step5_collectU_%s.R", 
                routdir, cis_window)
com = sprintf("R CMD BATCH '--args %s' %s %s",
                           cis_window, rinp, rout)
com2 = sprintf("sbatch -p %s -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                         queue, days, bout, mem, com)        
message(com2)
#system(com2)

q("no")
