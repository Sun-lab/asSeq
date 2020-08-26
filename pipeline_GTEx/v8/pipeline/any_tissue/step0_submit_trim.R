#usage
#R CMD BATCH '--args spec_file_name long' step0_submit_trim.R step0_submit_trim_long_AVO.Rout
#for example:
#R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt long' step0_submit_trim.R step0_submit_trim_long_AVO.Rout
args=(commandArgs(TRUE))
#args = "specifications_Adipose_Visceral_Omentum.txt"
args

seedval = 1565691
rt.dir = getwd()
specf = args[1]
model = args[2]
specs = unlist(read.table(specf, as.is=T))
specs
pref = specs[1]
nsam = specs[2]
queue = specs[3]
days = specs[4]
bmem = as.numeric(specs[5])
wrk.dir = sprintf("%s/%s", rt.dir, pref)
lib.dir = sprintf("%s/lib", rt.dir)
if(!file.exists(wrk.dir))dir.create(wrk.dir)
setwd(wrk.dir)
if(!is.na(seedval))specs[13]=seedval
specs[14] = wrk.dir
specs[15] = lib.dir
write.table(specs, "specifications.txt", row.names=F, col.names=F)

mem0 = sprintf("%sg", round(bmem/3))
mem1 = sprintf("%sg", bmem)
mem2 = sprintf("%sg", round(bmem*1.5))
c(mem1, mem2)

rinpdir = lib.dir
routdir = sprintf("%s/rout_%s", wrk.dir, pref)
boutdir = sprintf("%s/bout_%s", wrk.dir, pref)
if(!file.exists(routdir))dir.create(routdir)
if(!file.exists(boutdir))dir.create(boutdir)

mem = mem0
rinp = sprintf("%s/step0_trim_%s.R", rinpdir, model)
bout = sprintf("%s/step0_trim_%s.out", boutdir, model)
rout = sprintf("%s/step0_trim_%s.Rout", routdir, model)
com = sprintf("R CMD BATCH %s %s",
                          rinp, rout)
com2 = sprintf("sbatch -p %s -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                       queue, days, bout, mem, com)        
message(com2)
#system(com2)   

q("no")
