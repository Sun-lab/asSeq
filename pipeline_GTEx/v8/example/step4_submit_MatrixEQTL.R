#R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt long 5e5' step4_submit_MatrixEQTL.R step4_submit_MatrixEQTL_long.Rout
args=(commandArgs(TRUE))
#args = c("specifications_Muscle_Skeletal.txt", "short", "5e5")
args

seedval = NA
rt.dir = getwd()
specf = args[1]
model = args[2]
cis_window = as.numeric(args[3])

specs = unlist(read.table(specf, as.is=T))
pref = specs[1]
wrk.dir = sprintf("%s/%s", rt.dir, pref)
setwd(wrk.dir)
specs = unlist(read.table("specifications.txt", as.is=T))
specs
pref = specs[1]
nsam = specs[2]
queue = specs[3]
days = specs[4]
#rprog = specs[19]
rprog = "R"
lib.dir = specs[15]#sprintf("%s/lib", rt.dir)
lib.dir
setwd(wrk.dir)
mem = "4g"
if(length(specs)>20){
  mem = specs[21]
}
mem

chrmax = 9
#days = 2
days = 5
nsub = nsam
if(is.na(seedval))seedval = specs[13]

rinpdir = lib.dir
routdir = sprintf("%s/rout_%s", wrk.dir, pref)
boutdir = sprintf("%s/bout_%s", wrk.dir, pref)
if(!file.exists(routdir))dir.create(routdir)
if(!file.exists(boutdir))dir.create(boutdir)
for(chri in 9){
  rinp = sprintf("%s/step4_MatrixEQTL.R", rinpdir)
  rout = sprintf("%s/step4_MatrixEQTL_%s_%s_%s_%s_%s_%s.Rout",
                  routdir, chri, nsam, nsub, seedval, cis_window, model) 
  bout = sprintf("%s/step4_MatrixEQTL_%s_%s_%s_%s_%s_%s.out",
                  boutdir, chri, nsam, nsub, seedval, cis_window, model)           
  com = sprintf("%s CMD BATCH \"--args %s %s %s %s %s\" %s %s",
                 rprog, chri, nsub, seedval, cis_window, model, rinp, rout)

  com2 = sprintf("sbatch -p %s -t 0%s-0:00:00 -o %s --mem=%s --wrap='%s\'", 
                            queue, days, bout, mem, com)
  message(com2)
  system(com)
}

q("no")
