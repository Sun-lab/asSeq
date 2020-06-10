#R CMD BATCH  '--args specifications_Adipose_Visceral_Omentum.txt long 5e5' step2_submit_trecaseA.R step2_submit_trecaseA_long.Rout
args=(commandArgs(TRUE))
#args = c("specifications_Adipose_Visceral_Omentum.txt", "long", "5e5")
args

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
lib.dir = sprintf("%s/lib", rt.dir)
mem = "1g"

source(sprintf("%s/helpers.R", lib.dir))

nsub = nsam
int.dir = sprintf("%s_%s_%s", pref, nsub, cis_window)
fls = list.files(int.dir, "counti")
fls = fls[grep(model, fls)]

gns = gsub(".csv", "", t(sapply(fls, get_blockn, blocks=3:4)))
gns = aggregate(rep(1, nrow(gns)), by=list(gns[,1]), FUN=sum)
gns

rinpdir = lib.dir
routdir = sprintf("%s/rout_%s", wrk.dir, pref)
boutdir = sprintf("%s/bout_%s", wrk.dir, pref)
if(!file.exists(routdir))dir.create(routdir)
if(!file.exists(boutdir))dir.create(boutdir)

cnt = 0
for(i in 1:nrow(gns)){
  chri = gns[i, 1]
  for(indi in 1:gns[i,2]){ # step7_trecase_permute # 
    rinp = sprintf("%s/step2_trecase1.R", rinpdir)
    bout = sprintf("%s/step2_trecase1_%s_%s_%s_%s_%s.out", 
                    boutdir, chri, indi, nsub, cis_window, model)
    rout = sprintf("%s/step2_trecase1_%s_%s_%s_%s_%s.Rout", 
                    routdir, chri, indi, nsub, cis_window, model)
    com = sprintf("R CMD BATCH '--args %s %s %s %s' %s %s",
                               chri,indi, cis_window, model, rinp, rout)
    com2 = sprintf("sbatch -p %s -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", 
                             queue, days, bout, mem, com)        
    if(!file.exists(bout)){
      message(com2)
      #system(com2)
      cnt = cnt + 1
    }   
  }
  message(chri, " ", indi)
}
cnt

q("no")
