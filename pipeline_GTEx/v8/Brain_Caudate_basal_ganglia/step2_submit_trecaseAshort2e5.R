#root.dir = "/pine/scr/z/h/zhabotyn/R01"            cis_window
#gtex.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/_GTEx"
#cnt.dir = sprintf("%s/data_genotype_all/cnt", gtex.dir)

model = "short"
pref = "Brain_Caudate_basal_ganglia"
cis_window  = 2e5
queue = "general"
days = 7
chri = 6
dir.create("rout")
dir.create("bout")
mem = "1g"
nsub = 194
#nsub = 668
cnt = 0

get_blocks = function(x, split="_", blocks=2:3){
  unlist(strsplit(x, split=split))[blocks]
}

int.dir = sprintf("%s_%s_%s", pref, nsub, cis_window)
fls = list.files(int.dir, "counti")
fls = fls[grep(model, fls)]

gns = gsub(".csv", "", t(sapply(fls, get_blocks, blocks=3:4)))
gns = aggregate(rep(1, nrow(gns)), by=list(gns[,1]), FUN=sum)
gns[1:5,]

for(i in 1:21){
  chri = gns[i, 1]
  for(indi in 1:gns[i,2]){ # step7_trecase_permute # 
    qout = sprintf("bout/step2_trecase1_%s_%s_%s_%s_%s.out", chri, indi, nsub, cis_window, model)
        com = sprintf("R CMD BATCH '--args %s %s %s %s %s' step2_trecase1.R rout/step2_trecase1_%s_%s_%s_%s_%s.Rout",
                                            chri,indi,nsub, cis_window, model,              chri,indi,nsub, cis_window, model)
        com2 = sprintf("sbatch -p %s -N 1 -t 0%s-00:00:00 -o %s --mem=%s --wrap=\"%s\"", queue, days, qout, mem, com)        
        if(!file.exists(qout)){
          message(com2)
          system(com2)
          cnt = cnt + 1
        }   
  }
  message(chri, " ", indi)
}
cnt

q("no")
