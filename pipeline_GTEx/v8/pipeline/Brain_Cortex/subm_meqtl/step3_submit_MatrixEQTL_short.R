#root.dir = "/pine/scr/z/h/zhabotyn/R01"
#/pine/scr/z/h/zhabotyn/R01/GTEx/_GTEx/data_genotype_all/genotype
#root.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/_GTEx/"
#work.dir = sprintf("%s/2019_05_01", root.dir)
#setwd(work.dir)
#data.dir = sprintf("%s/omnidata", root.dir)
#cnt.dir = sprintf("%s/parplus_280", data.dir)


queue="general"
block = 1
chrmax = 22
days = 1
nsam = 205
nsub = 205
seedval = 1565691
cis_window = 5e5
model = "short"
chri = 22

boutdir = sprintf("boutM%s", nsam)
routdir = sprintf("routM%s", nsam)
if(!file.exists(boutdir))dir.create(boutdir)
if(!file.exists(routdir))dir.create(routdir)
for(chri in 1:chrmax){
#chri = 22
    if(chri>8){
      mem = "1g"
    }else{
      mem = "1g"
    }
    #expr = read.table(sprintf("%s/gene_Tcount_chr%s.dat", cnt.dir, chri), as.is=T)
    #info = read.table(sprintf("%s/Info_ex_chr%s.dat", cnt.dir, chri), header=T, as.is=T)  
    #ngen = nrow(info)

        suff = sprintf("%s", chri)
        com = sprintf("R CMD BATCH '--args %s %s %s %s %s %s' step3_MatrixEQTL.R %s/step3_MatrixEQTL_%s_%s_%s_%s_%s_%s.Rout",
                                           chri, nsam, nsub, seedval, cis_window, model,
                                   routdir,chri, nsam, nsub, seedval, cis_window, model)
#args = c("22", "194", "194", "1565691", "5e5", "long", "1")
        qout = sprintf("%s/step3_MatrixEQTL_%s_%s_%s_%s_%s_%s.out",
                                    boutdir, chri, nsam, nsub, seedval, cis_window, model)           
          com2 = sprintf("sbatch -p %s -N 1 -t 0%s-12:00:00 -o %s --mem=%s --wrap=\"%s\"", queue, days, qout, mem, com)        
          message(com2)
          system(com2)
        #}     
}

q("no")
