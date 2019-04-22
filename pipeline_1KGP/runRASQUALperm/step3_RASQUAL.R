#module add tabix;module add gcc;module add gsl
#export PATH=/nas/longleaf/home/zhabotyn/progs/rasqual-master/bin/:$PATH
args = commandArgs(trailingOnly = TRUE)
#args = c("1")
args
gnj = as.numeric(args[1])


#root.dir = "/pine/scr/z/h/zhabotyn/R01"
#pipe.dir = sprintf("%s/omnipipe", root.dir)
#assumed structure is
#root
#root/data
#root/pipe
data.dir = "../data"
cnt.dir = sprintf("%s/cnt", data.dir)

bam.dir = sprintf("%s/bam", data.dir)
geno.dir = "../datagen" 
vcfr.dir = sprintf("%s/vcfr", geno.dir)
info.dir = "../inf"
sub.dir = "out.dir"
if(!file.exists(sub.dir))dir.create(sub.dir)


hg38 = read.table(sprintf("%s/GTF_ex.gtf", info.dir), as.is=T)
win = 2e5

res = read.csv(sprintf("%s/number_of_snps_per_gene_win%s.csv", info.dir, win),as.is=T)

vcfngz = sprintf("%s/SNP_fix.vcf.gz", vcfr.dir)

info = sprintf("%s/Info_ex.dat", cnt.dir)
totc = sprintf("%s/Tcnt_ex.dat", cnt.dir)
xcnt = sprintf("%s/samples.dat", cnt.dir)
totc = read.table(totc, as.is=T)
xcnt = read.table(xcnt, as.is=T, header=T)
infc = read.table(info, as.is=T, header=T)

totb = sprintf("%s/Tcnt_ex.bin", sub.dir)
kbin = sprintf("%s/K_ex.bin", sub.dir)
xbin = sprintf("%s/X_ex.bin", sub.dir)

if(!file.exists(totb)){
  writeBin(c(t(totc)), totb)
}
if(!file.exists(xbin)){
  writeBin(c(t(xcnt[,-1])), xbin)
}
if(!file.exists(kbin)){
  Ks = xcnt[,2]
  Ks = Ks/median(Ks)
  writeBin(Ks, kbin)
}
  
np = 4
 
  resj = res[match(infc[gnj,1], res[,2]),]
  resj
  if(resj[,3]>0){
    fli[gnj,]
    sti = infc[gnj,3]-win
    eni = infc[gnj,4]+win
  
    nm = resj[3]
    nl = resj[3]
    
    exj = hg38[infc[gnj,2] == hg38[,10],]
    stj = paste(exj[,4], collapse=",")
    enj = paste(exj[,5], collapse=",")
    genj = sprintf("%s/%s.txt", sub.dir, resj[1,1])
    timj = sprintf("%s/%s_time.txt", sub.dir, resj[1,1])
    snps = sprintf("%s:%s-%s", chr, sti, eni)

    rasBASE = sprintf("rasqual -y %s -k %s -x %s -p %s -n %s", totb, kbin, xbin, np, nsub)
    rasPOSj = sprintf("-j %s -l %s -m %s -s %s -e %s", gnj, nl, nm, stj, enj)
    rasOPT = "-z --force"
    cmdi = sprintf("tabix %s %s | %s %s -f %s %s > %s",
      vcfgz, snps, rasBASE, rasPOSj, fli[gnj, 1], rasOPT, genj)
    system(cmdi)  
  
    timsta = proc.time()
      system(cmdi)
    timend = proc.time()  
    timend-timsta
    write.table(timend[3]-timsta[3], timj, row.names=F, col.names=F, quote=F)
  }
  message(gnj)
}



q("no")
