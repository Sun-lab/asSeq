#module add tabix;module add gcc;module add gsl
#export PATH=/nas/longleaf/home/zhabotyn/progs/rasqual-master/bin/:$PATH
args = commandArgs(trailingOnly = TRUE)
#args = c("11", "1", "5", "280")
#args = c("22", "151", "160", "280")
#args = c("11", "151", "160", "280")
args
chri = as.numeric(args[1])
stgn = as.numeric(args[2])
engn = as.numeric(args[3])
nsub = as.numeric(args[4])


#root.dir = "/pine/scr/z/h/zhabotyn/R01"
#pipe.dir = sprintf("%s/omnipipe", root.dir)
#assumed structure is
#root
#root/data
#root/pipe
data.dir = "../omnidata"
vcf.dir = "../data_genotype/465ind/rvcf"
tre.dir = sprintf("%s/parplus_280",data.dir)
ras.dir = sprintf("%s/rasqtmpupd",data.dir)
#output directory
if(!file.exists(ras.dir))dir.create(ras.dir)
if(nsub==280){
  sub.dir = sprintf("%s/all", ras.dir)
  out.dir = sprintf("%s/all", ras.dir)
}else{
  sub.dir = sprintf("%s/subnew%s", ras.dir, nsub)
  out.dir = sprintf("%s/subnew%s", ras.dir, nsub)
}
if(!file.exists(out.dir))dir.create(out.dir)

chr = sprintf("chr%s", chri)
fli = sprintf("%s/gene_info_%s.dat",tre.dir,chr)
fli = read.table(fli,header=T)
hg38 = read.table(sprintf("%s/exons_%s.gtf", sub.dir, chr), as.is=T)

res = read.csv("number_of_snps_per_gene.csv",as.is=T)
if(nsub==280){
  vcfn = sprintf("%s/all_hap_%s.vcf", sub.dir, chr)
  totb = sprintf("%s/all_tot_%s.bin", sub.dir, chr)
  kbin = sprintf("%s/all_K_%s.bin", sub.dir, chr)
  xbin = sprintf("%s/all_X.bin", sub.dir)
}else{
  vcfn = sprintf("%s/sub_hap_%s.vcf", sub.dir, chr)
  totb = sprintf("%s/sub_tot_%s.bin", sub.dir, chr)
  kbin = sprintf("%s/sub_K_%s.bin", sub.dir, chr)
  xbin = sprintf("%s/sub_X.bin", sub.dir)
}
vcfgz = sprintf("%s.gz", vcfn)
np = 4
 
win = 2e5

for(gnj in stgn:engn){
#gnj = 1
#gnj = gnj + 1
#  gnj = grep("ENSG00000149806.8", fli[,1]);gnj;fli[gnj,]
  resj = res[match(fli[gnj,1], res[,1]),]
  resj
  if(resj$within>=0){
    fli[gnj,]
    sti = fli[gnj,3]-win
    eni = fli[gnj,4]+win
  
    nm = resj[4]
    nl = resj[4]
    
    exj = hg38[fli[gnj,1] == hg38[,10],]
    stj = paste(exj[,4], collapse=",")
    enj = paste(exj[,5], collapse=",")
    genj = sprintf("%s/%s.txt", out.dir, resj[1,1])
    timj = sprintf("%s/%s_time.txt", out.dir, resj[1,1])
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
