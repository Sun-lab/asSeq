#module add tabix;module add gcc;module add gsl;module add r
#export PATH=/nas/longleaf/home/zhabotyn/progs/rasqual-master/bin/:$PATH
args = commandArgs(trailingOnly = TRUE)
#args = c("1")
#args = c("2")
args
gnj = as.numeric(args[1])

#to run the model on permuted data change this flag to TRUE
permute=FALSE


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
info.dir = "../inf"
if(permute){
  vcfr.dir = sprintf("%s/vcfp", geno.dir)
  sub.dir = "per.dir"
}else{
  vcfr.dir = sprintf("%s/vcfr", geno.dir)
  sub.dir = "out.dir"
}
if(!file.exists(sub.dir))dir.create(sub.dir)


hg38 = read.table(sprintf("%s/GTF_ex.gtf", info.dir), as.is=T)
win = 2e5

res = read.csv(sprintf("%s/number_of_snps_per_gene_win%s.csv", info.dir, win),as.is=T)

vcfgz = sprintf("%s/SNP_fix.vcf.gz", vcfr.dir)

info = sprintf("%s/Info_ex.dat", cnt.dir)
totc = sprintf("%s/Tcnt_ex.dat", cnt.dir)
xcnt = sprintf("%s/samples.dat", cnt.dir)
totc = read.table(totc, as.is=T)
xcnt = read.table(xcnt, as.is=T, header=T)
infc = read.table(info, as.is=T, header=T)

totb = sprintf("%s/Tcnt_ex.bin", cnt.dir)
kbin = sprintf("%s/K_ex.bin", cnt.dir)
xbin = sprintf("%s/X_ex.bin", cnt.dir)

#
#  kmatr = 10^(X$ttlCount)
#  kmatr = matrix(rep(kmatr/mean(kmatr),each=nrow(fli)), nrow=nrow(fli), ncol=nn)
#  kbin = sprintf("%s/all_K_chr%s.bin", sub.dir, i)
#  writeBin(as.double(c(t(kmatr))), con=file(kbin, "wb"), double())
#  writeBin(as.double(c(as.matrix(X[,c(2:5)]))), con=file(xbin,"wb"), double())
#  writeBin(as.double(t(as.matrix(tti))), con=file(totb, "wb"), double())


#these should alredy exist, but if not, create them
if(!file.exists(totb)){
  con = file(totb, "wb")
  writeBin(as.double(t(totc)), con=con, double())
  close(con)
}
if(!file.exists(xbin)){
  con = file(xbin, "wb")
  writeBin(as.double(c(as.matrix(xcnt[,c(2:5)]))), con=con, double())
  close(con)
}
if(!file.exists(kbin)){
  Ks = 10^xcnt[,2]
  Ks = matrix(rep(Ks/mean(Ks),each=nrow(totc)), nrow=nrow(totc))
  con = file(kbin,"wb")
  writeBin(as.double(c(t(Ks))), con=con, double())
  close(con)
}
  
#number of variables
np = 4
#number of samples
nsam = 25
 
resj = res[match(infc[gnj,1], res[,2]),]
resj
if(resj[,3]>0){
  chr = infc[gnj,2]
  sti = infc[gnj,3]-win
  eni = infc[gnj,4]+win
  
  nm = resj[1,3]
  nl = resj[1,3]
  kpe =  hg38[,11] %in% infc[gnj,1]
  exj = hg38[kpe,]
  stj = paste(exj[,5], collapse=",")
  enj = paste(exj[,6], collapse=",")
  genj = sprintf("%s/%s.txt", sub.dir, resj[1,2])
  timj = sprintf("%s/%s_time.txt", sub.dir, resj[1,2])
  snps = sprintf("%s:%s-%s", chr, sti, eni)

  rasBASE = sprintf("rasqual -y %s -k %s -x %s -p %s -n %s", totb, kbin, xbin, np, nsam)
  rasPOSj = sprintf("-j %s -l %s -m %s -s %s -e %s", gnj, nl, nm, stj, enj)
  rasOPT = "-z --force --posterior-genotype"# --population-only --min-coverage-depth 0.0 --as-only 
  cmdi = sprintf("tabix %s %s | %s %s -f %s %s > %s",
    vcfgz, snps, rasBASE, rasPOSj, infc[gnj, 1], rasOPT, genj)
  
  timsta = proc.time()
    system(cmdi)
  timend = proc.time()  
  timend-timsta
  write.table(timend[3]-timsta[3], timj, row.names=F, col.names=F, quote=F)
}
message(gnj)

if(file.size(genj)>100){
  resj = read.table(genj)
  nms = c("fID", "sID", "SNPchr", "SNPpos", "Ref", "Alt", "AF", "HWEChisq", 
          "impQual", "qval", "Chi2", "Pi", "DeltaErr", "PhiBias", "OD", "SNPID",
          "NfSNP", "NtestSNP", "niterNull", "niterAlt", "rndtieloc", "llh0",
          "conv0", "corfSNP", "corrSNP")
  colnames(resj) = nms
  resj[order(resj$qval),][1:4,]
  summary(10^resj$qval)
}

q("no")
