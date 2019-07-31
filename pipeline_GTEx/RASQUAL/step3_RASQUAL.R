# ml tabix;
# export CFLAGS="-I/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL/CLAPACK-3.2.1/INCLUDE -I/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL/CLAPACK-3.2.1/F2CLIBS -I/usr/include/gsl"
# export LDFLAGS="-L/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL/CLAPACK-3.2.1 -L/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL/CLAPACK-3.2.1/F2CLIBS -I/usr/include/gsl/lib"
# export PATH=/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL/rasqual/bin/:$PATH
setwd("/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL")
args = commandArgs(trailingOnly = TRUE)
#args = c("1", "22", "T")
#args = c("1201", "1", "F")
options(digits=22)
options(max.print=100000) 
args
gnj = as.numeric(args[1])
chri = as.numeric(args[2])

#to run the model on permuted data change this flag to TRUE
permute=as.logical(args[3])

#root.dir = "/pine/scr/z/h/zhabotyn/R01"
#pipe.dir = sprintf("%s/omnipipe", root.dir)
#assumed structure is
#root
#root/data
#root/pipe
data.dir = "/fh/scratch/delete90/sun_w/_GTEx/ncbi"
bam.dir = sprintf("%s/wb_raw_bamfiles", data.dir)
geno.dir = "../data_genotype_all" 
cnt.dir = sprintf("%s/cnt", geno.dir)
info.dir = ".."
if(permute){
  vcfr.dir = sprintf("%s/vcfp", geno.dir)
  sub.dir = "per.dir/detailed"
}else{
  vcfr.dir = sprintf("%s/vcfr", geno.dir)
  sub.dir = "out.dir/detailed"
}
if(!file.exists(sub.dir))dir.create(sub.dir)

get_block = function(str, split="\\.bam",block=1){
  unlist(strsplit(str, split=split))[block]
}

hg19 = read.table(sprintf("%s/gencode.v24lift37.annotation_nchr.gtf", info.dir), 
                  as.is=T, sep = '\t') 
hg19 = hg19[which(hg19$V1 == chri), ]
hg19[,10] = sapply(hg19[,9], get_block, split=";", block=1)
hg19[,11] = sapply(hg19[,10], get_block, split=" ", block=2)
hg19[1:2,]

win = 1e5

res = read.csv(sprintf("number_of_snps_per_gene_win1e+05_chr%s.csv", chri),as.is=T)

vcfgz = sprintf("%s/SNP_chr%s.vcf.gz", vcfr.dir, chri)

info = sprintf("%s/Info_ex_chr%s.dat", cnt.dir, chri)
# totc = sprintf("%s/Tcnt_ex_chr%s.dat", cnt.dir, chri)
# xcnt = sprintf("%s/samples.dat", cnt.dir)
# totc = read.table(totc, as.is=T)
# xcnt = read.table(xcnt, as.is=T, header=T)
infc = read.table(info, as.is=T, header=T)

totb = sprintf("%s/Tcnt_ex_chr%s.bin", cnt.dir, chri)
kbin = sprintf("%s/K_ex_chr%s.bin", cnt.dir, chri)
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
np = 1
#number of samples
nsam = 355
 
resj = res[match(infc[gnj,1], res[,2]),]
resj
if(resj[,3]>0){
  chr = infc[gnj,2]
  sti = infc[gnj,3]-win
  eni = infc[gnj,4]+win
  
  nm = resj[1,3]
  nl = resj[1,3]
  kpe =  hg19[,11] %in% infc[gnj,1]
  exj = hg19[kpe,]
  exj = exj[order(exj$V4, exj$V5), ]
  stj = paste(exj[,4], collapse=",")
  enj = paste(exj[,5], collapse=",")
  genj = sprintf("%s/%s.txt", sub.dir, resj[1,2])
  timj = sprintf("%s/%s_time.txt", sub.dir, resj[1,2])
  snps = sprintf("%s:%s-%s", chr, sti, eni)

  rasBASE = sprintf("rasqual -y %s -k %s -x %s -p %s -n %s", totb, kbin, xbin, np, nsam)
  rasPOSj = sprintf("-j %s -l %s -m %s -s %s -e %s", gnj, nl, nm, stj, enj)
  rasOPT = "-z --force --posterior-genotype" #--population-only --min-coverage-depth 0.0 --as-only 
  cmdi = sprintf("tabix %s %s | %s %s -f %s %s > %s",
    vcfgz, snps, rasBASE, rasPOSj, infc[gnj, 1], rasOPT, genj)
  
  timsta = proc.time()
    system(cmdi)
  timend = proc.time()  
  timend-timsta
  write.table(timend[3]-timsta[3], timj, row.names=F, col.names=F, quote=F)
}
print(cmdi)
message(gnj)

if(file.size(genj)>120){
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




# tabix ../data_genotype_all/vcfr/SNP_chr22.vcf.gz 22:50864181-51068485 | rasqual -y ../data_genotype_all/cnt/Tcnt_ex_chr22.bin -k ../data_genotype_all/cnt/K_ex_chr22.bin -x ../data_genotype_all/cnt/X_ex.bin -p 1 -n 376 -j 2 -l 3467 -m 3467 -s  50968485,50968428,50968428,50968183,50968138,50968138,50967767,50967767 -e 50968485,50968428,50968428,50968183,50968138,50968138,50967767,50967767 -f ENSG00000025708.13 -z --force --posterior-genotype 
