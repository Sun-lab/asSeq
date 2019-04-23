#module add tabix;module add gcc;module add gsl

set.seed(1565691)

data.dir = "../data"
cnt.dir = sprintf("%s/cnt", data.dir)
bam.dir = sprintf("%s/bam", data.dir)
geno.dir = "../datagen" 
vcf.dir = sprintf("%s/vcf", geno.dir)
#store permuted data to a separate folder from
vcfr.dir = sprintf("%s/vcfp", geno.dir)
if(!file.exists(vcfr.dir))dir.create(vcfr.dir)
info.dir = "../inf"
dir.rasc = sprintf("%s/gatkasc", data.dir)

snp = sprintf("%s/SNP_fix.vcf.gz", vcf.dir)
system(sprintf("gunzip -c %s | awk '{print $1,$2}' > tmpsnps.txt", snp))

genes = read.table(sprintf("%s/GTF_ex.gtf", info.dir), as.is=T)
geneid = unique(genes[,11])
to.rm = grep(":", geneid)
if(length(to.rm)>0)geneid=geneid[-to.rm]

res = data.frame(matrix(NA, nrow=length(geneid), ncol=3))
win = 2e5
for(i in 1:length(geneid)){
  #m = match(geneid[i], genes[,11])
  gnm = which(geneid[i]==genes[,11])
  chri = genes[gnm[1],2]
  system(sprintf("grep '%s ' tmpsnps.txt > chrisnps.txt", chri))  
  chrsnps = read.table("chrisnps.txt")
  posst = min(genes[gnm,5])
  posen = max(genes[gnm,6])
  kp = which(((posst-win)<=chrsnps[,2]) & ((posen+win)>chrsnps[,2]))
  res[i, 1] = chri
  res[i, 2] = geneid[i]
  res[i, 3] = length(kp)
  message(i)
}
write.csv(res, sprintf("%s/number_of_snps_per_gene_win%s.csv", info.dir, win), quote=F, row.names=F) 

#add vcf counts
vcf.hd = read.table(snp, as.is=T, nrow=4, comment.char=",", sep=",")
vcf.in = read.table(snp, as.is=T)

#AF in this file came from a bigger file, 
#replace it with frequency appropriate for these samples
count.ref = function(vec){mean(unlist(strsplit(vec, split="\\|"))==1)}
af = apply(vcf.in[,-(1:9)], 1, count.ref)
vcf.in[,8] = af

#if too many reads of wrong bases - may skip this count
#also add a permutation step for VCF
frac = 0.1
samples = unlist(strsplit(vcf.hd[4,1], split="\t"))[-(1:9)]
for(i in 1:length(samples)){
  alt = ref = numeric(nrow(vcf.in))
  sampli = samples[i]
  rasc = read.table(sprintf("%s/%s.txt",dir.rasc, sampli), header=T)
  m = match(rasc$variantID, vcf.in[,3])
  if(any(is.na(m))){warning(sum(is.na(m)), " unexpected snps")}
  ref[na.omit(m)] = rasc$refCount
  alt[na.omit(m)] = rasc$altCount
  asc = ref[na.omit(m)]+alt[na.omit(m)]
  rma = rasc$otherBases/asc>frac
  ref[na.omit(m)][rma] = alt[na.omit(m)][rma] = 0
  vcf.in[,i+9] = sprintf("%s:%s,%s", sample(vcf.in[,i+9]), ref, alt)
  message("individual ", i, " has been processed")
}

#save to rasqual count folder
out.vcf = sprintf("%s/SNP_fix.vcf", vcfr.dir)
write.table(vcf.hd, out.vcf, quote=F, row.names=F, col.names=F, sep="\t")
write.table(vcf.in, out.vcf, quote=F, row.names=F, col.names=F, sep="\t", append=T)
system(sprintf("bgzip -f %s", out.vcf))
out.vcf = sprintf("%s.gz", out.vcf)
system(sprintf("tabix %s", out.vcf))


#create other required files if they don't exist yet
totc = sprintf("%s/Tcnt_ex.dat", cnt.dir)
xcnt = sprintf("%s/samples.dat", cnt.dir)
totc = read.table(totc, as.is=T)
xcnt = read.table(xcnt, as.is=T, header=T)

totb = sprintf("%s/Tcnt_ex.bin", cnt.dir)
kbin = sprintf("%s/K_ex.bin", cnt.dir)
xbin = sprintf("%s/X_ex.bin", cnt.dir)

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

kin = readBin(kbin, n=200, double())
kin
tin = readBin(totb, n=200, double())
tin
xin = readBin(xbin, n=200, double())
xin


q("no")
