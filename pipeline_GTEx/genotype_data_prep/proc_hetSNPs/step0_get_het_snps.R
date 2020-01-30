args=(commandArgs(TRUE))
# args = c("21")
chri = args[1]

library(Matrix)
pars4gr = "/nas/longleaf/home/zhabotyn/progs/parser/parser4grNA"

base.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8"
vcf.dir = sprintf("%s/WGS_VCF", base.dir)
snp.dir = sprintf("%s/het_snps", vcf.dir)
if(!file.exists(snp.dir))dir.create(snp.dir)
snp.chr = sprintf("%s/chr%s", snp.dir, chri)
if(!file.exists(snp.chr))dir.create(snp.chr)

vcf.fli = sprintf("%s/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_chr%s.vcf.gz", 
vcf.dir, chri)


tmpfl = gzfile(vcf.fli)
lns = readLines(tmpfl, n=4)
close(tmpfl)
samvcf = gsub("#", "", unlist(strsplit(lns[4], split="\t")))
samvcf = samvcf[-(1:9)]

#get SNP information
comi = sprintf("zcat %s  | sed 's/ /\t/g' | sed /^#/d  | cut  -f '1-5,10'", vcf.fli)
p = pipe(comi)
info = read.table(p, fill=TRUE, row.names=NULL)
info[1:2,]
info[,1] = gsub("chr", "", info[,1])
info[,4] = as.character(info[,4])
info[,5] = as.character(info[,5])

info[1:2,]
info[1:5,]

#get SNP classificaitons
comi = sprintf("zcat %s  | sed 's/ /\t/g'  | sed /^#/d  | cut  -f '10-' | %s | cut -f '1-3'", vcf.fli, pars4gr)
p = pipe(comi)
geno = read.table(p, colClasses=c("integer"), fill=TRUE, row.names=NULL)
apply(geno, 2, max)
geno[1:4,]
table(geno[,3])
geno  = t(sparseMatrix(i=geno[,2], j=geno[,1], x=geno[,3]))
info = info[1:nrow(geno),]

dim(info)
dim(geno)

kp = nchar(info[,4])==1 & nchar(info[,5])==1; table(kp)
info = info[kp,]
geno = geno[kp,]

dim(info)
dim(geno)

for(indi in 1:ncol(geno)){
  kpi = geno[,indi] %in% 1:2
  cbind(info, geno[,indi])[kpi,][1:25,]
  #1 is 0|1 and 2 is 1|0
  #vcf file records 0 as reference (v4) and 1 as alternative)
  #i.e. for 1 we put ref|alt and for 2 we put alt|ref
  #otherwise - just take ref and alt values and for 2 switch them
  
  
  #final output 
  #four columns, chromosome, position, allele 1 and allele 2, without header
  #for intermediate I will also keep fifth column - snp id
  snps = info[kpi,c(1,2, 4:5, 3)]
  flip = info[kpi,c(1,2, 5:4, 3)]
  stat = geno[kpi, indi]
  cbind(snps, stat)[1:10,]
  cbind(flip, stat)[1:10,]
  snps[stat==2,] = flip[stat==2,]
  cbind(snps, stat)[1:10,]
  
  chriout = sprintf("%s/%s_chr%s.txt", snp.chr, samvcf[indi], chri)
  chriout
  write.table(snps, chriout, row.names=F, col.names=F, quote=F, sep="\t")
  if(indi %% 100 == 0) message(indi)
}

q("no")
