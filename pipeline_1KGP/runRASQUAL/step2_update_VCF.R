#module add tabix;module add gcc;module add gsl

data.dir = "../omnidata"
bam.dir = sprintf("%s/ufilt", data.dir)
info.dir = "../info"
ref = sprintf("%s/hg38.fa", info.dir)
vcf.dir = "../data_genotype/465ind/rvcf"
vcf.dirc = sprintf("%s/chrs", vcf.dir)
dir.out = sprintf("%s/gatkasc", data.dir)

tre.dir = sprintf("%s/parplus_conv",data.dir)
spe.dir = sprintf("%s/out2e5", pipe.dir)
ras.dir = sprintf("%s/rasq",data.dir)

vcf.dir = sprintf("%s/data_genotype/465ind/vcf/chrs", root.dir)
rvcf.dir = sprintf("%s/data_genotype/465ind/rvcf", root.dir)
if(!file.exists(rvcf.dir))dir.create(rvcf.dir)
#vcf.dir = sprintf("%s/data_genotype/465ind/vcfsub", root.dir)
pipe.dir = sprintf("%s/omnipipe", root.dir)
rasq.dir = sprintf("/nas02/home/z/h/zhabotyn/rasqual-master/bin")
#tre.dir = sprintf("%s/parplus_conv",data.dir)
spe.dir = sprintf("%s/ras", pipe.dir)
asc.dir = sprintf("%s/asc", spe.dir)
if(!file.exists(spe.dir))dir.create(spe.dir)
if(!file.exists(asc.dir))dir.create(asc.dir)
setwd(asc.dir)


count_within = function(i, spos, gns, off=1e5, off2=2e5){
  fl = (spos>=gns[i,4])&(spos<gns[i,5])
  fl2 = (spos>=(gns[i,4]-off))&(spos<(gns[i,5]+off))
  fl3 = (spos>=(gns[i,4]-off2))&(spos<(gns[i,5]+off2))
  c(sum(fl), sum(fl2), sum(fl3))
}
snp = sprintf("%s/SNP_VCF_tab.vcf", vcf.dir)
system(sprintf("awk '{print $1,$2}' %s > tmpsnps.txt", snp))

#hdr = read.table(snp, as.is=T, comment.char="", nrows=4, sep="@")
#snps = read.table(snp, as.is=T)

i = 1
for(i in 1:22){
  chri = sprintf("chr%s", i)
  system(sprintf("grep '%s ' tmpsnps.txt > chrisnps.txt", chri))  
  fli = sprintf("%s/gene_info_chr%s.dat",tre.dir,i)
  fli = read.table(fli, header=F, as.is=T)
  #resi = matrix(0, nrow=nrow(fli), ncol=3)
  posi = read.table("chrisnps.txt", as.is=T)
  #posi[1:4,]
  #table(posi[,1])
  posi = posi[posi[,1]==chri,2] 
  count_within(1, spos = posi, gns=fli)
  resi = sapply(1:nrow(fli), count_within, spos=posi, gns=fli)
#  colnames(resi) = fli[,1]
  
  if(i==1){
    res = data.frame(ids=as.character(fli[,1]), within=as.numeric(resi[1,]), w1e5=as.numeric(resi[2,]), w2e5=as.numeric(resi[3,]), chr=chri)
  }else{
    resi = data.frame(ids=as.character(fli[,1]), within=as.numeric(resi[1,]), w1e5=as.numeric(resi[2,]), w2e5=as.numeric(resi[3,]), chr=chri)
    res = rbind(res, resi)
  }
  message(i)
}
write.csv(res, "number_of_snps_per_gene.csv", quote=F, row.names=F) 


q("no")
