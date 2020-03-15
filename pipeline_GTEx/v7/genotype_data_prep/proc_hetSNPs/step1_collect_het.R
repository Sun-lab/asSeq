chri = 21
base.dir = "/pine/scr/z/h/zhabotyn/R01/GTEx/v8"
vcf.dir = sprintf("%s/WGS_VCF", base.dir)
snp.dir = sprintf("%s/het_snps", vcf.dir)
byindwid = sprintf("%s/byindwid", snp.dir)
byindnid = sprintf("%s/byindnid", snp.dir)
if(!file.exists(byindwid))dir.create(byindwid)
if(!file.exists(byindnid))dir.create(byindnid)

vcf.fli = sprintf("%s/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_chr%s.vcf.gz", 
vcf.dir, chri)


tmpfl = gzfile(vcf.fli)
lns = readLines(tmpfl, n=4)
close(tmpfl)
samvcf = gsub("#", "", unlist(strsplit(lns[4], split="\t")))
samvcf = samvcf[-(1:9)]

numdbl = numsnps2 = numsnps = rep(0, length(samvcf))

length(samvcf)
indi = 1

for(indi in 1:length(samvcf)){
  chriin = sprintf("%s/chr%s/%s_chr%s.txt", snp.dir, 1:22, samvcf[indi], 1:22)
  outwid = sprintf("%s/%s.txt", byindwid, samvcf[indi])
  outnid = sprintf("%s/%s.txt", byindnid, samvcf[indi])
  nsnp = sprintf("%s/nsnp_%s.txt", byindwid, samvcf[indi])
  chriin
  com = paste(c("cat", chriin, ">", outwid), collapse=" ")
  system(com)
  com = sprintf("sed -i 's/ /\t/g' %s", outwid)
  system(com)
  
  com = sprintf("wc -l %s > %s", outwid, nsnp)
  system(com)
  
  com = sprintf("cat %s | cut -f '1-4' > %s", outwid, outnid)
  system(com)

  if(indi %% 1 == 0)message(indi)
}


indin = list.files(path=byindwid, pattern="nsnp");indin[1:5]
setwd(byindwid)
nsnps = sprintf("%s/num_snps.txt", snp.dir) 
com = paste(c("cat", indin, ">", nsnps), collapse=" ")
system(com)

