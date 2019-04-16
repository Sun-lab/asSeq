#module add samtools; module remove java; module add java
options(digits=22)
data.dir = "../omnidata"
bam.dir = sprintf("%s/ufilt", data.dir)
info.dir = "../info"
ref = sprintf("%s/hg38.fa", info.dir)
vcf.dir = "../data_genotype/465ind/rvcf"
vcf.dirc = sprintf("%s/chrs", vcf.dir)
dir.out = sprintf("%s/gatkasc", data.dir)
if(!file.exists(dir.out))dir.create(dir.out)

spe.dir = "ras"
asc.dir = sprintf("%s/asc", spe.dir)
if(!file.exists(spe.dir))dir.create(spe.dir)
if(!file.exists(asc.dir))dir.create(asc.dir)

#need to put picard and GenomeAnalysisTK jar files next to work folder
pic = "../gatk/picard-2.8.2.jar"
gat = "../gatk/GenomeAnalysisTK.jar"

#a function to keep only base part of a sample - need to change split depending on a file name
get_block = function(str, split="\\.bam",block=1){
  unlist(strsplit(str, split=split))[block]
}

#keep only bam files
fls = list.files(bam.dir,pattern="\\.bam")
#if had previous preprocessing into allele-specific subfiles ignore them
if(length(grep("_hap", fls))>0)fls = fls[-grep("_hap", fls)]
nms = sapply(fls, get_block)
nms = sapply(nms, get_block, split="_", block=1)

#if we had previous TRECASE_MLE processing vcf files should already exist
snpt = sprintf("%s/SNP_VCF_tab.vcf", vcf.dirc)
for(i in 1:22){
  snp = sprintf("%s/SNP_VCF_chr%s.vcf", vcf.dirc, i)#input vcf
  #sometimes people write vcf files with spaces not tabs, but need to have as tabs
  hdr = read.table(snp, as.is=T, comment.char="", nrows=4, sep="@")
  hdr[4,1] = gsub(" ", "\t", hdr[4,1])
  snps = read.table(snp, as.is=T)
  #double-check if any of SNP positions is duplicated and remove it
  to.rm = names(table(snps[,2])[table(snps[,2])>1])
  if(length(to.rm)>0){
    m = which(is.na(match(snps[,2], to.rm)))
    dim(snps)
    snps = snps[m,]
    dim(snps)
  }
  #make sure it is ordered
  o = order(snps[,1], as.numeric(as.character(snps[,2])))
  snps = snps[o,]
  #write in the resulting file
  if(i==1)write.table(hdr, snp2, quote=F, row.names=F, col.names=F, sep="\t")
  write.table(snps, snpt, quote=F, row.names=F, col.names=F, sep="\t", append=T)

  message(i)
}

#getting sample ids for which we want to do allele-specific counts
vcfnm = read.table(snpt, comment.char="", nrows=1, skip=3, as.is=T)
samp = as.character(vcfnm[-(1:9)])
length(samp)

for(k in 1:length(fls)){
#k = 1
  if(nms[k] %in% samp){
    fil.inp = sprintf("%s/%s", bam.dir, fls[k])
    #first we had to do several fixes to bam files to make GATK run smoothly    
    fil.out = sprintf("%s/%s.bam", dir.out, nms[k])
    boup = gsub("\\.bam","\\.out", fil.out)
    rst = "SO=coordinate VALIDATION_STRINGENCY=SILENT RGLB=lib1 RGPL=illumina RGPU=unit1" 
    com = com1 = sprintf("java -jar %s AddOrReplaceReadGroups I=%s O=%s %s RGSM=%s", 
                                             pic, fil.inp, fil.out, rst, nms[k])
  
    rst = "VALIDATION_STRINGENCY=SILENT "
    fil.out2 = sprintf("%s/ob_%s", dir.out, fls[k])
    com2 = sprintf("java -jar %s ReorderSam I=%s O=%s R=%s %s", 
                        pic, fil.out, fil.out2, ref, rst)
    com3 = sprintf("mv %s %s", fil.out2, fil.out)
    com4 = sprintf("samtools index %s", fil.out)
    com = sprintf("%s && %s && %s && %s", com, com2, com3, com4)
      
    rst = "-U ALLOW_N_CIGAR_READS -minDepth 1 --minMappingQuality 10 --minBaseQuality 20"

    fil.out2 = gsub("\\.bam","\\.txt", fil.out)

    com5 = sprintf("java -Xmx7g -jar %s -T ASEReadCounter -R %s -I %s -o %s -sites %s %s", 
                   gat, ref, fil.out, fil.out2, snpt, rst)
    com = sprintf("%s && %s", com, com5)
  }
  #can add this command, but may start keeping the file at the beginning
  com6 = sprintf("rm -f %s*", oup)

  bcom = sprintf("sbatch --output %s --mem=16G --wrap='%s'",boup, com)
  system(bcom)
}

q("no")
