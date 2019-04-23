#module add samtools; module remove java; module add java; module add r;module add tabix
options(digits=22)
data.dir = "../data"
bam.dir = sprintf("%s/bam", data.dir)
geno.dir = "../datagen" 
vcf.dir = sprintf("%s/vcf", geno.dir)
info.dir = "../inf"
ref = sprintf("%s/hg38.fa", info.dir)
dic = sprintf("%s/hg38.dict", info.dir)
dir.out = sprintf("%s/gatkasc", data.dir)
if(!file.exists(dir.out))dir.create(dir.out)


#need to put picard and GenomeAnalysisTK jar files next to work folder
pic = "../gatk/picard-2.8.2.jar"
gat = "../gatk/GenomeAnalysisTK.jar"
#or set the path to the location
#may try exporting path
#export PATH=/nas02/home/z/h/zhabotyn/research/:$PATH
#pic = "picard-2.8.2.jar"
#gat = "GenomeAnalysisTK.jar"

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
snpt = sprintf("%s/SNP_ex.vcf.gz", vcf.dir)
snp2 = sprintf("%s/SNP_fix.vcf", vcf.dir)

#sometimes people write vcf files with spaces not tabs, but need to have as tabs
hdr = read.table(gzfile(snpt), as.is=T, comment.char="", nrows=4, sep="@")
hdr[4,1] = gsub(" ", "\t", hdr[4,1])
snps = read.table(gzfile(snpt), as.is=T)
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
snps[,-(1:9)] = apply(snps[,-(1:9)], 1:2, get_block, split=":")
#write in the resulting file
write.table(hdr, snp2, quote=F, row.names=F, col.names=F, sep="\t")
write.table(snps, snp2, quote=F, row.names=F, col.names=F, sep="\t", append=T)
system(sprintf("bgzip -f %s", snp2))
snpt = sprintf("%s.gz", snp2)
system(sprintf("tabix %s", snpt))

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
    fil.out1 = sprintf("%s/%sf.bam", dir.out, nms[k])
    boup = gsub("\\.bam","\\.out", fil.out)
    rst = "SO=coordinate VALIDATION_STRINGENCY=SILENT RGLB=lib1 RGPL=illumina RGPU=unit1" 
    com = com1 = sprintf("java -jar %s AddOrReplaceReadGroups I=%s O=%s %s RGSM=%s", 
                                             pic, fil.inp, fil.out, rst, nms[k])
    if(!file.exists(dic)){
      com0 = sprintf("java -jar %s CreateSequenceDictionary R=%s O=%s", pic, ref, dic)
      system(com0)
    }
    if(!file.exists(sprintf("%s.fai", ref))){
      com0 = sprintf("samtools faidx %s", ref)
      system(com0)
    }
    #preprocessing - fix mate information
    scom0 = sprintf("samtools sort -n -o %s %s", fil.out, fil.inp)
    scom1 = sprintf("samtools fixmate -r %s %s", fil.out, fil.out1)
    scom2 = sprintf("samtools sort -o %s %s", fil.out, fil.out1)    
    scom3 = sprintf("rm %s", fil.out1)
    com0 = sprintf("%s;%s;%s; %s", scom0, scom1, scom2, scom3)
    system(com0)
    
    #first need to ensure that fasta chromosomes are in the same order as bam
    rst = "VALIDATION_STRINGENCY=SILENT "
    #rst = "VALIDATION_STRINGENCY=STRICT "
    fil.out2 = sprintf("%s/ob_%s", dir.out, fls[k])
    com1 = sprintf("java -jar %s ReorderSam I=%s O=%s R=%s %s", 
                        pic, fil.out, fil.out2, ref, rst)
    #system(com1)

    #need to add "group name" - sample id and index
    rst = "SO=coordinate VALIDATION_STRINGENCY=SILENT RGLB=lib1 RGPL=illumina RGPU=unit1" 
    com2 = sprintf("java -jar %s AddOrReplaceReadGroups I=%s O=%s %s RGSM=%s", 
                                             pic, fil.out2, fil.out, rst, nms[k])
    com3 = sprintf("rm %s", fil.out2) 
    com4 = sprintf("samtools index %s", fil.out)
      

    #now we got to counting
    rst = "-U ALLOW_N_CIGAR_READS -minDepth 1 --minMappingQuality 10 --minBaseQuality 20"
    fil.out3 = gsub("\\.bam","\\.txt", fil.out)

    com5 = sprintf("java -Xmx7g -jar %s -T ASEReadCounter -R %s -I %s -o %s -sites %s %s", 
                   gat, ref, fil.out, fil.out3, snpt, rst)
    com = sprintf("%s && %s && %s && %s && %s", com1, com2, com3, com4, com5)
  }
  #can add this command, but may start keeping the file at the beginning
  com6 = sprintf("rm -f %s*", fil.out)
  #com = sprintf("%s;%s", com, com6)
  system(com)
  #note, test example is short so it is possible just to run them as one command
  #in real dataset it is better to submit each sample separately in a bcom style:
  #bcom = sprintf("sbatch --output %s --mem=16G --wrap='%s'",boup, com)
  #system(bcom)
  message(k)
}

q("no")
