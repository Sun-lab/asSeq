#ml samtools 
#ml picard
#ml GATK/3.7-Java-1.8.0_92 
setwd("/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL")
options(digits=22)
mem = 7
data.dir = "/fh/scratch/delete90/sun_w/_GTEx/ncbi"
bam.dir = sprintf("%s/wb_raw_bamfiles", data.dir)
geno.dir = "../data_genotype_all" 
vcf.dir = sprintf("%s/vcf", geno.dir)
info.dir = "."
ref = sprintf("%s/Homo_sapiens_assembly19.fasta", info.dir) 
dic = sprintf("%s/Homo_sapiens_assembly19.dict", info.dir)
dir.out = sprintf("%s/gatkasc", data.dir)
if(!file.exists(dir.out))dir.create(dir.out)


#need to put picard and GenomeAnalysisTK jar files next to work folder
# pic = "../gatk/picard-2.8.2.jar"
pic = "$EBROOTPICARD/picard.jar"
gat = "/app/easybuild/software/GATK/3.7-Java-1.8.0_92/GenomeAnalysisTK.jar"
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
snpt = sprintf("%s/GTEx_phased.vcf.gz", vcf.dir)
snp2 = sprintf("%s/plink.vcf", vcf.dir)

#
#block a
#
#sometimes people write vcf files with spaces not tabs, but need to have as tabs
# hdr = read.table(gzfile(snpt), as.is=T, comment.char="", nrows=217, sep="@")
# hdr[217,1] = gsub(" ", "\t", hdr[4,1])
# snps = read.table(gzfile(snpt), as.is=T)
#double-check if any of SNP positions is duplicated and remove it
# to.rm = names(table(snps[,2])[table(snps[,2])>1])
# if(length(to.rm)>0){
#   m = which(is.na(match(snps[,2], to.rm)))
#   dim(snps)
#   snps = snps[m,]
#   dim(snps)
# }


#getting sample ids for which we want to do allele-specific counts
vcfnm = read.table(snpt, comment.char="", nrows=1, skip=27, as.is=T)
samp = as.character(vcfnm[-(1:9)])
samp = gsub("_", "", substr(samp, 1, 10))
length(samp)
samp[1:5]

#match SRR with sample ID
proj_dir   = "/fh/fast/sun_w/licai/_GTEx/"
sraWB_file = "SraRunInfo_whole_blood.csv"

sraWB = read.csv(sprintf("%s%s", proj_dir, sraWB_file), as.is = T)
sraWB[1:2,]
sraWB$sam = gsub("-$","", sapply(sraWB$SampleName, substr, 1, 10))
dim(sraWB)
table(sraWB$Submission)
sraWB = sraWB[which(sraWB$Submission != "SRA807555"), ]
dim(sraWB)
table(sraWB$Submission)
ssr2kp = read.table("/fh/fast/sun_w/licai/_GTEx/R_batch4_whole_blood/stepC/output/step5_TReC_per_gene_filterIt/RunID_used_in_GTEx_wb.txt")
ssr2kp = ssr2kp$V1

samMatch = sraWB$sam[match(nms, sraWB$Run)]
table(samMatch %in% samp)

for(k in 1:length(fls)){
#k = 1
  if(nms[k] %in% ssr2kp){
  #k = 1
    fil.inp = sprintf("%s/%s", bam.dir, fls[k])
    #first we had to do several fixes to bam files to make GATK run smoothly    
    fil.out  = sprintf("%s/%s.bam", dir.out, samMatch[k])
    fil.out1 = sprintf("%s/%sf.bam", dir.out, samMatch[k])
    fil.outS = sprintf("%s/../wb_sort_bamfiles/%s_filtered_sorted_byQname.bam",
                       dir.out, nms[k])
    boup = gsub("\\.bam","\\.out", fil.out)

    if(!file.exists(dic)){
      file.remove(dic)
      com0 = sprintf("java -jar %s CreateSequenceDictionary R=%s O=%s", pic, ref, dic)
      system(com0)
    }
    if(!file.exists(sprintf("%s.fai", ref))){
      com0 = sprintf("samtools faidx %s", ref)
      system(com0)
    }
    #
    #b preprocessing - fix mate information
    #
    # scom0 = sprintf("samtools sort -n -o %s %s", fil.out, fil.inp) # have done it
    scom1 = sprintf("samtools fixmate -r %s %s", fil.outS, fil.out1)
    scom2 = sprintf("samtools sort %s %s", fil.out1, gsub(".bam$","",fil.out))
    scom3 = sprintf("rm %s; samtools index %s", fil.out1, fil.out)
    com0 = sprintf("%s && %s && %s", scom1, scom2, scom3)
    # system(com0)
    
    #
    #c. reorder bam file chromosom info to match reference FASTA file
    #first need to ensure that fasta chromosomes are in the same order as bam
    #
    rst = "VALIDATION_STRINGENCY=SILENT "
    #rst = "VALIDATION_STRINGENCY=STRICT "
    fil.out2 = sprintf("%s/ob_%s", dir.out, fls[k])
    if(file.exists(fil.out2))file.remove(fil.out2)
    com1 = sprintf("java -jar %s ReorderSam I=%s O=%s R=%s %s", 
                        pic, fil.out, fil.out2, ref, rst)
    # system(com1)

    #
    #d. need to add "group name" - sample id and index
    #
    rst = "SO=coordinate VALIDATION_STRINGENCY=SILENT RGLB=lib1 RGPL=illumina RGPU=unit1" 
    com2 = sprintf("java -jar %s AddOrReplaceReadGroups I=%s O=%s %s RGSM=%s", 
                                             pic, fil.out2, fil.out, rst, samMatch[k])
    # system(com2)
    com3 = sprintf("rm %s", fil.out2) 
    # system(com3)
    com4 = sprintf("samtools index %s", fil.out)
    # system(com4)
      

    #
    #e. now we got to counting
    #
    rst = "-U ALLOW_N_CIGAR_READS -minDepth 1 --minMappingQuality 10 --minBaseQuality 20"
    fil.out3 = gsub("\\.bam","\\.txt", fil.out)
    com5 = sprintf("java -Xmx%sg -jar %s -T ASEReadCounter -R %s -I %s -o %s -sites %s %s", 
                   mem, gat, ref, fil.out, fil.out3, snpt, rst)
    #system(com5)
    com = sprintf("%s && %s && %s && %s && %s && %s",com0, com1, com2, com3, com4, com5)
    # system(com)
    #can add this command, but may start keeping the file at the beginning
    #if(file.exists(fil.out))system(com6)  
    #note, test example is short so it is possible just to run them as one command
    #in real dataset it is better to submit each sample separately in a bcom style:
    bcom = sprintf("sbatch --exclusive --output %s --wrap='%s'",boup, com)
    message(bcom)
    system(bcom)
    # com6 = sprintf("rm -f %s*", fil.out)
    message(k)
    
  }

}

q("no")
