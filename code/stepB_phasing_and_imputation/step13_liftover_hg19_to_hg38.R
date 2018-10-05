setwd("/fh/fast/sun_w/licai/eQTL_KIRC/data_snp")
snps_dirs  = list.dirs()
snps_dirs[1:5]
# system('mkdir cgc_snpfile')
liftover = '/fh/fast/sun_w/bin/liftOverLinux/liftOver'
chain = "/fh/fast/sun_w/bin/liftOverLinux/hg19ToHg38.over.chain"

for(sam in snps_dirs[grep('TCGA',snps_dirs)]){
#creat bed file format
  system(sprintf("cd %s ; rm combined.txt ; for i in `seq 22` ; do cut -d' ' -f1-4 chr${i}.txt | awk '{print $1, $2-1, $2, $3, $4}' >> combined.txt ; done", sam)) 
#liftover to hg38
  system(sprintf('cd %s ; %s combined.txt %s combined_hg38.txt combined_unlifted.txt',sam, liftover, chain))
  system(sprintf("cd %s ;  cat combined_hg38.txt | cut -f1,3,4,5  >> %s", sam, paste0('../cgc_snpfile/',substr(sam, 3, 14),'.txt') ))
}

q(save ='no')