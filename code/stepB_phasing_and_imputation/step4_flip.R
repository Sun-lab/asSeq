#we use shapeit version 2 with hg19 aligned 1000G data
#this processing is done for each chromosome
#at this stage we check if the snps labeled as "strand issue" 
#can be fixed by flipping them with plink
setwd('/fh/fast/sun_w/licai/eQTL_KIRC/data')

# ----------------------------------------------------------------------
# check if the snps labeled as "strand issue" 
# can be fixed by flipping them with plink
# ----------------------------------------------------------------------
system('ml plink/1.90')
cancer_type = 'KIRC'
shapeloc = "/fh/fast/sun_w/licai/_tumor_eQTL/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"
plinkloc = "plink"
geno = "/fh/fast/sun_w/licai/eQTL_KIRC/data/chr"
gwas_alignments = "/fh/fast/sun_w/licai/eQTL_KIRC/data/chr/precheck_log"
refinfo = "/fh/fast/sun_w/licai/_tumor_eQTL/1000GP_Phase3"
proj_output = '/fh/fast/sun_w/licai/eQTL_KIRC/data/chr_flipped'
proj_log = '/fh/fast/sun_w/licai/eQTL_KIRC/data/chr_flipped/precheck_log'
checki<-matrix(NA,nrow=22,ncol=2)
rownames(checki)=sprintf("chr%s",1:22)

system(sprintf("mkdir -p %s", proj_output))
system(sprintf("mkdir -p %s", proj_log))


#as with each shapeit we do it on chromosome level
all1000G = "1000GP_Phase3" #common reference prefix
for(i in c(1:22)){
  #i = 22;
  
  #flip the snps of interest
  con = file(sprintf("%s/gwas.alignments_chr%d.snp.strand",gwas_alignments,i))
  lns = readLines(con)
  lns = lns[substr(lns,1,6)=="Strand"]
  strand = matrix(unlist(strsplit(lns,split="\t")),nrow=11)[4,]
  close(con)
  write.table(unique(strand),"tmp_strand.txt",row.names=F,col.names=F,quote=F)
  file.copy(from = sprintf("%s/%s.chr%s.map",geno,cancer_type,i), to = "tmp.map")
  file.copy(from = sprintf("%s/%s.chr%s.ped",geno,cancer_type,i), to = "tmp.ped")
  flipped = sprintf("%s/%s.chr%s_flip",proj_output,cancer_type,i)
  system(sprintf("%s --file tmp --flip tmp_strand.txt  --recode --out %s", plinkloc, flipped))
  cat("plink done ")
  
  #recheck using shapeit on flipped data
  pref = sprintf("%s -check", shapeloc)
  pedi = sprintf("--input-ped %s", flipped)
  mapr = sprintf("--input-map %s/genetic_map_chr%s_combined_b37.txt",refinfo,i)
  inpr = sprintf("--input-ref %s/%s_chr%s.hap.gz",refinfo,all1000G,i)
  legr = sprintf("%s/%s_chr%s.legend.gz",refinfo,all1000G,i)
  samr = sprintf("%s/%s.sample",refinfo,all1000G)
  out = sprintf("--output-log %s/%s.chr%s_flip_check", proj_log,cancer_type, i)
  com = sprintf("%s %s %s %s %s %s %s",
                pref, pedi, mapr, inpr, legr, samr, out)
  message(com) 
  system(com)
  
  message("processed ",i,"'th file")
  
  #remove temporary files
  system("rm -f tmp*")
  system("rm -f plink*")
}

sessionInfo()
q("no")
