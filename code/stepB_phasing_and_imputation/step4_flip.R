#we use shapeit version 2 with hg19 aligned 1000G data
#this processing is done for each chromosome

# ----------------------------------------------------------------------
# check if the snps labeled as "strand issue" 
# can be fixed by flipping them with plink
# ----------------------------------------------------------------------

# the location of shapeit and reference panel
shapeloc = "/fh/fast/sun_w/licai/_tumor_eQTL/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"
refinfo = "/fh/fast/sun_w/licai/_tumor_eQTL/1000GP_Phase3"
plinkloc = "plink"

# the location of chromoson level data
geno = "../../data/chr"
gwas_alignments = "../../data/chr/precheck_log"

# location to store output log file 
proj_log = "./log/step4_flip"
proj_output = '../../data/chr_flipped'
proj_log = '../../data/chr_flipped/precheck_log'

system(sprintf("mkdir -p %s", proj_output))
system(sprintf("mkdir -p %s", proj_log))

#as with each shapeit we do it on chromosome level
all1000G = "1000GP_Phase3" #common reference prefix

# change seq(22) tp 22 to run the example data
for(i in seq(22)){
  #i = 22;
  #flip the snps of interest
  con = file(sprintf("%s/gwas.alignments_chr%d.snp.strand",
                     gwas_alignments,i))
  lns = readLines(con)
  lns = lns[substr(lns,1,6)=="Strand"] # find the SNP with strand issue
  strand = matrix(unlist(strsplit(lns,split="\t")),nrow=11)[4,]
  close(con)
  write.table(unique(strand),"tmp_strand.txt",row.names=F,
              col.names=F,quote=F)
  file.copy(from = sprintf("%s/geno.chr%s.map",geno,i), 
            to = "tmp.map")
  file.copy(from = sprintf("%s/geno.chr%s.ped",geno,i), 
            to = "tmp.ped")
  flipped = sprintf("%s/geno.chr%s_flip",proj_output,i)
  system(sprintf("%s --file tmp --flip tmp_strand.txt  --recode --out %s",
                 plinkloc, flipped)) # plink to flip the strand
  cat("plink done ")
  
  #recheck using shapeit on flipped data
  pref = sprintf("%s -check", shapeloc)
  pedi = sprintf("--input-ped %s", flipped)
  mapr = sprintf("--input-map %s/genetic_map_chr%s_combined_b37.txt",
                 refinfo,i)
  inpr = sprintf("--input-ref %s/%s_chr%s.hap.gz",refinfo,all1000G,i)
  legr = sprintf("%s/%s_chr%s.legend.gz",refinfo,all1000G,i)
  samr = sprintf("%s/%s.sample",refinfo,all1000G)
  out = sprintf("--output-log %s/geno.chr%s_flip_check", proj_log, i)
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
