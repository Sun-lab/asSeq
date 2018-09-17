#we use shapeit version 2 with hg19 aligned 1000G data
#this processing is done for each chromosome
cancer_type = 'KIRC'
phased = "/fh/fast/sun_w/licai/eQTL_KIRC/data/phased"
shapeloc = "/fh/fast/sun_w/licai/_tumor_eQTL/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"
plinkloc = "plink"
geno = "/fh/fast/sun_w/licai/eQTL_KIRC/data/chr_flipped"
gwas_alignments = "/fh/fast/sun_w/licai/eQTL_KIRC/data/chr_flipped/precheck_log"
refinfo = "/fh/fast/sun_w/licai/_tumor_eQTL/1000GP_Phase3"
pipedir = '/fh/fast/sun_w/licai/eQTL_KIRC/R_batch1/stepB_phasing_and_imputation/pipline'
proj_log = "../log/step5_shapeit"

system(sprintf("mkdir -p %s", phased))
system(sprintf("chmod 777 %s", phased))
system(sprintf("mkdir -p %s", proj_log))

all1000G = "1000GP_Phase3" #common reference prefix
mem=32000
ncpu = 6
for(i in c(1:22)){
  #i = 16;
  pref = shapeloc
  pedi = sprintf("--input-ped %s/%s.chr%s_flip",geno,cancer_type,i)
  mapr = sprintf("--input-map %s/genetic_map_chr%s_combined_b37.txt",refinfo,i)
  #  inpr = sprintf("--input-ref %s/%s_chr%s_impute.hap.gz",refinfo,all1000G,i)
  #  legr = sprintf("%s/%s_chr%s_impute.legend.gz",refinfo,all1000G,i)
  inpr = sprintf("--input-ref %s/%s_chr%s.hap.gz",refinfo,all1000G,i)
  legr = sprintf("%s/%s_chr%s.legend.gz",refinfo,all1000G,i)
  samr = sprintf("%s/%s.sample",refinfo,all1000G)
  excl = sprintf("--exclude-snp %s/%s.chr%s_flip_check.snp.strand.exclude",gwas_alignments,cancer_type,i)
  out = sprintf("-O %s/phasing_chr%s",phased,i)
  # output_log = sprintf("--output-log %s/log/shapeit_chr%s",phased,i)
  oth = sprintf("--effective-size 20000 --seed 1234567 --thread %s",ncpu)
  log_output = sprintf("> %s/shapeit_chr%s.o 2> %s/shapeit_chr%s.e",
                       proj_log, i, proj_log, i)
  com = sprintf("time %s %s %s %s %s %s %s %s %s %s ",
                pref, pedi, mapr, inpr, legr, samr, excl, out,  oth, log_output)
  
  # qout = sprintf("--output=%s/step5_out_%s.out",pipedir,i)    
  com2 = sprintf("sbatch --mem=%s --wrap='%s'",mem,com)
  message(com)
  system(com2)
}
system(sprintf("mv shapeit* %s", shapeit_log))

sessionInfo()
q("no")

