#we use shapeit version 2 with hg19 aligned 1000G data
#this processing is done for each chromosome

# the location of shapeit and reference panel
shapeloc = "/fh/fast/sun_w/licai/_tumor_eQTL/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"
refinfo = "/fh/fast/sun_w/licai/_tumor_eQTL/1000GP_Phase3"
plinkloc = "plink"

# the location of chromoson level data after flipped
geno = "../../data/chr_flipped"
gwas_alignments = "../../data/chr_flipped/precheck_log"

# location to store output 
proj_log = "./log/step5_shapeit"
phased = "../../data/phased"

system(sprintf("mkdir -p %s", phased))
system(sprintf("chmod 777 %s", phased))
system(sprintf("mkdir -p %s", proj_log))

all1000G = "1000GP_Phase3" #common reference prefix
mem=32000
ncpu = 6
for(i in c(1:22)){
  #i = 16;
  pref = shapeloc
  pedi = sprintf("--input-ped %s/geno.chr%s_flip",geno,i)
  mapr = sprintf("--input-map %s/genetic_map_chr%s_combined_b37.txt",refinfo,i)
  inpr = sprintf("--input-ref %s/%s_chr%s.hap.gz",refinfo,all1000G,i)
  legr = sprintf("%s/%s_chr%s.legend.gz",refinfo,all1000G,i)
  samr = sprintf("%s/%s.sample",refinfo,all1000G)
  excl = sprintf("--exclude-snp %s/geno.chr%s_flip_check.snp.strand.exclude",gwas_alignments,i)
  out = sprintf("-O %s/phasing_chr%s",phased,i)
  oth = sprintf("--effective-size 20000 --seed 1234567 --thread %s",ncpu)
  log_output = sprintf("> %s/shapeit_chr%s.o 2> %s/shapeit_chr%s.e",
                       proj_log, i, proj_log, i)
  com = sprintf("time %s %s %s %s %s %s %s %s %s %s ",
                pref, pedi, mapr, inpr, legr, samr, excl, out,  oth, log_output)

  com2 = sprintf("sbatch --mem=%s --wrap='%s'",mem,com)
  message(com)
  system(com2)
}
system(sprintf("mv shapeit* %s", shapeit_log))

sessionInfo()
q("no")

