# the location of shapeit and reference panel
shapeloc = "/fh/fast/sun_w/licai/_tumor_eQTL/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"
refinfo = "/fh/fast/sun_w/licai/_tumor_eQTL/1000GP_Phase3"

# the location of chromoson level data
proj = "../../data/chr"

# location to store output log file 
proj_log = "./log/step1_check_ref1"

system(sprintf("mkdir -p %s/precheck_log", proj))
system(sprintf("mkdir -p %s", proj_log))

all1000G="1000GP_Phase3" #common reference prefix

for(i in seq(22)){
  pref = sprintf("time %s -check", shapeloc)
  pedi = sprintf("--input-ped %s/geno.chr%i", proj, i) 
  mapr = sprintf("--input-map %s/genetic_map_chr%s_combined_b37.txt",refinfo,i)
  inpr = sprintf("--input-ref %s/%s_chr%s.hap.gz",refinfo,all1000G,i)
  legr = sprintf("%s/%s_chr%s.legend.gz",refinfo,all1000G,i)
  samr = sprintf("%s/%s.sample",refinfo,all1000G)
  out = sprintf("--output-log %s/precheck_log/gwas.alignments_chr%s",proj,i)
  com = sprintf("%s %s %s %s %s %s %s",
                pref, pedi, mapr, inpr, legr, samr, out)
  message(com)
  system(com)
}

sessionInfo()
q("no")

