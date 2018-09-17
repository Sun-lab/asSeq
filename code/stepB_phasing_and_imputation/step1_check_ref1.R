shapeloc = "/fh/fast/sun_w/licai/_tumor_eQTL/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"
geno = "/fh/fast/sun_w/licai/eQTL_KIRC/data/chr"
refinfo = "/fh/fast/sun_w/licai/_tumor_eQTL/1000GP_Phase3"
cancer_type = 'KIRC'
proj = "/fh/fast/sun_w/licai/eQTL_KIRC/data/chr"
proj_log = "/fh/fast/sun_w/licai/eQTL_KIRC/R_batch1/stepB_phasing_and_imputation/log/step1_check_ref1"

system(sprintf("mkdir -p %s/precheck_log", geno))
system(sprintf("mkdir -p %s", proj_log))

all1000G="1000GP_Phase3" #common reference prefix

for(i in c(1:22)){
  pref = sprintf("time %s -check", shapeloc)
  pedi = sprintf("--input-ped %s/%s.chr%i", geno, cancer_type, i) # change cancer type
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

