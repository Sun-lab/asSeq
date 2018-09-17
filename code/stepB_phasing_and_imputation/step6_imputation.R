#imputation
#ensure that when we write positions to a file R doesn't switch them to scientific
options("scipen"=999,"digits"=4)

imputeloc = "/fh/fast/sun_w/licai/_tumor_eQTL/impute_v2.3.2_x86_64_dynamic/impute2"

root.folder = ".."
phased = "/fh/fast/sun_w/licai/eQTL_KIRC/data/phased"
geno = "/fh/fast/sun_w/licai/eQTL_KIRC/data/chr_flipped"
gwas_alignments = "/fh/fast/sun_w/licai/eQTL_KIRC/data/chr_flipped/precheck_log"
refinfo = "/fh/fast/sun_w/licai/_tumor_eQTL/1000GP_Phase3"
proj_output = proj_log = "../log/step6_imputation"
outdir = "/fh/fast/sun_w/licai/eQTL_KIRC/data/imputed"

system(sprintf("mkdir -p %s", outdir))
system(sprintf("chmod 777 %s", outdir))
system(sprintf("mkdir -p %s", proj_log))

#i = 22;j = 2
#this file will keep the split positions for each chromosome just in case they will be needed later
# Actually, this is almost an empty file?
write.table(data.frame("chr","ind"),sprintf("%s/step6_indicies.txt",proj_log),
            append=F,quote=F,row.names=F,col.names=F,sep="\t")
#we still use the same 1000G reference:
#all1000G="ALL_1000G_phase1integrated_v3" #common reference prefix
all1000G = "1000GP_Phase3" #common reference prefix
#processing  each chromosome, impute runs on 5MB chunks (no more than 7MB), so first split into chunks
ncpu = 1
get_block = function(str,split=" ",block=2){
  unlist(strsplit(str,split=split))[block]
}
# mem = 32000
for(i in 2:22){
  #i = 22
  #using the appropriate positions from both data and 1000G we find the range in which we will do imputation
  khf = sprintf("%s/phasing_chr%s.haps",phased,i)
  dati = read.table(khf,as.is=T)
  rng = range(dati[,3])
  
  #  con = gzfile(sprintf("%s/ALL_1000G_phase1integrated_v3_chr%s_impute.legend.gz",refinfo,i))
  chri = sprintf("chr%s",i)
  con = gzfile(sprintf("%s/%s_%s.legend.gz",refinfo,all1000G,chri))
  refi = readLines(con)
  close(con)
  #  rng2 = range(as.numeric(matrix(unlist(strsplit(refi[-1]," ")),nrow=11)[2,]))
  rng2 = range(as.numeric(sapply(refi[-1],get_block)))
  rng[1] = min(rng[1],rng2[1])
  rng[2] = max(rng[2],rng2[2])
  message(paste(rng2-rng,collapse=" "))
  
  mb = 1e6
  indj = sprintf("%se6",c(seq(floor(rng[1]/1e6),floor(rng2[2]/mb),by=5),ceiling(rng[2]/mb)))
  write.table(data.frame(i,paste(indj,collapse=";")),sprintf("%s/indicies.txt",proj_log),
              append=T,quote=F,row.names=F,col.names=F,sep="\t")
  #since we prephased using shapeit we need to use -use_prephased_g and -known_haps_g
  pref = sprintf("%s -use_prephased_g", imputeloc)
  m = sprintf("-m %s/genetic_map_chr%s_combined_b37.txt",refinfo,i)
  #  h = sprintf("-h %s/%s_chr%s_impute.hap.gz",refinfo,all1000G,i)
  #  l = sprintf("-l %s/%s_chr%s_impute.legend.gz ",refinfo,all1000G,i)
  h = sprintf("-h %s/%s_chr%s.hap.gz",refinfo,all1000G,i)
  l = sprintf("-l %s/%s_chr%s.legend.gz ",refinfo,all1000G,i)
  khg = sprintf("-known_haps_g %s",khf)
  oth = "-align_by_maf_g -Ne 20000 -seed 12345"
  
  #here we will run by the 5MB increment from the rng[1] until we exaust the chromosome
  for(j in 2:length(indj)){
    #j = length(indj)
    out = sprintf("-phase -o %s/phased_imputed_chr%s_%s_%s",outdir,i,indj[j-1],indj[j])    
    pos = sprintf("-int %s %s",indj[j-1],indj[j])
    
    log_output = sprintf("> %s/impute_chr%s.o 2> %s/impute_chr%s.e",
                         proj_log, i, proj_log, i)
    com = sprintf("time %s %s %s %s %s %s %s %s %s",pref,m,h,l,khg,oth,pos,out, log_output)
    #system(com)
    # com2 = sprintf("bsub -n %s -R 'span[hosts=1]' -M %s -q %s '%s'",ncpu,mem,queue,com)    
    #  qout = sprintf("--output=%s/rout/step5_out_%s.out",pipedir,i)    
    com2 = sprintf("sbatch --exclusive -t 0-5 --wrap='%s'",com)
    message(com2)
    system(com2)
  }
}

sessionInfo()
q("no")
