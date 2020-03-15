# ml tabix;
# export CFLAGS="-I/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL/CLAPACK-3.2.1/INCLUDE -I/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL/CLAPACK-3.2.1/F2CLIBS -I/usr/include/gsl"
# export LDFLAGS="-L/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL/CLAPACK-3.2.1 -L/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL/CLAPACK-3.2.1/F2CLIBS -I/usr/include/gsl/lib"
# export PATH=/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL/rasqual/bin/:$PATH

# ------------------------------------------------------------
# submit RUASQUL by gene
# ------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
#args = c("22", "FALSE")
args
chri = as.numeric(args[1])
permute = as.logical(args[2])
chri
permute

if(permute){
  outF = "permuted_log"
}else{
  outF = "fixed_log"
}
outF
if(!dir.exists(outF)){
  dir.create(outF)
}
dir.create(paste0(outF, "/chr", chri))

#routrb  # permutation
#boutrb  # 
queue="general"
hours = 12
mem = 16

geno.dir = "../data_genotype_all" 
cnt.dir = sprintf("%s/cnt", geno.dir)


Ngene = system(sprintf("wc -l %s/Info_ex_chr%s.dat", cnt.dir, chri), intern = T)
Ngene = as.numeric(strsplit(Ngene, " ")[[1]][1]) - 1
for(gni in 1:Ngene){
  com = sprintf("R CMD BATCH '--args %s %s %s' step3_RASQUAL.R %s/chr%s/step3_RASQUAL_%s.Rout",
                gni, chri, permute, outF, chri,  gni)
  qout = sprintf("./%s/chr%s/out_%s.out",
                 outF,   chri,  gni)           
  com2 = sprintf("sbatch -t 0-%s --output %s --mem=%sg --wrap=\"%s\"", hours, qout, mem, com) 
  # --partition=restart --requeue
  message(com2)
  system(com2)     
}

q("no")
# grep Error step3*


# ------------------------------------------------------------
# submit RUASQUL by Chromosome
# ------------------------------------------------------------

permute = T
for(chri in 1:4){
  com = sprintf("sbatch -t 0-12 R CMD BATCH '--args %s %s' step3_submit.R step3_submit_permute%s_%s.Rout",
                chri, permute, permute, chri)
  message(com)
}


# ------------------------------------------------------------
# submit failed jobs RUASQUL results
# ------------------------------------------------------------
permute =T
out.dir = "/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL/out.dir"
per.dir = "/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL/per.dir"
if(permute){
  outF = "permuted_log"
  fileP = sprintf("%s/genes_to_check_per_%s.txt", 
                  dir1, chri)
}else{
  outF = "fixed_log"
  fileP = sprintf("%s/genes_to_check_%s.txt", 
                  dir1, chri)
}
dir1 = ifelse(permute, per.dir, out.dir)
for(chri in 22:1){
  #chri = 21
  if(T){
    rerun = read.table(file = fileP, as.is =T, header = T)
    rerun_gene = rerun$geneID[which(is.na(rerun$filesize) | rerun$filesize ==0)]
    for(gni in rerun_gene){
      com = sprintf("R CMD BATCH '--args %s %s %s' step3_RASQUAL.R %s/chr%s/step3_RASQUAL_%s.Rout",
                    gni, chri, permute, outF, chri,  gni)
      qout = sprintf("%s/chr%s/out_%s.out",
                     outF,   chri,  gni)           
      com2 = sprintf("sbatch --exclusive -t 2-0 --output %s --wrap=\"%s\"\n", 
                     qout, com) #  --partition=largenode -c 6 -t 2-0 --mem=368000
      cat(com2)
      # system(com2)     
    }
  }
}
  


# # ------------------------------------------------------------
# # check all RUASQUL results
# # ------------------------------------------------------------
# out.dir = "/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL/out.dir/detailed/"
# per.dir = "/fh/fast/sun_w/licai/_GTEx/R_batch5_wb_RASQUAL/per.dir/detailed/"
# cnt.dir = "../data_genotype_all/cnt"
# geno_dir = "../data_genotype_all/genotype"
# 
# flr = list.files(out.dir)
# flt = list.files(out.dir, pattern = "_time.txt")
# flr = setdiff(flr, flt)
# length(flt)
# flp = list.files(per.dir)
# flt = list.files(out.dir, pattern = "_time.txt")
# flp = setdiff(flp, flt)
# length(flp)
# 
# geneMissed = which(! genes %in% gsub(".txt","", flr))
# if(length(geneMissed) != 0){
#   for(gni in geneMissed){
#     com = sprintf("R CMD BATCH '--args %s %s' step3_RASQUAL.R routrb/chr%s/step3_RASQUAL_%s.Rout",
#                   gni, chri, chri,  gni)
#     qout = sprintf("./routrb/chr%s/out_%s.out",
#                    chri,  gni)
#     com2 = sprintf("sbatch -t 0-1 --output %s --mem=16g --wrap=\"%s\"",  qout, com)
#     message(com2)
#   }
# }
# 
# 
# # check the time between fixed and permuted 
# out.dir = "out.dir/"
# per.dir = "per.dir/"
# 
# timefileP = list.files(per.dir, pattern = "_time.txt")[sample(1:15000, 200)]
# table(timefileP %in% list.files(out.dir, pattern = "_time.txt"))
# timefileP = gsub("_time.txt","", timefileP)
# times = matrix(NA, nrow = length(timefileP), ncol = 3)
# rownames(times) = timefileP
# colnames(times) = c("per", "fixed", "diff")
# for(fli in timefileP){
#   # fli = timefileP[1]
#   times[fli, 1] = as.numeric(read.table(paste0(per.dir, fli, "_time.txt"), as.is =T)$V1)
#   times[fli, 2] = as.numeric(read.table(paste0(out.dir, fli, "_time.txt"), as.is =T)$V1)
#   times[fli, 3] = times[fli, 2] - times[fli, 1]
# }
# summary(times)

