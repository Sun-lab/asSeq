# -----------------------------------------------------------------------
# handle chromosome 1 only
# -----------------------------------------------------------------------

chrId = 1

# -----------------------------------------------------------------------
# prepare position files that will be used in ligage the MACH results
# -----------------------------------------------------------------------

setwd("~/lustre/TCGA/MACH_output_EA")

pos1 = 1:110001
write.table(pos1, file = "pos1", quote =F, col.names = F, row.names =F)

for (p in 2:34){
  print(c(((p-1)*100000-10000 +1), (p*100000+10000 + 1)))
  pos = ((p-1)*100000-10000 +1):(p*100000+10000 + 1)
  write.table(pos, file = sprintf("pos%d", p),quote =F, col.names = F, row.names =F)
}

# ----------------------------------------------------------------------
# read in the files names of MACH output
# ----------------------------------------------------------------------

setwd("~/lustre/TCGA/MACH_output_EA")

mcdir = "../data_EA/"

ffs = list.files(pattern="out.gz")
ff1 = gsub(".out.gz", "", ffs, fixed=TRUE)
ff1 = gsub("chr", "", ff1, fixed=TRUE)

chr = matrix(unlist(strsplit(ff1,"_")), ncol=2, byrow = TRUE)
chr = data.frame(ffs, as.numeric(chr[,1]), as.numeric(chr[,2]))
names(chr) = c("file", "chr", "part")

od = order(chr$chr, chr$part)
chr = chr[od,]

chri = chr[chr$chr == chrId,]

if (any(chri$part != seq(nrow(chri)))){
  stop("Not all parts are finished")
}

wp = 1
cmd1 = "~/luster/TCGA/1000G/software/ligateHaplotypes_V004/ligateHaplotypes.pl"

# -----------------------------------------------------------------------
# Deal with the first split file. This is the command we want to generate
#   ligateHaplotypes.pl -ah chr22_1.out.gz -ap  pos1 -bh chr22_2.out.gz \
#   -bp pos2 -o chr22.cum1_2
# -----------------------------------------------------------------------

if(wp == 1){
  cmd = sprintf("%s -ah %s -ap pos1 -bh %s -bp pos2", cmd1, chri$file[wp], chri$file[wp+1])
  cum = sprintf("chr%d.cum_%d_%d", i, chri$part[wp], chri$part[wp+1])
  cmd = sprintf("%s -o %s > %s.log", cmd, cum, cum)
  system(cmd)
  wp = wp + 1
  
  message(wp, " ", date())
}

# -----------------------------------------------------------------------
# Deal with the split files in the middle
# -----------------------------------------------------------------------

while(wp < (max(chri$part)-1)){
  cum = sprintf("chr%d.cum_%d_%d", i, chri$part[wp-1], chri$part[wp])
  cmd = sprintf("%s -ah %s.hap -ap %s.pos -bh %s -bp pos%d", cmd1, cum, cum, chri$file[wp+1], wp+1)
  cum = sprintf("chr%d.cum_%d_%d", i, chri$part[wp], chri$part[wp+1])
  cmd = sprintf("%s -o %s > %s.log", cmd, cum, cum)
  system(cmd)
  wp = wp + 1
  
  message(wp, " ", date())
}

# -----------------------------------------------------------------------
# Deal with the last split file
# -----------------------------------------------------------------------

if(wp == (max(chri$part)-1)){
# Creat new position files
  nsnp = system(sprintf("wc -l %schr%d.%d.snps", mcdir, i, wp+1), intern=TRUE)
  nsnp = as.numeric(unlist(strsplit(nsnp, split=" "))[1])
  pos1 = (wp*100000-10000 +1)
  pos2 = pos1 + nsnp - 1

  fpos = sprintf("chr%d_%d_end.pos", i, wp+1)
  write.table(pos1:pos2, file=fpos, quote=F, col.names=F, row.names=F)

  cum = sprintf("chr%d.cum_%d_%d", i, chri$part[wp-1], chri$part[wp])
  cmd = sprintf("%s -ah %s.hap -ap %s.pos -bh %s -bp %s", cmd1, cum, cum, chri$file[wp+1], fpos)
  cum = sprintf("chr%d.cum_%d_%d", i, chri$part[wp], chri$part[wp+1])
  cmd = sprintf("%s -o chr%d > %s.log", cmd, i, cum)
  system(cmd)
  wp = wp + 1

  message(wp, " ", date())
}

