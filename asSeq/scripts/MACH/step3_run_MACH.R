# -d: data file...marker names you included in your genotype files

#The data file should look like this:
#
# M marker1
# M marker2
# ...

# -p: pedigree file in Merlin format.

#The pedigree file should list one individual per row. Each row
#should start with an family id and individual id, followed by a
#father and mother id (which should both be 0, 'zero', since
#mach1 assumes individuals are unrelated), and sex. These initial
#columns are followed by a series of marker genotypes, each with
#two alleles. Alleles can be coded as 1, 2, 3, 4 or A, C, G, T.
#Missing genotypes can be represented as "M M" or "0 0"
#
#For example:
#
# FAM1001   ID1234  0   0   M   1 1   1 2   2 2
# FAM1002   ID1234  0   0   F   1 2   2 2   3 3
#
#Or:
#
#  FAM1001   ID1234  0   0   M  A A   A C   C C
#  FAM1002   ID1234  0   0   F  A C   C C   G G
#  
#-s: reference snp files. obtatined by splitting the 1000Genome snps
#-h: reference haplotype files. obtatined by splitting the 1000Genome hap files
#--r: how may round you want to run. larger rounds gives finer haplotypes. Recommended to be higher than 20.
#--autoflip: this is used when your input pedigree file are not necessarily on the same "+"/"-" strand. 
#            If you can preprare your pedigree file to be on the same strand, this should not be needed
#--states: this specifies the number of halotypes used in the imputation process, larger number 
#          leads to heavier computational cost

setwd("/lustre/scr/w/e/weisun/TCGA/MACH_sh")

sh   = "_run_MACH_EA_new.sh"
cat("", file=sh)

for(c1 in 1:22){
  
  ffs  = list.files(path="../data_EA", pattern=sprintf("chr%d.[[:digit:]]+.hap", c1))
  ks   = length(ffs)
  
  for(k1 in 1:ks){
    fHap = sprintf("chr%d.%d.hap",  c1, k1)
    fSNP = sprintf("chr%d.%d.snps", c1, k1)

    ww1 = which(ffs==fHap)
    if(length(ww1) != 1){ stop("ah... cannot find this file...\n") }
        
    cmd  = "../software/mach.1.0.18.Linux/executables/mach1"
    cmd  = sprintf("%s -d ../data_EA/genotype_marker_chr%d.txt", cmd, c1)
    cmd  = sprintf("%s -p ../data_EA/genotype_ped_chr%d_forward.txt", cmd, c1)
    cmd  = sprintf("%s -s ../data_EA/chr%d.%d.snps", cmd, c1, k1)
    cmd  = sprintf("%s -h ../data_EA/chr%d.%d.hap",  cmd, c1, k1)
    cmd  = sprintf("%s --r 100 --phase --autoflip --greedy --compact", cmd)
    cmd  = sprintf("%s --states 200 -o ../MACH_output_EA/chr%d_%d.out", cmd, c1, k1)
    cmd  = sprintf("%s --forceImputation > EA_chr%d_%d.log\n",   cmd, c1, k1)
    
    sh1  = sprintf("EA_chr%d_%d.sh", c1, k1)
    cat(cmd, file=sh1)
    
    cat(sprintf("bsub -q week bash %s\n", sh1), file=sh, append=TRUE)
  }
  
}
