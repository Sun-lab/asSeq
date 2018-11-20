
# ----------------------------------------------------------------------
# gerate a list of heterzygous SNPs for each sample, which will be 
# used by function "extractAsReads"
#
# the output file is a tab-delimiated file with four columns, 
# chromosome, position, allele 1 and allele 2, without header. 
#
# the 1000G reference data were downloaded from 
# ftp://share.sph.umich.edu/1000genomes/fullProject/2012.02.14/
#   v2.20101123.autosomes.snps.tgz
#   v2.20101123.snpAnnotations.tgz
# ----------------------------------------------------------------------

setwd("~/1000G_data/")

for(chrId in 1:22){
  chr1  = sprintf("chr%d", chrId)

  message(chr1, " ", date())

  # ----------------------------------------------------------------------
  # read in the SNP annotation 
  # ----------------------------------------------------------------------

  annoFile = sprintf("anno_out/%s.annotation.txt", chr1)
  anno1    = read.table(annoFile, header=TRUE, sep="\t", 
              comment.char="", as.is=TRUE)

  snpFile = sprintf("snps/v2.20101123.%s.snps", chr1)
  snp1    = scan(snpFile, what=character())

  ## some of the markers in snpFile are indels
  ## try to exclude them for now

  wsnp = match(anno1$SNP, snp1)
  
  if(any(any(snp1[wsnp] != anno1$SNP))){
    stop("mismatch ;[\n")
  }

  # ----------------------------------------------------------------------
  # read the haplotypes imputed by MACH, sample by sample
  # ----------------------------------------------------------------------

  setwd("~/MACH_output_EA/")

  ff1 = sprintf("%s.hap", chr1)
  cmd = sprintf("wc -l %s", ff1)
  wc1 = system(cmd, intern=TRUE)
  wc1 = unlist(strsplit(wc1, split="\\s+"))[2]
  ## nn1 is the sample size
  nn1 = wc1/2 
  
  for(k in 1:nn1){
    message("  ", k, " ", date())
    
    # --------------------------------------------------------------------
    # use sed command to extrat two rows of the large file,
    # which are the hapltypes of an individual
    # --------------------------------------------------------------------
    
    ff2  = sprintf("%s.%d.%d.tmp", ff1, 2*k-1, 2*k)
    cmd1 = sprintf("sed -n '%d,%d p' %s > %s", 2*k-1, 2*k, ff1, ff2)
    system(cmd1)

    # --------------------------------------------------------------------
    # read in haplotypes and do some checking
    # --------------------------------------------------------------------
    
    dat1 = scan(ff2, what=character())
    
    cmd1 = sprintf("rm %s", ff2)
    system(cmd1)

    if(dat1[1] != dat1[4]){ stop("sample names do not match\n") }
    
    if(dat1[2] != "HAPLO1"){ 
      stop(sprintf("expect HAPLO1 but see %s\n", dat1[2])) 
    }
    
    if(dat1[5] != "HAPLO2"){ 
      stop(sprintf("expect HAPLO2 but see %s\n", dat1[5])) 
    }
    
    hap1 = unlist(strsplit(dat1[3], split=""))
    hap2 = unlist(strsplit(dat1[6], split=""))
    
    if(length(hap1) != length(snp1)){
      stop("length of haplotype 1 does not match annotation\n")
    }
    
    if(length(hap2) != length(snp1)){
      stop("length of haplotype 2 does not match annotation\n")
    }
    
    hap1 = toupper(hap1[wsnp])
    hap2 = toupper(hap2[wsnp])
    
    alleles = c("A", "C", "G", "T")
    
    # --------------------------------------------------------------------
    # extract SNPs with heterozygous genotypes
    # --------------------------------------------------------------------
    
    wdiff  = which(hap1 != hap2 & hap1 %in% alleles & hap2 %in% alleles)
    
    hap1A  = hap1[wdiff]
    hap1B  = hap2[wdiff]
    pos1   = anno1$coordinate[wdiff]
    
    hetDat = data.frame(chr=rep(chr1, length(pos1)), pos=pos1, hap1A, hap1B)
    hetDat = hetDat[order(pos1),]
    
    # --------------------------------------------------------------------
    # write out results
    # --------------------------------------------------------------------
    
    sam1   = unlist(strsplit(dat1[1], split="->"))[2]
    ff2    = sprintf("hetSNP/hetSNP_%s.txt", sam1)
    
    if(chr1 == "chr1"){
      app = FALSE
    }else{
      app = TRUE
    }
    
    write.table(hetDat, file = ff2, append = app, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)

  }
}
