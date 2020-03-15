#
#c
#create other required files if they don't exist yet
#

# ----------------------------------------------------------------------
# TReC
# ----------------------------------------------------------------------

step5    =  "../R_batch4_whole_blood/stepC/output/step5_TReC_per_gene_filterIt/"
geno.dir = "../data_genotype_all" 
cnt.dir = sprintf("%s/cnt", geno.dir)
vcfr = sprintf("%s/vcfr/SNP_chr22.vcf.gz", geno.dir)

vcfr.hd = read.table(vcfr,as.is=T, nrow = 4, comment.char=",", sep=",")
samples = unlist(strsplit(vcfr.hd[4,1], split="\t"))[-(1:9)]
samples = gsub("_", "", substr(samples, 1, 10))
length(samples)


xcnt = sprintf("%s/whole_blood_covs_eQTL.txt",geno.dir)
xcnt = read.table(xcnt, as.is=T, header=T)
xcnt = xcnt[samples,, drop=F]
xcnt = data.frame("ID" = rownames(xcnt), xcnt)
xcnt$rd = log(xcnt$rd)
head(xcnt)


for(chri in 1:22){
  # chri = 22
  
  step5_file = sprintf("gene_level_counts_filter_out_low_expressed_chr%s.txt",chri)
  totc = sprintf("%s/%s", step5, step5_file)
  totc = read.table(totc, as.is=T, check.names = F)
  totc = totc[, samples]
  totc[1:5,1:5]
  totb = sprintf("%s/Tcnt_ex_chr%s.bin", cnt.dir, chri)
  kbin = sprintf("%s/K_ex_chr%s.bin", cnt.dir, chri)
  xbin = sprintf("%s/X_ex.bin", cnt.dir)
  
  geneInfo_chr = read.csv(sprintf("%s/Info_ex_chr%s.dat",cnt.dir, chri), 
                          header = T, as.is = T, sep ='\t')
  
  if(file.exists(totb)){
    con = file(totb, "wb")
    writeBin(as.double(t(totc)), con=con, double())
    close(con)
  }
  if(file.exists(kbin)){
    Ks = exp(xcnt[,2])
    Ks = matrix(rep(Ks/mean(Ks),each=nrow(totc)), nrow=nrow(totc))
    Ks[1:5,1:5]
    con = file(kbin,"wb")
    writeBin(as.double(c(t(Ks))), con=con, double())
    close(con)
  }
  
}
if(file.exists(xbin)){
  con = file(xbin, "wb")
  writeBin(as.double(c(as.matrix(xcnt[,c(2:7),drop=F]))), con=con, double())
  close(con)
}

kin = readBin(kbin, n=200, double())
kin
tin = readBin(totb, n=200, double())
tin
xin = readBin(xbin, n=200, double())
xin

