args = commandArgs(TRUE)
#args = c("22")
chri = args[1]

#ensure that when we write positions to a file R doesn't switch them to scientific
options("scipen"=999,"digits"=4)

# For dosage, I mean AA = 0, AB=1 and BB =2. The genotype output from IMPUTE2 has posterior prob of the three genotypes and I just calculate the expectation.
# For file format, see http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html
# File location example:
# /lustre/scr/c/h/chongjin/COAD/stepA_phasing_and_imputation/output/step8_submit_impute_kill/imputed/phased_imputed_chr1_10e6_15e6
# Example:  Suppose you want to create a genotype for 2 individuals at 5 SNPs whose genotypes are
# 
# SNP 1 : AA AA
# SNP 2 : GG GT
# SNP 3 : CC CT
# SNP 4 : CT CT
# SNP 5 : AG GG
# 
# The correct genotype file would be
# 
# SNP1 rs1 1000 A C 1 0 0 1 0 0
# SNP2 rs2 2000 G T 1 0 0 0 1 0
# SNP3 rs3 3000 C T 1 0 0 0 1 0
# SNP4 rs4 4000 C T 0 1 0 0 1 0
# SNP5 rs5 5000 A G 0 1 0 0 0 1
# 
# So, at SNP3 the two alleles are C and T so the set of 3 probabilities for each indvidual correspond to the genotypes CC, CT and TT respectively.

root.folder = "."
phased = "/fh/fast/sun_w/licai/eQTL_KIRC/data/phased"
outdir = "/fh/fast/sun_w/licai/eQTL_KIRC/data_genotype"
if(!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
inputdir = "/fh/fast/sun_w/licai/eQTL_KIRC/data/imputed"
setwd(inputdir)

#the list of samples for which we will get snps
samples = read.table(sprintf("%s/phasing_chr%s.sample",phased,chri),header=T,as.is=T)
samples = samples[-1,2]
#define snps (I'll exclude indels with matching to the base
bases=c("A","C","G","T")

#get all the appropriate chunks for the given chromosome
fls = list.files(pattern=sprintf("chr%s_",chri))
to.rm = union(union(union(grep(".txt",fls),grep("info",fls)),grep("warnings",fls)),grep("summary",fls))
fls = fls[-to.rm]
hps = fls[grep("_haps",fls)]
alp = fls[grep("_allele_probs",fls)]
fls = setdiff(fls,union(hps,alp))
fls
hps
alp

append=F
#
j = 1;off=5
probs = c(.9,.95,.99,.995,.999)

sample_count = length(samples)
summ_het = matrix(0,nrow=sample_count,ncol=length(probs))
rownames(summ_het) = samples
colnames(summ_het)=sprintf("p%s",probs)
tot = rep(0,sample_count)

#phe
for(j in 1:length(fls)){   
  #j=match("phased_imputed_chr2_85e6_90e6",fls)+1
  #j = 1
  
  flj = fls[j]
  snpj = read.table(flj,as.is=T)
  #hpj = read.table(hps[j],as.is=T)
  #alj = read.table(alp[j],as.is=T)
  # snpj$V4 = as.character(snpj$V4)
  # snpj$V5 = as.character(snpj$V5)
  
  rm.indels =  which(snpj$V4%in%bases & snpj$V5%in%bases)
  snpj = snpj[rm.indels,]
  #alj = alj[rm.indels,]
  #hpj = hpj[rm.indels,]
  #dim(hpj)
  #dim(alj)
  dim(snpj)
  
  for(ind in 1:sample_count){
    #ind = 1
    #ind=ind+1
    nmind=sprintf("%s",samples[ind])
    out = sprintf("%s/%s",outdir,nmind)
    if(!file.exists(out))dir.create(out, recursive=TRUE)
    indloc = off +(ind-1)*3+1
    snpind = snpj[,c(1:5,indloc,indloc+1,indloc+2)]
    snpind[,1] = sprintf("chr%s",chri)
    # choose the most likely genotype if its probability is larger than 0.8
    snpind[,6] = apply(data.matrix(snpind[,6:8]), 1, function(x){
                                                      geno = which(x >= 0.8)
                                                      ifelse(length(geno) != 0L, geno-1, NA)})
    outfile = sprintf("%s/chr%s.txt",out,chri)
    app=file.exists(outfile)
    write.table(snpind[,c(1,3,4,5,2,6)],outfile,row.names=F,col.names=F,quote=F,append=app)
    #}
  }
  
  message(flj)
}

#summary of 
#outfile = sprintf("%s/chr_%s_num_snps_1.txt",outdir,chri)
#write.table(cbind(summ_het,tot),outfile,row.names=T,col.names=T,quote=F)

sessionInfo()
q("no")
