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
phased = "../../data/phased"
outdir = "../../data/data_genotype"
out = sprintf("%s",outdir)
if(!file.exists(out))dir.create(out, recursive=TRUE)
outfile = sprintf("%s/chr%s.txt",out,chri)
if( file.exists(outfile) ) file.remove(outfile)
library(data.table)

# Load cpp function
library(Rcpp); library(RcppArmadillo)
cpp_fn = file.path("../../code/stepB_phasing_and_imputation","dose2geno.cpp")
sourceCpp(cpp_fn,showOutput = TRUE)

if(!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
inputdir = "../../data/imputed"
setwd(inputdir)

#the list of samples for which we will get snps
samples = read.table(sprintf("%s/phasing_chr%s.sample",phased,chri),header=T,as.is=T)
samples = samples[-1,2]
#define snps (I'll exclude indels with matching to the base
bases=c("A","C","G","T")

#get all the appropriate chunks for the given chromosome
fls = list.files(pattern=sprintf("chr%s_",chri))
to.rm = union(union(union(grep(".txt",fls),grep("info_by_sample",fls)),grep("warnings",fls)),grep("summary",fls))
fls = fls[-to.rm]
hps = fls[grep("_haps",fls)]
alp = fls[grep("_allele_probs",fls)]
info = fls[grep("info$",fls)]
fls = setdiff(fls,union(union(hps,alp), info))
fls
info
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


for(jj in 1:length(fls)){   
	#j=match("phased_imputed_chr2_85e6_90e6",fls)+1
	#j = 1

	# flj = fls[j]
	# snpj = read.table(flj,as.is=T)
	# infoj = read.table(info[j], as.is = T, header = T)
	# hpj = read.table(hps[j],as.is=T)
	# alj = read.table(alp[j],as.is=T)
	# snpj$V4 = as.character(snpj$V4)
	# snpj$V5 = as.character(snpj$V5)

	cat(paste0("Progress: ",jj," out of ",length(fls)," - Reading in ",fls[jj],"...\n"))
	snpj_fn 	= file.path(inputdir,fls[jj])
	infoj_fn 	= file.path(inputdir,paste0(fls[jj],"_info"))
	hpsj_fn 	= file.path(inputdir,paste0(fls[jj],"_haps"))

	# Import dosages, haplotypes, QC metrics
	cat("\tReading in dosages, haplotypes, QC metrics...\n")
	snpj = data.table::fread(snpj_fn,sep = " ",
		header = FALSE,stringsAsFactors = FALSE,
		quote = "",data.table = FALSE,showProgress = FALSE)
	infoj = data.table::fread(infoj_fn,sep = " ",
		header = TRUE,stringsAsFactors = FALSE,
		quote = "",data.table = FALSE,showProgress = FALSE)
	hpj = data.table::fread(hpsj_fn,sep = " ",
		header = FALSE,stringsAsFactors = FALSE,
		quote = "",data.table = FALSE,showProgress = FALSE)
  
	rm.indels =  which(snpj$V4%in%bases & snpj$V5%in%bases)
	snpj = snpj[rm.indels,]
	infoj = infoj[rm.indels,]
	#alj = alj[rm.indels,]
	hpj = hpj[rm.indels,]
	#dim(hpj)
	#dim(alj)
	dim(snpj)

	R2kp = which(infoj$info > 0.3)
	MAF2kp = which(infoj$exp_freq_a1 > 0.05 & infoj$exp_freq_a1 <0.95)
	rm.MAF.R2 = intersect(R2kp, MAF2kp)
	snpj = snpj[rm.MAF.R2,]
	infoj = infoj[rm.MAF.R2,]
	hpj = hpj[rm.MAF.R2,]
	dim(snpj)
	dim(hpj)
	dim(infoj)

	genotype_mat = matrix(NA, nrow=nrow(hpj), ncol=sample_count)
	# colnames(genotype_mat) = samples
	  
	for(ind in 1:sample_count){
		#ind = 1
		#ind=ind+1
		# nmind=sprintf("%s",samples[ind])
		# out = sprintf("%s/%s",outdir,nmind)

		haploc = off+(ind-1)*2+1
		indloc = off +(ind-1)*3+1
		snpind = snpj[,c(1:5,indloc,indloc+1,indloc+2)]
		hapind = hpj[,c(1:5,haploc, haploc+1)]
		snpind[,1] = sprintf("chr%s",chri)
		# # choose the most likely genotype if its probability is larger than 0.8
		# snpind[,6] = apply(data.matrix(snpind[,6:8]), 1, function(x){
		#           geno = which(x >= 0.8)
		#           ifelse(length(geno) != 0L, geno-1, NA)})

		# choose the most likely genotype if its probability is larger than 0.8. and 
		# code it to 1:AA, 2:AB, 3:BA, 4:BB.  (A-V4, B-V5)
		# create a inputed genotype file for all samples per chromosom
		# geno_ind = apply(data.matrix(snpind[,6:8]), 1,function(x){
		#	geno = which(x >= 0.8)
		#	ifelse(length(geno)!=0, geno-1, NA)
		#	})
		geno_ind = Rcpp_dose2geno(
			doses = as.matrix(snpj[,snp_index]),
			geno_thres = 0.8)[,1]
		geno_ind[which(geno_ind == 5)] = NA
		geno_ind[which(geno_ind == 2)] = 3
		geno_hetind = which(geno_ind == 1)
		head(hapind)
		flip = hapind[geno_hetind,6]==1
		geno_ind[geno_hetind[flip]] = 2
		table(geno_ind)
		# outfile = sprintf("%s/chr%s.txt",out,chri)
		# app=file.exists(outfile)
		#}
		genotype_mat[,ind] = geno_ind
	}
	app=file.exists(outfile)
	genotype_mat = cbind(snpind[,c(1,3,4,5,2)],genotype_mat)
	write.table(genotype_mat, outfile,row.names=F,col.names=F,quote=F,append=app)
	message(flj)
}

#summary of 
#outfile = sprintf("%s/chr_%s_num_snps_1.txt",outdir,chri)
#write.table(cbind(summ_het,tot),outfile,row.names=T,col.names=T,quote=F)

sessionInfo()
q("no")
