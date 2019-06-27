#module add tabix;module add gcc;module add gsl;module add r
#export PATH=/nas/longleaf/home/zhabotyn/progs/rasqual-master/bin/:$PATH

#For this script it is assumed that eigenMT is located in the next folder

normscore = function(vec) {
    len  = length(na.omit(vec))+1
    rank = rank(na.omit(vec))
    ties = (rank - floor(rank)) > 0
    new.vec = vec[!is.na(vec)] 
    new.vec[!ties]=qnorm(rank[!ties]/len)
    new.vec[ties] =0.5*(qnorm((rank[ties]+0.5)/len)+qnorm((rank[ties]-0.5)/len))
    vec[!is.na(vec)] = new.vec
    vec
}

norm = F
#norm = T
data.dir = "../data"
cnt.dir = sprintf("%s/cnt", data.dir)

bam.dir = sprintf("%s/bam", data.dir)
geno.dir = "../datagen" 
info.dir = "../inf"
eMT.dir = "../eMT"
vcf.dir = sprintf("%s/vcf", geno.dir)
gen.dir = sprintf("%s/gen", geno.dir)
if(!file.exists(gen.dir))dir.create(gen.dir)

snp = sprintf("%s/SNP_ex.vcf.gz", vcf.dir)

SNP_file_name = sprintf("%s/SNP.txt", gen.dir)
covariates_file_name = sprintf("%s/Covariates.txt", cnt.dir) 

#
#a.
#reformat from VCF to #alleles
#
vcf.hd = read.table(snp, as.is=T, nrow=4, comment.char=",", sep=",")
samp = unlist(strsplit(vcf.hd[[4,1]], split="\t"))[-(1:9)]
vcf.in = read.table(snp, as.is=T)
gen.in = matrix(substr(unlist(vcf.in[,-(1:9)]),1,3), ncol=length(samp))
gen.in[gen.in=="1|1"] = 2
gen.in[gen.in=="0|1"] = 1
gen.in[gen.in=="1|0"] = 1
gen.in[gen.in=="0|0"] = 0
colnames(gen.in) = samp
gen.in = cbind(snpid=vcf.in[,3], gen.in)
write.table(gen.in, SNP_file_name, row.names=F, col.names=T, quote=F, sep="\t")

#
#b.
#need to transpose covariate matrix
#
covar = read.table(sprintf("%s/samples.dat", cnt.dir), as.is=T)
covar[1,1] = "id"
write.table(t(covar[,1:5]), covariates_file_name, row.names=F, col.names=F, quote=F, sep="\t")

#
#c.
#reformat from counts to log10 expression
#
info = read.table(sprintf("%s/Info_ex.dat", cnt.dir), header=T, as.is=T)
expr = read.table(sprintf("%s/Tcnt_ex.dat", cnt.dir), as.is=T)
#will add normalization
depth = 10^as.numeric(covar[-1,2])
depth = depth/median(depth)
expr = log1p(as.matrix(expr)%*%diag(depth))
if(norm){
  output_file_name = "output_norm.txt"
  expression_file_name = sprintf("%s/GE_norm.dat", cnt.dir)
  expr = t(apply(expr, 1, normscore))
}else{
  output_file_name = "output_unnorm.txt"
  expression_file_name  = sprintf("%s/GE_unnorm.dat", cnt.dir)
}
colnames(expr) = samp
expr = cbind(geneid=info$id, expr)
write.table(expr, expression_file_name, row.names=F, col.names=T, quote=F, sep="\t")

snpspos = vcf.in[,2]
genepos = round((info[,3]+info[,4])/2)

library(MatrixEQTL)
cisDist = 2e5

useModel = modelLINEAR; 
pvOutputThreshold = 1e-1;
errorCovariance = numeric();
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name );


gene = SlicedData$new();
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile( expression_file_name )

cvrt = SlicedData$new();
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
cvrt$LoadFile( covariates_file_name )

me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    snpspos = snpspos,
    genepos = genepos, 
    cisDist = cisDist,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);


#adding eigenMT
#this will output eigenMT corrected p-values for minimum p-values from matrixEQTL output
eig.py = sprintf("%s/eigenMT.py", eMT.dir)
chri1 = genepos[1, 2]   #we process one chromosome at a time, so first value works for all
qtl1 = output_file_name #input to eigenMT is output from matrixEQTL
gen1 = SNP_file_name    #reusing SNP information from matrixEQTL input
genpos1 = snpspos_file_name#positions for the SNP_file_name, also used above
phepos1 = genepos_file_name#positions for expression file, also used above
out1 = gsub("output_", "eigenMT_", qtl1)
com = sprintf("python %s --CHROM %s --QTL %s --GEN %s --GENPOS %s --PHEPOS %s --OUT %s", 
                  eig.py, chri1, qtl1, gen1, genpos1, phepos1, out1)
system(com)



q("no")
