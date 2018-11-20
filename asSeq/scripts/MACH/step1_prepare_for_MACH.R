
setwd("/Volumes/Moon/TCGA_BRCA/genotype_normal/")

# ------------------------------------------------------------
# read data
# ------------------------------------------------------------

date()
dat = read.table("birdseed-v2.calls.txt", header=TRUE, na.strings="-1")
date()

dim(dat)
dat[1:2,1:5]

table(dat[,2])
table(is.na(dat[,2]))

# ------------------------------------------------------------
# read in the samples to be used
# ------------------------------------------------------------

sam = read.table("../info/brca_samples2use_after_qc_female_caucasian.txt", 
  header=TRUE, sep="\t", as.is=TRUE)
dim(sam)
sam[1,]

mat1 = match(sam$DNAnorml_arrayFile, names(dat))

if(any(is.na(mat1))){
  stop("some sampls in the sample file are not in the data\n")
}

dat = dat[,c(1,mat1)]
dim(dat)
dat[1:2,1:5]

if(any(sam$DNAnorml_arrayFile != names(dat)[-1])){
  stop("sample name mismatch\n")
}

# ------------------------------------------------------------
# read in SNP annotation 
# ------------------------------------------------------------

path = "/Users/suninsky/research/data/Affy/Affy6/anno/"
date()
info = read.csv(sprintf("%s/GenomeWideSNP_6.na32.annot.csv", path), 
                as.is=TRUE, comment.char="#")
date()

names(info)
table(info$Chromosome)
table(info$Strand)

dim(info)
info[1:2,1:4]
## some SNPs have more than more probe.Set
length(unique(info$Probe.Set.ID))
length(unique(info$dbSNP.RS.ID))

dim(dat)
dat[1:2,1:4]
length(unique(dat$probeset_id))

# ------------------------------------------------------------
# extract data 
# ------------------------------------------------------------

if(! all(dat$probeset_id %in% info$Probe.Set.ID)){
  stop("some probeset_id in data file are not recognized\n")
}

info = info[match(dat$probeset_id, info$Probe.Set.ID),]
dim(info)
all(dat$probeset_id == info$Probe.Set.ID)
table(info$Chromosome)

tb1 = table(info$dbSNP.RS.ID)
table(tb1)

snp2 = names(tb1)[tb1>1]
length(snp2)
snp2[1:5]

w2rm1 = which(info$Chromosome == "---" | info$Strand == "---")
w2rm2 = setdiff(which(info$dbSNP.RS.ID %in% snp2), match(snp2, info$dbSNP.RS.ID))

w2rm  = union(w2rm1, w2rm2)
length(w2rm1)
length(w2rm2)
length(w2rm)

if(length(w2rm) > 0){
  info = info[-w2rm,]
  dat  = dat[-w2rm,]
}

dim(info)
dim(dat)
table(info$Chromosome)
length(unique(info$Probe.Set.ID))
length(unique(info$dbSNP.RS.ID))

# ------------------------------------------------------------
# sort by SNP location 
# ------------------------------------------------------------

chrs = info$Chromosome
chrs[which(chrs=="X")]   = "23"
chrs[which(chrs=="Y")]   = "24"
chrs[which(chrs=="MT")]  = "25"

table(chrs)
chrs = as.numeric(chrs)
table(chrs)

pos  = as.numeric(info$Physical.Position)
summary(pos)

od   = order(chrs, pos)

info = info[od,]
dat  = dat[od,]

dim(info)
info[1:5,1:4]

dim(dat)
dat[1:2,1:2]

if(! all(dat$probeset_id == info$Probe.Set.ID)){
  stop("probeset_id mismatch\n")
}

# ------------------------------------------------------------
# remove the SNPs that have more than 5% missing values 
# ------------------------------------------------------------

pdata = data.matrix(dat[,-1])
dim(pdata)

nNA = rowSums(is.na(pdata))
summary(nNA)

w2rm = which(nNA > 0.05*ncol(pdata))
length(w2rm)

if(length(w2rm) > 0){
  dat   = dat[-w2rm,]
  info  = info[-w2rm,]
  pdata = pdata[-w2rm,]
}

dim(dat)
dim(info)
dim(pdata)

# ------------------------------------------------------------
# switch to forward strand
# ------------------------------------------------------------

table(info$Strand)
table(info$Strand.Versus.dbSNP)
table(info$Strand, info$Strand.Versus.dbSNP)

table(info$Allele.A)
table(info$Allele.B)

switchStrand <- function(xx){
	wA = which(xx == "A")
	wC = which(xx == "C")
	wG = which(xx == "G")
	wT = which(xx == "T")
  
  if(length(wA) > 0){ xx[wA] = "T" }
  if(length(wT) > 0){ xx[wT] = "A" }
  if(length(wC) > 0){ xx[wC] = "G" }
  if(length(wG) > 0){ xx[wG] = "C" }
  
  xx
}

wRev = which(info$Strand == "-")

info$Allele.A.Forward = info$Allele.A
info$Allele.A.Forward[wRev] = switchStrand(info$Allele.A[wRev])

info$Allele.B.Forward = info$Allele.B
info$Allele.B.Forward[wRev] = switchStrand(info$Allele.B[wRev])

table(info$Allele.A, info$Allele.A.Forward)
table(info$Allele.B, info$Allele.B.Forward)

# ------------------------------------------------------------
# further check MAF
# ------------------------------------------------------------

MAF = 0.5*rowSums(pdata, na.rm=TRUE)/rowSums((!is.na(pdata)), na.rm=TRUE)
MAF = pmin(MAF, 1 - MAF)
summary(MAF)

xx = strsplit(info$Minor.Allele.Frequency, split=" // ", fixed=TRUE)
table(sapply(xx, length))
xx = matrix(unlist(xx), ncol=5, byrow=TRUE)
xx = as.numeric(xx[,1])

cor(MAF, xx)^2

png("../figures/MAF_data_vs_affy_anno_Caucasian.png", width=5, height=5, 
  res=200, units="in")
par(mar=c(5,4,1,1))
smoothScatter(xx, MAF, xlab="MAF @ Caucasian", ylab="MAF @ Caucasian sample")
dev.off()

# ------------------------------------------------------------
# generate genotype data in terms of nucleotides for all 
# SNPs for MACH
# ------------------------------------------------------------

uchrs  = unique(info$Chromosome)

for(chr1 in uchrs){
  message("\n", chr1, " ", date())

  wchr   = which(info$Chromosome == chr1)
  
  if(any(diff(wchr) != 1)){
    stop("SNPs are not ordered\n")  
  }
    
  snps   = info$dbSNP.RS.ID[wchr]
  pdata1 = pdata[wchr,]
  info1  = info[wchr,]
  datPed = matrix("0", nrow=ncol(pdata1), ncol=2*nrow(pdata1))

  samID  = sam$DNAnorml_patientID
  famID  = 1:nrow(datPed)
  father = rep(0, nrow(datPed))
  mother = rep(0, nrow(datPed))
  sex    = rep("f", nrow(datPed))
  
  for(i in 1:nrow(pdata1)){
    
    if(i %% 10000 == 0){
      message("  ", i, " ", date())
    }
    
    w0 = which(pdata1[i,] == 0)  
    w1 = which(pdata1[i,] == 1)  
    w2 = which(pdata1[i,] == 2)
    
    if(length(w0)>0){
      datPed[w0,c(2*i-1, 2*i)] = info1$Allele.A.Forward[i]
    }
    
    if(length(w1)>0){
      datPed[w1,2*i-1] = info1$Allele.A.Forward[i]
      datPed[w1,2*i]   = info1$Allele.B.Forward[i]
    }
    
    if(length(w2)>0){
      datPed[w2,c(2*i-1, 2*i)] = info1$Allele.B.Forward[i]
    }
  }
  datPed = cbind(famID, samID, father, mother, sex, datPed)
  
  ff1 = sprintf("../data_EA/genotype_ped_chr%s_forward.txt", chr1)
  ff2 = sprintf("../data_EA/genotype_marker_chr%s.txt", chr1)

  write.table(datPed, file = ff1, append = FALSE, quote = FALSE, 
              sep = " ", row.names = FALSE, col.names = FALSE)
  
  datMak = cbind(rep("M", length(snps)), snps)

  write.table(datMak, file = ff2, append = FALSE, quote = FALSE, 
              sep = " ", row.names = FALSE, col.names = FALSE)

}


