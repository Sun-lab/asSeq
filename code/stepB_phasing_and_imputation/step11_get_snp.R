args = commandArgs(TRUE)
#args = c("22")
chri = args[1]

#ensure that when we write positions to a file R doesn't switch them to scientific
options("scipen"=999,"digits"=4)
#note that we also need to convert these snps to build 38
#I have converted by liftover positions of those snps
#also these snps build37 (under R302) build38 (under R312)
#source("http://bioconductor.org/biocLite.R");biocLite("SNPlocs.Hsapiens.dbSNP.20120608")
#library(SNPlocs.Hsapiens.dbSNP.20120608)
#library(SNPlocs.Hsapiens.dbSNP141.GRCh38)

root.folder = "."
phased = "/fh/fast/sun_w/licai/eQTL_KIRC/data/phased"
outdir = "/fh/fast/sun_w/licai/eQTL_KIRC/data_snp"
if(!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
inputdir ="/fh/fast/sun_w/licai/eQTL_KIRC/data/imputed"
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

#phe
for(j in 1:length(fls)){   
  #j=match("phased_imputed_chr2_85e6_90e6",fls)+1
  #j = 1
  
  flj = fls[j]
  snpj = read.table(flj,as.is=T)
  hpj = read.table(hps[j],as.is=T)
  alj = read.table(alp[j],as.is=T)
  infoj = read.table(info[j], as.is = T, header = T)
  snpj$V4 = as.character(snpj$V4)
  snpj$V5 = as.character(snpj$V5)
  
  rm.indels =  which(snpj$V4%in%bases & snpj$V5%in%bases)
  snpj = snpj[rm.indels,]
  alj = alj[rm.indels,]
  hpj = hpj[rm.indels,]
  infoj = infoj[rm.indels,]
  dim(hpj)
  dim(infoj)
  dim(alj)
  dim(snpj)
  
  # remove R2 < 0.3  "'info' is similar to the r-squared metrics reported by other programs like MaCH and Beagle"
  R2kp = which(infoj$info > 0.3)
  snpj = snpj[R2kp,]
  alj = alj[R2kp,]
  hpj = hpj[R2kp,]
  infoj = infoj[R2kp,]
  dim(hpj)
  dim(infoj)
  dim(alj)
  dim(snpj)
  
  
  for(ind in 1:sample_count){
    #ind = 1
    #ind=ind+1
    nmind=sprintf("%s",samples[ind])
    out = sprintf("%s/%s",outdir,nmind)
    if(!file.exists(out))dir.create(out)
    hetloc = off+(ind-1)*3+2
    indloc = off +(ind-1)*2+1
    for(probk in 1:length(probs)){
      summ_het[ind,probk] = summ_het[ind,probk] + sum(snpj[,hetloc] >= probs[probk])
    }
    tot[ind] = tot[ind] + nrow(snpj)
    #I'll use only the snps with probability>0.999
    flag = snpj[,hetloc]>probs[5]
    if(sum(flag)>0){
      snpind = snpj[flag,]
      alind = alj[flag,c(1:5,indloc,indloc+1)]
      hpind = hpj[flag,c(1:5,indloc,indloc+1)]
      infoind = infoj[flag, ]
      flip = (hpind[,6]==1)
      hpind[,1] = sprintf("chr%s",chri)
      hpind[flip,4:5] = hpind[flip,5:4]     
      outfile = sprintf("%s/chr%s.txt",out,chri)
      app=file.exists(outfile)
      write.table(hpind[,c(1,3,4,5,2)],outfile,row.names=F,col.names=F,quote=F,append=app)
    }
  }
  
  message(flj)
}

#summary of 
outfile = sprintf("%s/chr_%s_num_snps_1.txt",outdir,chri)
write.table(cbind(summ_het,tot),outfile,row.names=T,col.names=T,quote=F)

sessionInfo()
q("no")
