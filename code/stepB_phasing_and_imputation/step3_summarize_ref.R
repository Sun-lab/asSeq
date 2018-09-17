# getting a chromosome level summary of how many snps are to be removed 
# due to missing values or to strand issue

geno =  "/fh/fast/sun_w/licai/eQTL_KIRC/data/chr"
gwas_alignments = "/fh/fast/sun_w/licai/eQTL_KIRC/data/chr/precheck_log"
cancer_type = 'KIRC'
for(suff in c("")){
  
  chrs = 1:22
  summ = matrix(NA,nrow=length(chrs),ncol=4)
  colnames(summ) = c("total","not.comp","strand","missing")
  rownames(summ) = paste("chr",chrs,sep="")
  for(i in chrs){
    #i = 16
    #genodat = read.table(sprintf("%s/geno_chr%d.map",geno,i),as.is=T)
    #summ[i,1] = nrow(genodat)
    
    genodat = try(system(sprintf("cat %s/%s.chr%s.map | wc -l", geno,cancer_type, i), intern=TRUE))
    summ[i,1] = as.numeric(genodat)
    
    dat = read.table(sprintf("%s/gwas.alignments_chr%d.snp.strand",
                             gwas_alignments,i),header=F,as.is=T, sep =',')
    dat = dat[-1, ]
    summ[i,2] = length(dat)
    tbl = table(matrix(unlist(strsplit(dat,"\t")),ncol=length(dat))[1,])
    summ[i,3] = tbl[["Strand"]]
    summ[i,4] = tbl[["Missing"]]
    message(i)
  }
  message(suff)
}

summ

sessionInfo()
q("no")
