# getting a chromosome level summary of how many snps are to be removed
# due to missing values or to strand issue

geno =  "../../data/chr"
gwas_alignments = "../../data/chr/precheck_log"

for (suff in c("")) {
  chrs = seq(22)
  summ = matrix(NA, nrow = 22, ncol = 4)
  colnames(summ) = c("total", "not.comp", "strand", "missing")
  rownames(summ) = paste("chr", 1:22, sep = "")
  for (i in chrs) {
    genodat = try(system(sprintf("cat %s/geno.chr%s.map | wc -l", geno, i), 
                         intern =TRUE))
    summ[i, 1] = as.numeric(genodat)
    dat = read.table(sprintf("%s/gwas.alignments_chr%d.snp.strand",gwas_alignments, i),
                     header = F, as.is = T, sep = ',')
    dat = dat[-1, ]
    summ[i, 2] = length(dat)
    tbl = table(matrix(unlist(strsplit(dat, "\t")), ncol = length(dat))[1,])
    summ[i, 3] = tbl[["Strand"]]
    summ[i, 4] = tbl[["Missing"]]
    message(i)
  }
  message(suff)
}

summ

sessionInfo()
q("no")
