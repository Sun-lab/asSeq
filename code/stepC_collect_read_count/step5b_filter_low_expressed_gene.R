setwd("/fh/fast/sun_w/licai/eQTL_KIRC/data_cloud/step5_TReC_per_gene_filterIt")
TReC_files = list.files(pattern = 'gene_level_counts_filterIt_total.txt')
# list.files(pattern = 'TReC_liberal_filterIt_total.txt')

nms = system(sprintf('cut -f1 %s', TReC_files[1]), intern = T)
nms = nms[-1]
nms[1:5]
length(nms)
geneE = matrix(NA, nrow = length(nms), ncol = length(TReC_files))
colnames(geneE) = substr(TReC_files,1,16)
rownames(geneE) = nms
for(i in TReC_files){
  # i = TReC_files[1]
  ffi = read.table(i, header = T, as.is = T, check.names = F)
  head(ffi)
  if(all(nms== rownames(ffi))){
        geneE[,colnames(ffi)] = ffi[,1]
  }
  else
    geneE[rownames(ffi),colnames(ffi)] = ffi[,1]
}

dim(geneE)
geneE[1:5,1:5]

# ----------------------------------------------------------------------------
# filter out low expressed genes 
# gene with 75% of samples having TReC less than 20 was removed 
# ----------------------------------------------------------------------------

r75  = apply(geneE, 1, quantile, probs = 0.75)
w2kp = which(r75 >= 20)
length(w2kp)

geneE = geneE[w2kp, ]
dim(geneE)

write.table(geneE, file = sprintf("gene_level_counts_filterIt_total_after_filter_low_expressed_gene.txt"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", row.names = TRUE,col.names = TRUE)
