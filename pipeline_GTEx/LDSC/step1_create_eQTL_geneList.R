# ------------------------------------------------------------
# create annotation file & create ld score for each cell type
# ------------------------------------------------------------

# create gene_List file 
# use most significant eQTL location as start end point
inputdir = "/fh/fast/sun_w/licai/_GTEx/wb_eQTL_output/"
genodir  = "/fh/fast/sun_w/licai/_GTEx/data_genotype_all/cnt/"
outdir   = "/fh/fast/sun_w/licai/LDSC/data_eQTL"
qval_cutoff = 0.05

append = F
for(chri in 1:22){
  geneInfo = read.table(sprintf("%sInfo_ex_chr%s.dat",genodir, chri), header = T,
                        as.is = T)
  snpInfo  = read.table(sprintf("%s/../genotype/snpInfo_chr%s.txt", genodir, chri),
                        header = T, as.is = T)
 
  Matrix_eQTL_file = sprintf("%swb_Matrix_eQTL_output_chr%s.txt", 
                             inputdir, chri) 
  Matrix_eQTL = read.table(Matrix_eQTL_file, header = T, as.is = F)
  number_of_indtest = read.table(sprintf("%s/wb_Matrix_eQTL_eigenMT_chr%s.txt", 
                                         inputdir, chri), header = T, as.is = F)
  for(geni in unique(Matrix_eQTL$gene)){
    mm2 = Matrix_eQTL[which(Matrix_eQTL$gene == geni), ]
    eQTL_snp = mm2$SNP[which.min(mm2$p.value)]
    geneInfo$start[geni] = geneInfo$end[geni] = eQTL_snp
  }
  geneInfo = geneInfo[complete.cases(geneInfo), ]
  qval = number_of_indtest$BF
  genes   =  geneInfo$gene[number_of_indtest$gene[which(qval<qval_cutoff)]]
  control =  setdiff(geneInfo$gene, genes)
  write.table(genes, file = sprintf("%s/MatrixEQTL.txt", outdir), 
              col.names = F, row.names = F,
              sep ='\t', quote =F, append = append)
  write.table(geneInfo, file = sprintf("%s/geneList_MatrixEQTL.txt", outdir), 
              col.names = F, row.names = F, sep ='\t', 
              quote =F, append = append)
  write.table(control, file = sprintf("%s/control_MatrixEQTL.txt", outdir),
              col.names = F, row.names = F, sep ='\t',
              quote =F, append = append)
  append = T
}

q("no")
