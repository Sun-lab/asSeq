
# add header (GENE  CHR START END) using vim
# ------------------------------------------------------------
# creat annotation file and calculate ld score
# ------------------------------------------------------------

file_ct = sprintf("anno_ldsc.sh")
cat("#!/bin/bash\n", file= file_ct)
cat("ct=$1\n", file= file_ct, append = T)
cat("chr=$2\n", file= file_ct, append = T)
cat("control=$3\n", file= file_ct, append = T)
cat("cd ../ldsc/\n", file= file_ct, append = T)
cat("echo 'cell type:' ${ct}\n", file= file_ct, append = T)
cat("echo 'chr:' ${chr}\n", file= file_ct, append = T)
cat("echo 'control': ${control}\n", file= file_ct, append = T)

cat(sprintf('module load anaconda2\n'), file= file_ct, append = T)
cat(sprintf('source activate ldsc\n'), file= file_ct, append = T)
cat(sprintf('ml bedtools\n'), file= file_ct, append = T)

cat("python make_annot.py \\
    --gene-set-file ../data_eQTL/${control}${ct}.txt \\
    --gene-coord-file ../data_eQTL/geneList_${ct}.txt \\
    --windowsize 100000 \\
    --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \\
    --annot-file ../data_eQTL/${control}${ct}.${chr}.annot.gz  
    \n\n",  file=file_ct, append = T)
cat("python ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \\
    --ld-wind-cm 1 --annot ../data_eQTL/${control}${ct}.${chr}.annot.gz \\
    --thin-annot \\
    --out ../data_eQTL/${control}${ct}.${chr} \\
    --print-snps hapmap3_snps/hm.${chr}.snp \n
    ",  file= file_ct, append = T)
cat("\nsource deactivate\n", file= file_ct, append = T)

file_ct 
ct='MatrixEQTL'
control="control_"#"" #
for(chri in 22:1){
  message(sprintf("sbatch %s %s %s %s", file_ct, ct, chri, control))
}

# ------------------------------------------------------------
# analysis
# ------------------------------------------------------------

cat("python ldsc.py \\
--h2-cts ../data/clozuk_pgc2.meta.sumstats_with_rsid.gz \\
--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \\
--out ../data/clozuk_pgc2.meta_ct6 \\
--ref-ld-chr-cts ../data/ct6.ldcts \\
--w-ld-chr weights_hm3_no_hla/weights.")
# ---

cat("#!/bin/bash\n", file= "ldscore_analysis.sh")
cat("ct=$1\n", file= "ldscore_analysis.sh", append = T)
cat("echo 'cell type:' ${ct}\n", file= "ldscore_analysis.sh", append = T)
cat("cd ../ldsc/\n", file= "ldscore_analysis.sh", append = T)
cat(sprintf('module load anaconda2\n'), file= "ldscore_analysis.sh", append = T)
cat(sprintf('source activate ldsc\n'), file= "ldscore_analysis.sh", append = T)
cat(sprintf('ml bedtools\n'), file= "ldscore_analysis.sh", append = T)
cat("python ldsc.py \\
	--h2 ../data/clozuk_pgc2.meta.sumstats_with_rsid.gz \\
	--w-ld-chr weights_hm3_no_hla/weights. \\
	--ref-ld-chr ../data_eQTL/${ct}.,../data_eQTL/control_${ct}. \\
	--frqfile-chr 1000G.mac5eur. \\
	--out ../output_eQTL/${ct} \\
	--print-coefficients",  file="ldscore_analysis.sh", append = T)
cat("\nsource deactivate\n", file="ldscore_analysis.sh", append = T)

q('no')