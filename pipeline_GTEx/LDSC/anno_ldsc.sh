#!/bin/bash
ct=$1
chr=$2
control=$3
cd ../ldsc/
echo 'cell type:' ${ct}
echo 'chr:' ${chr}
echo 'control': ${control}
module load anaconda2
source activate ldsc
ml bedtools
python make_annot.py \
    --gene-set-file ../data_eQTL/${control}${ct}.txt \
    --gene-coord-file ../data_eQTL/geneList_${ct}.txt \
    --windowsize 100000 \
    --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
    --annot-file ../data_eQTL/${control}${ct}.${chr}.annot.gz  
    

python ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
    --ld-wind-cm 1 --annot ../data_eQTL/${control}${ct}.${chr}.annot.gz \
    --thin-annot \
    --out ../data_eQTL/${control}${ct}.${chr} \
    --print-snps hapmap3_snps/hm.${chr}.snp 

    
source deactivate
