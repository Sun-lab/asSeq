#!/bin/bash
ct=$1
cd ../ldsc/
module load anaconda2
source activate ldsc
ml bedtools
for chr in {1..22}
    do
    python make_annot.py \
    --gene-set-file ../data_eQTL/control_${ct}.txt \
    --gene-coord-file ../data_eQTL/geneList_${ct}.txt \
    --windowsize 100000 \
    --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
    --annot-file ../data_eQTL/control_${ct}.${chr}.annot.gz
    done

for chr in {1..22}
    do
    python ../ldsc/ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
    --ld-wind-cm 1 --annot ../data_eQTL/control_${ct}.${chr}.annot.gz --thin-annot \
    --out ../data_eQTL/control_${ct}.${chr} \
    --print-snps hapmap3_snps/hm.${chr}.snp
    done
source deactivate
