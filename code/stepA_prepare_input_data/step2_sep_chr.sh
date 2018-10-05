#!/bin/bash
#SBATCH --job-name=sep_chr
cd ../../data
mkdir -p chr
ml plink/1.90

for chr in $(seq 1 22);
do
plink --file geno \
--chr $chr \
--recode \
--out ./chr/geno.chr$chr ;
done
