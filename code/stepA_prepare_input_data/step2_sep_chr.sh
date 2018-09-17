#!/bin/bash
#SBATCH --job-name=sep_chr

ml plink/1.90

cd /fh/fast/sun_w/licai/eQTL_KIRC/data/
mkdir -p chr

for chr in $(seq 1 22); do
plink --file KIRC \
--chr $chr \
--recode \
--out ./chr/KIRC.chr$chr ;
done

