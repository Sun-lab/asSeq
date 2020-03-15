#!/bin/bash
#SBATCH --job-name=sep_chr

cd ../../data_genotype_all
mkdir -p chr
ml plink/1.90

# for chr in $(seq 1 22);
# do
# plink --file GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC_RMDUP_biallelic \
# --chr $chr \
# --recode \
# --out ./chr/geno.chr$chr ;
# done


sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 1 --recode --out ./chr/geno.chr1"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 2 --recode --out ./chr/geno.chr2"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 3 --recode --out ./chr/geno.chr3"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 4 --recode --out ./chr/geno.chr4"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 5 --recode --out ./chr/geno.chr5"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 6 --recode --out ./chr/geno.chr6"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 7 --recode --out ./chr/geno.chr7"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 8 --recode --out ./chr/geno.chr8"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 9 --recode --out ./chr/geno.chr9"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 10 --recode --out ./chr/geno.chr10"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 11 --recode --out ./chr/geno.chr11"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 12 --recode --out ./chr/geno.chr12"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 13 --recode --out ./chr/geno.chr13"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 14 --recode --out ./chr/geno.chr14"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 15 --recode --out ./chr/geno.chr15"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 16 --recode --out ./chr/geno.chr16"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 17 --recode --out ./chr/geno.chr17"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 18 --recode --out ./chr/geno.chr18"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 19 --recode --out ./chr/geno.chr19"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 20 --recode --out ./chr/geno.chr20"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 21 --recode --out ./chr/geno.chr21"
sbatch --wrap="plink --file GTEx_Analysis_max2allel_rmDup --chr 22 --recode --out ./chr/geno.chr22"
