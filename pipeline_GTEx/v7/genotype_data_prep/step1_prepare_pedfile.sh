#1. copy vcf file to working (not necessary)
cp /fh/fast/sun_w/GTEx/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz .

# only keep allel = 2   https://www.biostars.org/p/141156/
sbatch --wrap="bcftools view -m2 -M2 -v snps  GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz > GTEx_Analysis_max2allel.vcf"

# rm duplicate position
sbatch --wrap="plink --vcf GTEx_Analysis_max2allel.vcf --out GTEx_Analysis_max2allel --recode"

## duplicated positions need to be removed -f2 ignored the first two field| -D list all the duplcated one | filter out missing rate > 10%
uniq -f2 -D GTEx_Analysis_max2allel.map | cut -f2 > dupeSNP.txt
plink --file GTEx_Analysis_max2allel -exclude dupeSNP.txt --geno 0.1 --recode --out GTEx_Analysis_max2allel_rmDup
