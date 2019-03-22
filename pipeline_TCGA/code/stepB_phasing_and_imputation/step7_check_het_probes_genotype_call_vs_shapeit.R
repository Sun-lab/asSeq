# Checking the pipeline input genotype calls vs. the phased output

sample_id = "sample1"  # usually the format is "TCGA-B0-5099-01A"
genotype_file = '../../data/sample_geno_call_22.txt' # change to the name of genotype files 
anno_file = '/fh/fast/sun_w/licai/_tumor_eQTL/GenomeWideSNP_6-na35-annot-csv/GenomeWideSNP_6.na35.annot.csv'

# 1.1 Load the SNP 6 csv file
anno = read.table(anno_file, sep = ',', header = T, as.is = T)
snpcsv = anno[,c(1:4, 9:10)]
rm(anno)
head(snpcsv)

# 2. Query using genotype_calls.txt
genotype_calls = read.table(genotype_file, 
                            sep = '\t',header = T, as.is = T, check.names = F)
dim(genotype_calls)
genotype_calls[1:5,1:5]

# 2.1. Get the column corresponding to sample_id
sample_cel_file = substr(sample_id, 1, 12) 
sample_cel_file 
sample_birdseed_call = genotype_calls[,c(which(sample_cel_file == colnames(genotype_calls))), drop = F]
head(sample_birdseed_call)

# 2.2 Merge birdseed calls with SNP 6 reference
snpcsv$Physical.Position = as.numeric(snpcsv$Physical.Position)
genotype_before_phasing = merge(x=snpcsv, y=sample_birdseed_call,
                                by.x = "Probe.Set.ID", by.y = "row.names", sort=FALSE)
nrow(genotype_before_phasing)
genotype_before_phasing = genotype_before_phasing[genotype_before_phasing$Chromosome!="---",]
nrow(genotype_before_phasing)
genotype_before_phasing$Chromosome = factor(genotype_before_phasing$Chromosome, levels=c(as.character(1:22), "X", "Y"))
head(genotype_before_phasing)
nrow(genotype_before_phasing)

# 2.3 Sort
genotype_before_phasing = genotype_before_phasing[order(genotype_before_phasing$Chromosome, genotype_before_phasing$Physical.Position), ]
head(genotype_before_phasing)

# 3. Load the genotype after phasing using shapeit
phased = "../../data/phased"
sample_table = read.table(sprintf("%s/phasing_chr22.sample", phased), as.is=TRUE)
sample_index = which(sample_table[,2] == gsub('.CEL','', sample_cel_file))

#for (i in 1:22) {
library(doParallel)
library(foreach)
#setup parallel backend to use 6 processors
cl = makeCluster(6, outfile="")
registerDoParallel(cl)
getDoParWorkers()
haps_table = foreach(i=1:22, .combine=rbind) %dopar% {
  haps = read.table(sprintf("%s/phasing_chr%s.haps", phased, i), as.is=TRUE)
  haps[,c(1, 2, 3, 4, 5, 2*sample_index, 2*sample_index + 1)]
}
stopCluster(cl)
#}
#haps_table = do.call("rbind", haps_table_list)
colnames(haps_table) = c("Chromosome", "dbSNP.RS.ID", "Physical.Position", "Allele.A", "Allele.B", "hap1", "hap2")
head(haps_table)
# write.table(haps_table, file="genotype_shapeit.txt", row.names=FALSE, quote=FALSE, sep="\t")

# 4. Merge the two lists
genotype_merged = merge(x=genotype_before_phasing, y=haps_table, by=c("dbSNP.RS.ID", "Chromosome", "Physical.Position"))
genotype_merged = genotype_merged[order(genotype_merged$Chromosome, genotype_merged$Physical.Position),]
nrow(genotype_merged)
head(genotype_merged)
table(genotype_merged[,sample_cel_file])
table(genotype_merged$Allele.A.x, genotype_merged$Allele.A.y)
table(genotype_merged[,sample_cel_file] , (genotype_merged$hap1 +genotype_merged$hap2))

#write.table(genotype_merged, file="step7_genotype_calls_vs_shapeit.txt", row.names=FALSE, quote=FALSE, sep="\t")

sessionInfo()

q(save = 'no')