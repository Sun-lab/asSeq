# Checking the pipeline input genotype_calls.txt vs. the imputed output

sample_id = "TCGA-B0-5099-01A"
# 1.1 Load the SNP 6 csv file
anno = read.table('/fh/fast/sun_w/licai/_tumor_eQTL/GenomeWideSNP_6-na35-annot-csv/GenomeWideSNP_6.na35.annot.csv',
                  sep = ',', header = T, as.is = T)
snpcsv = anno[,c(1:4, 9:10)]
rm(anno)
head(snpcsv)

# 2. Query using COAD_genotype_calls.txt
genotype_calls = read.table('/fh/fast/sun_w/licai/eQTL_KIRC/data/KIRC_genotype_calls.txt', 
                            sep = '\t',header = T, as.is = T, check.names = F)
dim(genotype_calls)
genotype_calls[1:5,1:5]

# 2.1. Get the column corresponding to TCGA-G4-6320-10A-01D-1718-01
# dict = read.table(file = '/fh/fast/sun_w/licai/_tumor_eQTL/data/COAD_clinical_meta.txt',
#                   sep = '\t', as.is = T, header = T)
# sample_cel_file = dict$file_name[dict$associated_entities.entity_submitter_id == sample_id]
# sample_cel_file = gsub('.birdseed.data.txt','', sample_cel_file)
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
# output="/lustre/scr/c/h/chongjin/COAD/stepA_phasing_and_imputation/output/step9_check_het_probes_genotype_calls_vs_shapeit"
# system(sprintf("mkdir %s", output))
# setwd(output)
# write.table(genotype_before_phasing, file="genotype_calls.txt", row.names=FALSE, quote=FALSE, sep="\t")

# 3. Load the combined SNPs list
snps_dir  = dir(path="/fh/fast/sun_w/licai/eQTL_KIRC/data_snp",
                pattern=sample_cel_file,
                full.names=TRUE)
snps_dir  # Should contain ind380_TCGA-G4-6320-10A-01D-1718-01, else something is wrong
setwd(snps_dir)
system(sprintf("cd %s ; rm combined.txt ; for i in `seq 22` ; do cut -d' ' -f1-4 chr${i}.txt | sed -n 's/^chr//p' >> combined.txt ; done", snps_dir)) # rbind snp files for all chr
snps = read.table(paste0(snps_dir, "/combined.txt"), header=FALSE)
colnames(snps) = c("Chromosome", "Physical.Position", "Allele.A", "Allele.B")
snps$Physical.Position = as.numeric(as.character(snps$Physical.Position))
nrow(snps)
genotype_imputed = snps
genotype_imputed$Chromosome = factor(genotype_imputed$Chromosome, levels=c(as.character(1:22), "X", "Y"))
genotype_imputed = genotype_imputed[order(genotype_imputed$Chromosome, genotype_imputed$Physical.Position), ]
head(genotype_imputed)

# 4. Merge the two lists
genotype_merged = merge(x=genotype_before_phasing, y=genotype_imputed, by=c("Chromosome", "Physical.Position"))
nrow(genotype_merged)
head(genotype_merged)
table(genotype_merged[,sample_cel_file])
table(genotype_merged$Allele.A.y == genotype_merged$Allele.B.y, genotype_merged[,sample_cel_file])

#write.table(genotype_merged, file="/fh/fast/sun_w/licai/_tumor_eQTL/R_batch1/stepB_phasing_and_imputation/pipeline/genotype_calls_vs_imputed_snps.txt", row.names=FALSE, quote=FALSE, sep="\t")


sessionInfo()
