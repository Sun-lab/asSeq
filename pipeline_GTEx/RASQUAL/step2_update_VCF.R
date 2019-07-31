args = commandArgs(TRUE)
#args = c("22")
chri = args[1]

#module add tabix;module add gcc;module add gsl

data.dir = "/fh/scratch/delete90/sun_w/_GTEx/ncbi"
bam.dir = sprintf("%s/wb_raw_bamfiles", data.dir)
geno.dir = "../data_genotype_all" 
cnt.dir = sprintf("%s/cnt", geno.dir)
vcf.dir = sprintf("%s/vcf", geno.dir)
vcfr.dir = sprintf("%s/vcfr", geno.dir)
if(!file.exists(vcfr.dir))dir.create(vcfr.dir)
info.dir = ".."
dir.rasc = sprintf("%s/gatkasc", data.dir)

snp = sprintf("%s/GTEx_phased_%s.vcf.gz", vcf.dir, chri)
if(!file.exists(sprintf("chrisnps_%s.txt", chri)))
  system(sprintf("gunzip -c %s | awk '{print $1,$2}' > chrisnps_%s.txt", snp,chri))

step5    =  "../R_batch4_whole_blood/stepC/output/step5_TReC_per_gene_filterIt/"
step5_file = sprintf("gene_level_counts_filter_out_low_expressed_chr%s.txt",chri)

genes = read.table(sprintf("%s/gencode.v24lift37.annotation_nchr.gtf",
                           info.dir), as.is=T, skip = 5, sep = "\t")
genes = genes[which(genes[,1] == chri), ]
geneid = sapply(strsplit(genes[,9], ";"), "[", 1)
geneid = gsub("gene_id ", "", geneid)
geneid[1:10]
genes[,10] = geneid
# geneid = unique(geneid)
# to.rm = grep(":", geneid)
# if(length(to.rm)>0)geneid=geneid[-to.rm]
geneid = system(sprintf("cut -f1 %s%s", step5, step5_file), intern = T)[-1]
length(geneid)
sam2kp = read.table(sprintf("%s/whole_blood_samples_eQTL.txt", geno.dir), 
                    as.is = T)
sam2kp = as.character(sam2kp$V1)
length(sam2kp)

#
#a. did in permute_VCF
#count number of SNPs within a window - will be used as input for RASQUAL
#  
# for(chri in 1:23){
#   system(sprintf("grep '^%s ' tmpsnps.txt > chrisnps_%s.txt", chri,chri))
# }
# win = 1e5
# snp.cnt = sprintf("%s/number_of_snps_per_gene_win%s.csv", info.dir, win)
# if(!file.exists(snp.cnt)){
#   res = data.frame(matrix(NA, nrow=length(geneid), ncol=3))
#   for(i in 1:length(geneid)){
#     #m = match(geneid[i], genes[,11])
#     gnm = which(geneid[i]==genes[,10])
#     chri = genes[gnm[1],1]
#     chrsnps = read.table(sprintf("chrisnps_%s.txt", chri))
#     posst = min(genes[gnm,4])
#     posen = max(genes[gnm,5])
#     kp = which(((posst-win)<=chrsnps[,2]) & ((posen+win)>chrsnps[,2]))
#     res[i, 1] = chri
#     res[i, 2] = geneid[i]
#     res[i, 3] = length(kp)
#     message(i)
#   }
#   write.csv(res, snp.cnt, quote=F, row.names=F) 
# }

#
#b.
#add snp level allele-specific counts to vcf
#
vcf.hd = read.table(snp, as.is=T, nrow=4, comment.char=",", sep=",")
vcf.hd[3,1] = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
vcf.in = read.table(snp, as.is=T)

#AF in this file came from a bigger file, 
#replace it with frequency appropriate for these samples
count.ref = function(vec){mean(unlist(strsplit(vec, split="\\|"))==1)}
af = apply(vcf.in[,-(1:9)], 1, count.ref)
summary(af)
vcf.in[,8] = af
vcf.in[,9] = "GT:AS"
#if too many reads come from wrong bases - may skip this count
frac = 0.1
samples = unlist(strsplit(vcf.hd[4,1], split="\t"))[-(1:9)]
samples = gsub("_", "", substr(samples, 1, 10))
sam2rm = NULL
for(i in 1:length(samples)){
  alt = ref = numeric(nrow(vcf.in))
  sampli = samples[i]
  if(sampli %in% sam2kp ){
    # system(sprintf("grep ^%s %s/%s.txt > tmpRASE.txt",
    #                chri,dir.rasc, sampli))
    # rasc = try(read.table("tmpRASE.txt",
                      # header=F, as.is = T))# chri
    rasc = read.table(sprintf("%s/%s.txt",dir.rasc, sampli),
                      header = T, as.is = T)
    rasc = rasc[which(rasc$contig == chri), ]
    # if(class(rasc) == "try-error"){
    #   rasc = NULL
    # }
    m = match(rasc[,3], vcf.in[,3])
    if(any(is.na(m))){warning(sum(is.na(m)), " unexpected snps")}
    ref[na.omit(m)] = rasc[,6]
    alt[na.omit(m)] = rasc[,7]
    asc = ref[na.omit(m)]+alt[na.omit(m)]
    rma = rasc[,12]/asc>frac
    ref[na.omit(m)][rma] = alt[na.omit(m)][rma] = 0
    vcf.in[,i+9] = sprintf("%s:%s,%s", vcf.in[,i+9], ref, alt)
    message("individual ", i, " has been processed")
  }else{
    sam2rm = c(sam2rm, i) # rm samples that is not whole boold 
  }
}
length(sam2rm)
if(length(sam2rm) > 0){
  line4 = c(unlist(strsplit(vcf.hd[4,1], split="\t"))[(1:9)], 
            samples[-sam2rm])
  line4 = paste0(line4, collapse = "\t")
  vcf.hd[4,1] = line4
  vcf.in = vcf.in[, -(sam2rm+9)]
}

#save to rasqual count folder
out.vcf = sprintf("%s/SNP_chr%s.vcf", vcfr.dir, chri)
write.table(vcf.hd, out.vcf, quote=F, row.names=F, col.names=F, sep="\t")
write.table(vcf.in, out.vcf, quote=F, row.names=F, col.names=F, sep="\t", append=T)
system(sprintf("bgzip -f %s", out.vcf))
out.vcf = sprintf("%s.gz", out.vcf)
system(sprintf("tabix %s", out.vcf))

q("no")


