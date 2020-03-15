args = commandArgs(TRUE)
#args = c("22")
chri = as.numeric(args[1])


get_block = function(str,split="_",block=1){
  unlist(strsplit(str,split=split))[block]
}

options("scipen"=999,"digits"=4)

root.dir = "/fh/fast/sun_w/licai/_GTEx/"
data.dir = sprintf("%s/data_genotype_all/", root.dir)
# info.dir = sprintf("%s/info", root.dir)
pipe.dir = sprintf("%s/R_batch5_wb_RASQUAL", root.dir)
shape.dir = sprintf("%s/shapeit", root.dir)

geno.dir = sprintf("%s/data_genotype_all/", root.dir)
phas.dir = sprintf("%s/phased",geno.dir)
# refinfo = sprintf("%s/1000GP_Phase3",info.dir)
# impu.dir = sprintf("%s/imputed",geno.dir)

# out.dir = sprintf("%s/snps_omni",geno.dir)
# if(!file.exists(out.dir))dir.create(out.dir)
# lift.dir = sprintf("%s/liftover",root.dir)
vcf.dir = sprintf("%s/vcf",data.dir)
if(!file.exists(vcf.dir))dir.create(vcf.dir)

fls = list.files(phas.dir)
hap = fls[grep("hap",fls)]

alid = read.table(sprintf("%s/phasing_chr1.sample",phas.dir),
                  as.is=T,skip=2)[,2]

# will filter out sample in step2
# samf = "../R_batch4_whole_blood/stepC/output/step5_TReC_per_gene_filterIt/"
# smpls = read.table(sprintf("%s/read_depth.txt",samf), 
#                    as.is=T,header=T)
# smpls = rownames(smpls)
# smpli = match(smpls,alid)
# all(smpli==sort(smpli))

bases=c("A","C","G","T")
# maf = 0.01 will keep all sample

hps = fls[grep(".haps",fls)]
hps
off = 5
hpj = read.table(sprintf("%s/phasing_chr%s.haps",phas.dir,chri),as.is=T)
hpj$V4 = as.character(hpj$V4)
hpj$V5 = as.character(hpj$V5)

rm.indels =  which(hpj$V4%in%bases & hpj$V5%in%bases)
hpj = hpj[rm.indels,]
dim(hpj)

inf = hpj[, 1:5]   
dat = hpj[,-(1:5)]

rm(hpj)

lif2 = cbind(inf[,c(1,3,2,4,5)], ".", "PASS", 0, "GT")
anyDuplicated(lif2[,3])
table(lif2[,2] == sort(lif2[,2]))
smpli = 1:length(alid)
indloc = (smpli-1)*2+1
indloc2 = (smpli-1)*2+2   
vcf = matrix(".", nrow = nrow(dat), ncol = length(alid))
for(indi in smpli){
  vcf[, indi] = paste(dat[,indloc[indi]], dat[,indloc2[indi]], sep = "|")
}
out = cbind(lif2, vcf)
str = matrix(c("##fileformat=VCFv4.0",
               "##fileDate=20150416",
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"))
str[4] = paste(str[4],paste(alid,collapse = "\t"),sep = "\t")
write.table(str,sprintf("%s/GTEx_phased_%s.vcf",vcf.dir,chri),
            quote=F,row.names=F,col.names=F, append = F)

write.table(out,sprintf("%s/GTEx_phased_%s.vcf",vcf.dir,chri),
            quote=F,row.names=F,col.names=F, append = T, sep="\t")

rm(out)
rm(vcf)
gc()
message(chri)

q("no")


for(i in (22:1)){
  com<-sprintf("sbatch --exclusive R CMD BATCH --no-save --no-restore '--args %s' step0_prepare_vcf.R step0_prepare_vcf_%s.Rout",i,i)
  message(com)
  # system(com)
}


