rm(list = ls())
proj = "/fh/fast/sun_w/licai/eQTL_LUAD/"
setwd(proj)
system('mkdir -p ./data')
json_file = "./metadata.cart.2018-09-26.json"

# ----------------------------------------------------------------------
# check samples
# ----------------------------------------------------------------------
library(jsonlite)
library(XML)

# list2df and df2df are two functions to process the JSON file to data.frame
list2df <- function(clinic){
  nms = NULL
  for(i in 1:length(clinic)){
    ci = unlist(clinic[[i]])
    nms = union(nms, names(ci))
  }
  length(nms)
  
  nms = sort(nms)
  
  clinicDF = NULL
  
  for(i in 1:length(clinic)){
    
    ci = unlist(clinic[[i]])
    vi = rep(NA, length(nms))
    vi[match(names(ci), nms)] = ci
    
    clinicDF = rbind(clinicDF, vi)
  }
  
  colnames(clinicDF) = nms
  clinicDF = as.data.frame(clinicDF, row.names=1:nrow(clinicDF),
                           stringsAsFactors=FALSE)
  clinicDF
  
}

df2df <- function(meta){
  
  nms = NULL
  for(i in 1:nrow(meta)){
    ci = unlist(meta[i,])
    nms = union(nms, names(ci))
  }
  length(nms)
  
  meta1 = NULL
  
  for(i in 1:nrow(meta)){
    ci = unlist(meta[i,])
    vi = rep(NA, length(nms))
    ni = intersect(names(ci), nms)
    vi[match(ni, nms)] = ci[match(ni, names(ci))]
    
    meta1 = rbind(meta1, vi)
  }
  
  colnames(meta1) = nms
  meta1 = as.data.frame(meta1, row.names=1:nrow(meta1))
  meta1
  
}

meta = fromJSON(json_file)
meta[1:2,]

meta1 = df2df(meta)
meta1[1:2, ]

barcode = strsplit(as.character(meta1$associated_entities.entity_submitter_id), '-')
barcode = do.call(rbind, barcode)
colnames(barcode) = c('study',"tss", "participant", "sample", "portion", "plate", "center")
barcode[1:2,]  

length(unique(barcode[, 'participant']))
table(barcode[, 'sample'])
table(barcode[, 'plate'])

meta1 = cbind.data.frame(meta1, barcode, stringsAsFactors = F)

# only take normal tissue 
# 10	Blood Derived Normal
# 11	Solid Tissue Normal
w2kp = which(meta1$sample %in% c(paste0(10, LETTERS), paste0(11, LETTERS)))
meta2 = meta1[w2kp, ]
dim(meta2)

tp = table(meta2$participant)
table(tp)

if(any(tp != 1)){
  sam2check = names(tp)[tp > 1]
  ww2check = which(meta2$participant %in% sam2check)
  meta2check = meta2[ww2check,]
  
  sams = unique(meta2check$participant)
  w2kp = NULL
  
  for(sam1 in sams){
    ww1 = which(meta2check$participant == sam1)
    ww2 = ww1[which(meta2check$sample[ww1] %in% paste0(10, LETTERS))]
    if(length(ww2) > 1L){
      ww2 = ww2[which.max(meta2check$archive.file_size[ww2])]
    }
    if(length(ww2) == 0L){
      ww2 = ww1[which.max(meta2check$archive.file_size[ww1])]
    }
    w2kp = c(w2kp, ww2)
  }
  meta3 = rbind(meta2check[w2kp, ], meta2[-ww2check,])
  
}else{
  meta3 = meta2
}
dim(meta3)
tp = table(meta3$participant)
table(tp)
table(meta3$sample)

write.table(meta3, file = './data/LUAD_clinical_meta.txt',
            sep = '\t', quote = F, row.names = F)

barcode = substr(as.character(meta3$associated_entities.entity_submitter_id), 1, 12)
length(unique(barcode))
# ----------------------------------------------------------------------
# Merge all the genotype files
# ----------------------------------------------------------------------
# -1 = NN, 0 = AA, 1 = AB, 2 = BB

setwd('./data_raw')
files = list.files(pattern = 'birdseed.data.txt$',
                   recursive = T)
files[1:5]

samples = sapply(strsplit(files, "/", fixed = T), '[[', 2)
sam2kp = which(samples %in% meta3$file_name)

SNPs_all = system(paste0('cut -f1 ', files[sam2kp]), intern = T)
SNPs_all = SNPs_all[-c(1:2)]

genotype_calls = matrix(NA, nrow =length(SNPs_all), ncol = length(sam2kp))
colnames(genotype_calls) = gsub('.birdseed.data.txt', '', samples[sam2kp])
rownames(genotype_calls) = SNPs_all

for(f in files[sam2kp]){
  sam = strsplit(readLines(f, n =1L), '\t')[[1]][2]
  calls = read.table(f, sep = '\t',header = T, as.is = T, skip = 1)
  if(all(SNPs_all == calls$Composite.Element.REF)){
    genotype_calls[, sam] = calls$Call
  }else{
    genotype_calls[calls$Composite.Element.REF, sam] = calls$Call
  }
  
}

clnm = barcode[match(samples[sam2kp],meta3$file_name)]
colnames(genotype_calls) = clnm

dim(genotype_calls)
genotype_calls[1:5,1:5]

write.table(genotype_calls, file = '../data/genotype_calls.txt',
            sep ='\t', quote = F)

genotype_calls[1:5,1:5]

# ----------------------------------------------------------------------
# create ped/map file - plink format
# ----------------------------------------------------------------------
# there are some missing chr, dbSNP, pos with '---'
# removed them in the ped/map file 
# MAP file 

anno_file = "/fh/fast/sun_w/licai/_tumor_eQTL/GenomeWideSNP_6-na35-annot-csv/"
anno = read.table(paste0(anno_file, 'GenomeWideSNP_6.na35.annot.csv'),
                  sep = ',', header = T, as.is = T)

dim(anno)
table(SNPs_all %in% anno$Probe.Set.ID)

anno$Genetic.Position = 0
map_anno = anno[match(SNPs_all,anno$Probe.Set.ID), c('Chromosome', 'dbSNP.RS.ID','Genetic.Position', 'Physical.Position','Probe.Set.ID', 'Allele.A','Allele.B')]
map_anno$Chromosome = sub("---", NA, map_anno$Chromosome)
map_anno$Physical.Position = sub("---", NA, map_anno$Physical.Position)
map_anno$dbSNP.RS.ID = sub("---", NA, map_anno$dbSNP.RS.ID)
map_anno = map_anno[complete.cases(map_anno),]
dim(map_anno)
map_anno[1:5,]  

map_anno$Physical.Position = as.numeric(map_anno$Physical.Position)
map_anno = map_anno[with(map_anno, order(Chromosome, Physical.Position)),]
map_anno[1:20,]  
table(map_anno$Chromosome)

# remove duplicated cases
dup.cases = duplicated(map_anno[,1:4])
map_anno[dup.cases,]

if(any(dup.cases))
  map_anno = map_anno[!dup.cases,]


write.table(map_anno[,1:4], file ='../data/geno.map', row.names = F, col.names = F, quote = F)

# PED file format 
# Family ID [string]
# Individual ID [string]
# Father ID [string]
# Mother ID [string]
# Sex [integer]
# Phenotype [float]

IndID = colnames(genotype_calls)
FamID = paste0('FamID', seq(ncol(genotype_calls)))
FID = paste0('FID', seq(ncol(genotype_calls)))
MID = paste0('MID', seq(ncol(genotype_calls)))
Sex = rep(-9, ncol(genotype_calls))
Phenotype = rep(-9, ncol(genotype_calls))

genotype_calls = genotype_calls[as.character(map_anno$Probe.Set.ID), ]
table(map_anno$Probe.Set.ID == rownames(genotype_calls))

# convert 0/1/2 to T/C/G/A  
# -1 = NN, 0 = AA, 1 = AB, 2 = BB (Allele.A, Allele.B)
genotype.tcga = matrix(NA, nrow = nrow(genotype_calls), ncol = ncol(genotype_calls))
colnames(genotype.tcga) = colnames(genotype_calls)


for(i in 1:nrow(genotype_calls)){
  ref = as.character(map_anno[i, 'Allele.A'])
  alt = as.character(map_anno[i, 'Allele.B'])
  genotype.tcga[i,] = sapply(genotype_calls[i, ], function(x){if(x == 0)
    return(paste0(ref, ' ' ,ref))
    if(x == 1)
      return(paste0(ref, ' ' ,alt))
    if(x == 2)
      return(paste0(alt, ' ' , alt))
    if(x == -1)
      return('0 0')})
}

geno.ped = cbind(FamID, IndID, FID,MID, Sex, Phenotype, t(genotype.tcga))
geno.ped[1:5,1:10]
table(geno.ped[,'IndID'] == colnames(genotype.tcga))

write.table(geno.ped, file ='../data/geno.ped',
            row.names = F, col.names = F, quote = F, sep ='\t')
q('no')