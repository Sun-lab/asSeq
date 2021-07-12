
## when writing data into text file, it may use scientific format
## when you read it into c, and using atoi. it will make mistakes
## say 97000000 is written as 9.7e+07, and c think it is 9
## options("scipen") can control write out behavior

options(scipen=20)

# --------------------------------------------------------- 
# organize the annotation in exon level
# ---------------------------------------------------------

setwd("~/research/data/human/")

# --------------------------------------------------------- 
# ensemble genes
# ---------------------------------------------------------

ff  = "Homo_sapiens.GRCh37.66.exon.gtf"
date()
inf = read.table(ff, sep="\t", as.is=TRUE, quote="")
date()

dim(inf)
inf[1:2,]

names(inf) = c("chr", "source", "feature", "start", "end", 
"score", "strand", "frame", "anno")

sort(table(inf$source))
table(inf$chr)
table(inf$feature)
table(inf$score)
table(inf$strand)
table(inf$frame)

# --------------------------------------------------------- 
# whether any exon share the same start and end locations
# ---------------------------------------------------------

id  = paste(inf$chr, inf$start, inf$end, sep=":")
uid = unique(id)

length(id)
length(uid)

tid = table(id)
table(tid)
sort(tid, decreasing=TRUE)[1:3]

# --------------------------------------------------------- 
# sort the exons
# ---------------------------------------------------------

od     = order(inf$chr, inf$start, inf$end)
inf    = inf[od,]

# --------------------------------------------------------- 
# obtain gene_id
# --------------------------------------------------------- 

reg1   = regexpr('gene_id\\s"(\\S+)";', inf$anno, perl=TRUE)
len1   = attributes(reg1)[[1]]
nadd   = length(unlist(strsplit("gene_id", split=""))) + 2
geneId = substr(inf$anno, reg1+nadd, reg1+len1-3)

# --------------------------------------------------------- 
# obtain transcript_id
# --------------------------------------------------------- 

reg1   = regexpr('transcript_id\\s"(\\S+)";', inf$anno, perl=TRUE)
len1   = attributes(reg1)[[1]]
nadd   = length(unlist(strsplit("transcript_id", split=""))) + 2
tranId = substr(inf$anno, reg1+nadd, reg1+len1-3)

# --------------------------------------------------------- 
# obtain gene_name
# --------------------------------------------------------- 

reg1   = regexpr('gene_name\\s"(\\S+)";', inf$anno, perl=TRUE)
len1   = attributes(reg1)[[1]]
nadd   = length(unlist(strsplit("gene_name", split=""))) + 2
geneNm = substr(inf$anno, reg1+nadd, reg1+len1-3)

# --------------------------------------------------------- 
# obtain transcript_name
# --------------------------------------------------------- 

reg1   = regexpr('transcript_name\\s"(\\S+)";', inf$anno, perl=TRUE)
len1   = attributes(reg1)[[1]]
nadd   = length(unlist(strsplit("transcript_name", split=""))) + 2
tranNm = substr(inf$anno, reg1+nadd, reg1+len1-3)

# ---------------------------------------------------------
# find those duplicated exons
# ---------------------------------------------------------

nn = nrow(inf)
nn

wSame = (inf$chr[-nn] == inf$chr[-1])
wSame = wSame & (inf$start[-nn] == inf$start[-1])
wSame = wSame & (inf$end[-nn] == inf$end[-1])

## the last entry will never been dropped
wSame = c(wSame, FALSE)
length(wSame)

table(wSame)

# ---------------------------------------------------------
# collapse annoations for the deleted row
# ---------------------------------------------------------

# suppose there are 6 records altogether, and 
# the 4th and 5th records are the same
# wSame = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)
# so I want to combine the 4th and 5th records
# rev(as.numeric(!wSame)) = (1 1 0 1 1 1)
# cumsum(rev(as.numeric(!wSame))) = (1 2 2 3 4 5)
# rev(cumsum(rev(as.numeric(!wSame)))) = (5 4 3 2 2 1)

id2kp      = rev(cumsum(rev(as.numeric(!wSame))))
pasteUniqu = function(v){paste(unique(v),collapse=":")}
geneId2use = tapply(geneId, id2kp, pasteUniqu)

if(!is.character(geneId2use)){
  stop("I expecte geneID2use to be a character vector :(\n")
}

message("there are ", length(geneId2use), " unique exons.")

xx = grep(":", geneId2use)
if(length(xx) > 0){
  message(length(xx), " exons belong to more than one gene.")
  geClusters = strsplit(geneId2use[xx], split=":")
  t1         = table(sapply(geClusters, length))
  message("their distributuion is")
  print(t1)
}

tranId2use = tapply(tranId, id2kp, pasteUniqu)
xx = grep(":", tranId2use)
if(length(xx) > 0){
  message(length(xx), " exons belong to more than one transcript.")
  trClusters = strsplit(tranId2use[xx], split=":")
  t1         = table(sapply(trClusters, length))
  message("their distributuion is")
  print(t1)
}


geneNm2use = tapply(geneNm, id2kp, pasteUniqu)
tranNm2use = tapply(tranNm, id2kp, pasteUniqu)

# ---------------------------------------------------------
# drop duplicated exons
# ---------------------------------------------------------

infNew = list()
nms = names(inf)
nms[1:8]

for(i in 1:8){
  nm1 = nms[i]  
  cat(i, nm1, "\n")

  if(nm1 == "source" || nm1 == "strand"){
    it1 = tapply(inf[[nm1]], id2kp, pasteUniqu)
  }else{
    it1 = tapply(inf[[nm1]], id2kp, unique)
  }
  
  if(mode(it1) == "list") { stop("hm... non unique ", nm1, "\n") }
  
  infNew[[nm1]] = it1
}

geneId2use = paste("gene_id \"", geneId2use, "\";", sep="")
tranId2use = paste("transcript_id \"", tranId2use, "\";", sep="")
geneNm2use = paste("gene_name \"", geneNm2use, "\";", sep="")
tranNm2use = paste("transcript_name \"", tranNm2use, "\";", sep="")

infNew$anno = paste(geneId2use, tranId2use, geneNm2use, tranNm2use, sep=" ")

infNew = as.data.frame(infNew)
dim(infNew)
infNew[1:2,]

# --------------------------------------------------------- 
# sort the exons
# ---------------------------------------------------------

od     = order(infNew$chr, infNew$start, infNew$end)
any(diff(od) < 0)
infNew = infNew[od,]

# ---------------------------------------------------------
# write out
# ---------------------------------------------------------

id  = paste(infNew$chr, infNew$start, infNew$end, sep=":")
uid = unique(id)

length(id)
length(uid)

write.table(infNew, file = "Homo_sapiens.GRCh37.66.unique.exon.gtf", 
  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, 
  col.names = FALSE)

