
annoVersion = "Mus_musculus.NCBIM37.67"

setwd("~/research/data/mouse/")

# --------------------------------------------------------- 
# read in ensemble gene annotation
# ---------------------------------------------------------

ref  = sprintf("%s.exon.gtf", annoVersion)
info = read.table(ref, sep="\t", as.is=TRUE)
dim(info)
info[1:2,]

names(info) = c("chr", "source", "feature", "start", "end", 
"score", "strand", "frame", "anno")

table(info$chr)
table(info$source)
table(info$feature)

info   = info[info$feature == "exon",]
dim(info)

# --------------------------------------------------------- 
# obtain gene_id
# --------------------------------------------------------- 

reg1   = regexpr('gene_id\\s(\\S+);', info$anno, perl=TRUE)
len1   = attributes(reg1)[[1]]
nadd   = length(unlist(strsplit("gene_id", split=""))) + 1
info$geneId = substr(info$anno, reg1+nadd, reg1+len1-2)

# --------------------------------------------------------- 
# obtain transcript_id
# --------------------------------------------------------- 

reg1   = regexpr('transcript_id\\s(\\S+);', info$anno, perl=TRUE)
len1   = attributes(reg1)[[1]]
nadd   = length(unlist(strsplit("transcript_id", split=""))) + 1
info$tranId = substr(info$anno, reg1+nadd, reg1+len1-2)

# --------------------------------------------------------- 
# obtain transcript length
# --------------------------------------------------------- 

infoTLen = tapply(info$end - info$start+1, info$tranId, sum)

# --------------------------------------------------------- 
# calculate the number of exons per gene and 
# the nubmer of isoforms per gene
# --------------------------------------------------------- 

nT = tapply(info$tranId, info$geneId, function(v){length(unique(v))} )
nE = tapply(info$tranId, info$geneId, length )

all(names(nT) == names(nE))
nTE = data.frame(geneId=names(nE), nE = nE, nT=nT)
dim(nTE)

# --------------------------------------------------------- 
# read in information of unique exons 
# --------------------------------------------------------- 

bedF = sprintf("~/research/data/mouse/%s.nonoverlap.exon.bed", annoVersion)

bedD = read.table(bedF, sep="\t", as.is=TRUE)
dim(bedD)
bedD[1:2,]

names(bedD) = c("chr", "start", "end", "name", "score", "strand")

ids  = matrix(unlist(strsplit(bedD$name, split="|", fixed=TRUE)), byrow=TRUE, ncol=3)
dim(ids)
ids  = unique(ids[,1:2], MARGIN=1)
dim(ids)

colnames(ids) = c("clustID", "geneID")
ids  = as.data.frame(ids, stringsAsFactors=FALSE)
ids[1:2,]

# --------------------------------------------------------- 
# map between transcript ID and gene ID 
# --------------------------------------------------------- 

gID2 = strsplit(ids$geneID, split=":")
cID2 = rep(ids$clustID, times=sapply(gID2, length))
ids2 = data.frame(clustID=cID2, geneID=unlist(gID2), stringsAsFactors=FALSE)

dim(ids2)
ids2 = unique(ids2, MARGIN=1)
dim(ids2)

tb1  = table(ids2$clustID)
table(tb1)

length(unique(ids$clustID))
length(unique(ids2$clustID))

mat1 = match(nTE$geneId, ids2$geneID)

if(!(all(ids2$geneID[mat1] == nTE$geneId))){
  stop("geneID not mapped\n")  
}

nTE$clustID =ids2$clustID[mat1]
rownames(nTE) = NULL

nTE$geneId  = as.character(nTE$geneId)
nTE$clustID = as.character(nTE$clustID)

save(nTE, file=sprintf("%s.nTE.RData", annoVersion))
