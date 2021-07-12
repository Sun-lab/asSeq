
annoVersion = "Mus_musculus.NCBIM37.67"

## when writing data into text file, it may use scientific format
## when you read it into c, and using atoi. it will make mistakes
## say 97000000 is written as 9.7e+07, and c think it is 9
## options("scipen") can control write out behavior

options(scipen=20)

# --------------------------------------------------------- 
# read in data
# ---------------------------------------------------------

setwd("~/research/data/mouse/")

ff  = sprintf("%s.nonoverlap.exon.gtf", annoVersion)

date()
inf = read.table(ff, sep="\t", as.is=TRUE, header=FALSE, quote="")
date()

names(inf) = c("chr", "source", "feature", "start", "end", 
"score", "strand", "frame", "anno")

dim(inf)
inf[1:2,]

table(inf$chr)
table(inf$strand)

summary(inf$end - inf$start)

table(inf$end - inf$start == 0)

## for bed format, the first base in a chromosome is numbered 0.
## while in gtf format, the first base in a chromosome is numbered 1.

inf$start = inf$start - 1

# --------------------------------------------------------- 
# obtain clust_id
# --------------------------------------------------------- 

reg1   = regexpr('clustId\\s"(\\S+)";', inf$anno, perl=TRUE)
len1   = attributes(reg1)[[1]]
nadd   = length(unlist(strsplit("clustId", split=""))) + 2
clustId = substr(inf$anno, reg1+nadd, reg1+len1-3)

# --------------------------------------------------------- 
# obtain gene_id
# --------------------------------------------------------- 

reg1   = regexpr('gene_id\\s"(\\S+)";', inf$anno, perl=TRUE)
len1   = attributes(reg1)[[1]]
nadd   = length(unlist(strsplit("gene_id", split=""))) + 2
geneId = substr(inf$anno, reg1+nadd, reg1+len1-3)

# --------------------------------------------------------- 
# obtain exon_id
# --------------------------------------------------------- 

reg1   = regexpr('exon_id\\s"(\\S+)";', inf$anno, perl=TRUE)
len1   = attributes(reg1)[[1]]
nadd   = length(unlist(strsplit("exon_id", split=""))) + 2
exonId = substr(inf$anno, reg1+nadd, reg1+len1-3)

# --------------------------------------------------------- 
# construct bed file
# --------------------------------------------------------- 

names = paste(clustId, geneId, exonId, sep="|")
score = rep("666", length(names))
bed   = cbind(inf$chr, inf$start, inf$end, names, score, inf$strand)

# ---------------------------------------------------------
# write out results
# ---------------------------------------------------------

setwd("~/research/data/mouse/")

outFile = sprintf("%s.nonoverlap.exon.bed", annoVersion)

write.table(bed, col.names = FALSE, append = FALSE, 
file = outFile, quote = FALSE, sep = "\t", row.names = FALSE)

