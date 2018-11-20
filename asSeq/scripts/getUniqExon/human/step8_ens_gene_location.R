
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

# --------------------------------------------------------- 
# obtain gene_id
# --------------------------------------------------------- 

reg1   = regexpr('gene_id\\s"(\\S+)";', inf$anno, perl=TRUE)
len1   = attributes(reg1)[[1]]
nadd   = length(unlist(strsplit("gene_id", split=""))) + 2
geneId = substr(inf$anno, reg1+nadd, reg1+len1-3)

chrs   = tapply(inf$chr, geneId, unique)
start  = tapply(inf$start, geneId, min)
end    = tapply(inf$end, geneId, max)

all(names(chrs) == names(start))
all(names(chrs) == names(end))

info = data.frame(id=names(chrs), chr=chrs, start=start, end=end)
dim(info)
info[1:2,]

write.table(info, file = "ens_gene_loci.txt", append = FALSE,  
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

