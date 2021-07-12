library(isoform)

setwd("~/research/data/human/")

# ---------------------------------------------------------
# load information of known isoforms
# ---------------------------------------------------------

nonOverlapExonFile  = "Homo_sapiens.GRCh37.66.nonoverlap.exon.gtf"
isoAll = knownIsoforms(nonOverlapExonFile)

length(isoAll)
nIso = sapply(isoAll, ncol)
summary(nIso)
table(nIso)

save(isoAll, file="Homo_sapiens.GRCh37.66.nonoverlap.exon.knownIsoforms.RData")
