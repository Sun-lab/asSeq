
annoVersion = "Mus_musculus.NCBIM37.67"

library(isoform)

setwd("~/research/data/mouse/")

# ---------------------------------------------------------
# load information of known isoforms
# ---------------------------------------------------------

nonOverlapExonFile  = sprintf("%s.nonoverlap.exon.gtf", annoVersion)
isoAll = knownIsoforms(nonOverlapExonFile)

length(isoAll)
nIso = sapply(isoAll, ncol)
summary(nIso)
table(nIso)

outFile = sprintf("%s.nonoverlap.exon.knownIsoforms.RData", annoVersion)

save(isoAll, file=outFile)
