
library(data.table)
library(readxl)
library(stringr)

# --------------------------------------------------------------------
# read in resutls
# --------------------------------------------------------------------

tp53_file = "long_quasi_Whole_Blood_tp53.csv"
tp53 = fread(file.path("../results/GTEx8_Whole_Blood_summary/", tp53_file))
dim(tp53)
tp53[1:2,]

summary(tp53$p_cond)
pi0 = 2*mean(tp53$p_cond > 0.5)
pi0

qval = nrow(tp53)*pi0*tp53$p_cond/rank(tp53$p_cond)
table(tp53$p_cond < 0.01)
table(qval < 0.01)
table(qval < 0.05)
table(qval < 0.10)
table(qval < 0.15)
table(qval < 0.20)

tp53$qval = qval
tp53$ens_id = str_extract(tp53$id, '(\\S+)(?=\\.\\d+)')

dim(tp53)
tp53[1:2,]

# --------------------------------------------------------------------
# read in TP53 targets
# --------------------------------------------------------------------

tp53_anno = read_excel("../Reference/TP53/onc2016502x4.xlsx")
dim(tp53_anno)
tp53_anno[1:2,]

table(tp53_anno$`ensembl ID` %in% tp53$ens_id)
tp53_anno[!tp53_anno$`ensembl ID` %in% tp53$ens_id,]

tp53_target = tp53$ens_id %in% tp53_anno$`ensembl ID`
table(tp53_target)

table(qval < 0.05, tp53_target)
table(qval < 0.10, tp53_target)
table(qval < 0.20, tp53_target)

chisq.test(qval < 0.05, tp53_target)
chisq.test(qval < 0.10, tp53_target)
chisq.test(qval < 0.20, tp53_target)

gc()
sessionInfo()
q(save = "no")
