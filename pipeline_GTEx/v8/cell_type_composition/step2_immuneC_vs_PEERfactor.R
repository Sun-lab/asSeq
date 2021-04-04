rm(list = ls())
setwd("/fh/fast/sun_w/licai/_tumor_eQTL/GitHub/asSeq/pipeline_GTEx/v8/cell_type_composition/")

# -----------------------------------------------------------------
# Read in immune cell composition data
# -----------------------------------------------------------------
library(reshape2)
library(ggplot2)

imc = read.table("data/CIBERSORTx_Results.txt", header = T, sep="\t", 
                 as.is = T, check.names = F)
imc[1:5,1:5]
dim(imc)
colnames(imc)

rownames(imc) = imc$Mixture
imc = imc[, -c(1, 24:26)]
dim(imc)

# ----------------------------------------------------------------------
# correlation test between Immune cell compositon and peer factor
# ----------------------------------------------------------------------

Covariates = read.table("data/Whole_Blood.v8.covariates.txt", header = T, 
                        sep ='\t', as.is = T, check.names = F)
dim(Covariates)
Covariates[1:3,1:2]

rownames(Covariates) = Covariates$ID
Covariates = Covariates[, -1]

#colnames(Covariates)
peerf = grep('InferredCov', rownames(Covariates))
Covariates = Covariates[peerf, ]
rownames(Covariates)

# match samples 
table(rownames(imc)  == colnames(Covariates))

refic = which.max(apply(imc, 2, median))
refic

ref = imc[, refic]
log_imc = log((imc[,-refic] + 10^-6)/(ref + 10^-6)) 
log_imc[1:5,]
log_imc = data.frame(log_imc)

fit1 = lm(t(Covariates)[,1] ~ ., data = log_imc)
summary(fit1)

df1 = data.frame(fitted=fitted(fit1), peer1=t(Covariates)[,1])
dim(df1)
df1[1:2,]

saveRDS(df1, file="data/peer1_fitted_values.rds")


df1 = readRDS("data/peer1_fitted_values.rds")
dim(df1)
df1[1:2,]

pdf("./figures/PEERfactor1.vs.Fitted_whole_blood.pdf", width = 3, height = 3)
par(mar=c(5,4,1,1), bty="n")
plot(df1$fitted, df1$peer1, xlab = 'Fitted values',
     ylab = "PEER factor 1", main = '', pch = 19, cex=0.5)
dev.off()


pdf("./figures/PEERfactor1.vs.Neutrophils.pdf", width = 3, height = 3)
par(mar=c(5,4,1,1), bty="n")
plot(ref, df1$peer1, xlab = 'Neutrophils',
     ylab = "PEER factor 1", main = '', pch = 19, cex=0.5)
dev.off()

q("no")
