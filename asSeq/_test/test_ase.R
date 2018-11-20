
library(asSeq)
source("~/research/R/asSeq/_test/aseR.R")

data(KLK1)
attach(KLK1)

ase(Y1, Y2, Z, output.tag="tes_ase", p.cut=0.5, local.only=FALSE)

rest = read.table("tes_ase_eqtl.txt", header=TRUE, sep="\t")
rest

system("rm tes_ase_eqtl.txt")
system("rm tes_ase_freq.txt")

zeta = as.numeric(abs(Z-2)<1.5)
yT = Y1+Y2
y2 = Y2
y2[Z==3] = Y1[Z==3]

a1 = aseR(y2, yT, zeta)
p1 = 1 - pchisq(2.0*a1$logLikH1 - 2.0*a1$logLikH0, 1)
