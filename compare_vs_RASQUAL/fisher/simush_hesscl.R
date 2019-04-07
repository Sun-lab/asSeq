args = commandArgs(trailingOnly = TRUE)
percase = 0.1
args
b0 = as.numeric(args[1])
th = as.numeric(args[2])
dv = as.numeric(args[3])
niter = as.numeric(args[4])
ss = as.numeric(args[5])
nsnp = as.numeric(args[6])
b1 = 0
mn = 100
set.seed(12345)

phi = th/dv
rho = phi/(1+phi)
c(rho, th)

library(numDeriv)
library(rcppreqtl)
ls("package:rcppreqtl")
source("fishercl.R")

Xm = makeXmatr(ss=ss)
betas = c(3,.2,.05,.5)
logiti = function(x){1/(1+exp(-x))}
#simu8 the data according to provided covariates
#in addition to one allele-specific count per individual it produces
#the allele-specific reads equally split among 8, 4 and 2 SNPS

  dat8 = simu8(niter, Xmatr=Xm$Xmatr, haplotype=Xm$thp, totmean=mn, dblcnt = 0, 
     percase = percase, phiNB = th, phiBB = th/dv, b0 = b0, b1 = 0, betas = betas)

     #dat8$asn - one count per individual
     #dat8$asn2 - data split in 1 SNPs
     #dat8$asn4 - data split in 4 SNPs
     #dat8$asn8 - data split in 8 SNPs

alpha = 0.05
if(!file.exists("resnew"))dir.create("resnew")
suff = sprintf("add%s_phi%s_div%s_iter%s_ss%s_nsnp", 
               b0, th, dv, niter, ss)
suff

#fitting the data in 8-SNP scenario
dat0 = dat8
dat0$asn = dat0$asn8S
dat0$asnp = dat0$asnp8S
dat0$haplotype = dat0$haplotype8S
nSNP = 8
st3 = proc.time()
en3 = proc.time()
en3[3]-st3[3]

#est hess
i = 1
#i = i + 1
st3 = proc.time()
resa = res8=res4=res2=res = matrix(0, nrow=niter, ncol=8)
for(i in 1:niter){
#i = 3
par0 = runif(2)
o1 = optim(par0, nllBB , n = dat8$asn[i,], y=dat8$asnp[i,], method="L-BFGS-B",hessian=T,lower=c(1e-10,1e-10), upper=(1-c(1e-10,1e-10)))
#hm = hessian(nllBB , x=o1$par, n = dat8$asn[i,], y=dat8$asnp[i,])
hm = fishBB(n = dat8$asn[i,], y=dat8$asnp[i,], rho=o1$par[1], pi1=o1$par[2])
res[i, 1:2] = o1$par
er = qr.solve(o1$hessian)
res[i, 3:4] = diag(er)
res[i,5] = er[1,2]
er = qr.solve(hm)
res[i, 6:7] = diag(er)
res[i,8] = er[1,2]


o1 = optim(par0, nllBB , n = dat8$asn2[i,], y=dat8$asnp2[i,], method="L-BFGS-B",hessian=T,lower=c(1e-10,1e-10), upper=(1-c(1e-10,1e-10)))
#hm = hessian(nllBB , x=o1$par, n = dat8$asn2[i,], y=dat8$asnp2[i,])
hm = fishBB(n = dat8$asn2[i,], y=dat8$asnp2[i,], rho=o1$par[1], pi1=o1$par[2])
res2[i, 1:2] = o1$par
er = qr.solve(o1$hessian)
res2[i, 3:4] = diag(er)
res2[i,5] = er[1,2]
er = qr.solve(hm)
res2[i, 6:7] = diag(er)
res2[i,8] = er[1,2]

o1 = optim(par0, nllBB , n = dat8$asn4[i,], y=dat8$asnp4[i,], method="L-BFGS-B",hessian=T,lower=c(1e-10,1e-10), upper=(1-c(1e-10,1e-10)))
#hm = hessian(nllBB , x=o1$par, n = dat8$asn4[i,], y=dat8$asnp4[i,])
hm = fishBB(n = dat8$asn4[i,], y=dat8$asnp4[i,], rho=o1$par[1], pi1=o1$par[2])
res4[i, 1:2] = o1$par
er = qr.solve(o1$hessian)
res4[i, 3:4] = diag(er)
res4[i,5] = er[1,2]
er = qr.solve(hm)
res4[i, 6:7] = diag(er)
res4[i,8] = er[1,2]

o1 = optim(par0, nllBB , n = dat8$asn8[i,], y=dat8$asnp8[i,], method="L-BFGS-B",hessian=T,lower=c(1e-10,1e-10), upper=(1-c(1e-10,1e-10)))
#hm = hessian(nllBB , x=o1$par, n = dat8$asn8[i,], y=dat8$asnp8[i,])
hm = fishBB(n = dat8$asn8[i,], y=dat8$asnp8[i,], rho=o1$par[1], pi1=o1$par[2])
res8[i, 1:2] = o1$par
er = qr.solve(o1$hessian)
res8[i, 3:4] = diag(er)
res8[i,5] = er[1,2]
er = qr.solve(hm)
res8[i, 6:7] = diag(er)
res8[i,8] = er[1,2]
if(i%%25==0)message(i)
}
en3 = proc.time()
en3[3]-st3[3]
res[1:4,]

write.table(res,  sprintf("hess_1SNP_%s.csv",mn), row.names=F, col.names=F, quote=F, sep=",")
write.table(res2, sprintf("hess_2SNP_%s.csv",mn), row.names=F, col.names=F, quote=F, sep=",")
write.table(res4, sprintf("hess_4SNP_%s.csv",mn), row.names=F, col.names=F, quote=F, sep=",")
write.table(res8, sprintf("hess_8SNP_%s.csv",mn), row.names=F, col.names=F, quote=F, sep=",")

#mn = 500;rho=.33;b0=0;logiti = function(x){1/(1+exp(-x))}
#mn = 100;rho=.33;b0=0;logiti = function(x){1/(1+exp(-x))}
res =  read.csv(sprintf("hess_1SNP_%s.csv",mn))
res2 = read.csv(sprintf("hess_2SNP_%s.csv",mn))
res4 = read.csv(sprintf("hess_4SNP_%s.csv",mn))
res8 = read.csv(sprintf("hess_8SNP_%s.csv",mn))

png("SSE_vs_modelSE.png", height=8, width=8, res=300, units="in")
par(mfrow=c(2,2))
plot(density(res[,1]))
lines(density(rnorm(1e5, mean=rho, sd=sqrt(median(res2[,3])))), lty=2)
plot(density(res2[,1]), col="green")
lines(density(rnorm(1e5, mean=rho, sd=sqrt(median(res2[,3])))), col="green", lty=2)
plot(density(res4[,1]), lty=1, col="blue")
lines(density(rnorm(1e5, mean=rho, sd=sqrt(median(res4[,3])))), col="blue", lty=2)
plot(density(res8[,1]), lty=1, col="red")
lines(density(rnorm(1e5, mean=rho, sd=sqrt(median(res8[,3])))), col="red", lty=2)
dev.off()

cols = c("black", "blue", "goldenrod", "red")
png("SSE_vs_modelSE_an.png", height=4, width=8, res=300, units="in")
par(mfrow=c(1,2))
plot(density(res[,1]), bty="n", main="OD", xlab="observed rho")
abline(v=rho, lty=3)
lines(density(res2[,1]), col=cols[2])
lines(density(res4[,1]), lty=1, col=cols[3])
lines(density(res8[,1]), lty=1, col=cols[4])
i = 6
modelSE=c(sqrt(median(res[,i])), sqrt(median(res2[,i])),
          sqrt(median(res4[,i])),sqrt(median(res8[,i])))
obserSE=c(sd(res[,1]), sd(res2[,1]),
          sd(res4[,1]),sd(res8[,1]))
legend("topleft", legend=c("sd(mod)", round(modelSE,3)), bty="n", text.col=c(cols[1], cols))
legend("topright", legend=c("sd(obs)", round(obserSE,3)), bty="n", text.col=c(cols[1], cols))

plot(density(res[,2]), bty="n", main="eQTL", xlab="observed pi")
abline(v=logiti(b0), lty=3)
lines(density(res2[,2]), col=cols[2])
lines(density(res4[,2]), lty=1, col=cols[3])
lines(density(res8[,2]), lty=1, col=cols[4])

i = 7
modelSE=c(sqrt(median(res[,i])), sqrt(median(res2[,i])),
          sqrt(median(res4[,i])),sqrt(median(res8[,i])))
obserSE=c(sd(res[,2]), sd(res2[,2]),
          sd(res4[,2]),sd(res8[,2]))
legend("topleft", legend=c("sd(mod)", round(modelSE,3)), bty="n", text.col=c(cols[1], cols))
legend("topright", legend=c("sd(obs)", round(obserSE,3)), bty="n", text.col=c(cols[1], cols))
dev.off()

c(var(res[,1]), median(res[,3]))
c(var(res2[,1]), median(res2[,3]))
c(var(res4[,1]), median(res4[,3]))
c(var(res8[,1]), median(res8[,3]))





cols = c("black", "blue", "goldenrod", "red")
png("SSE_vs_modelSE_an2.png", height=4, width=8, res=300, units="in")
par(mfrow=c(1,2))
plot(density(res[,1]), bty="n", main="OD", xlab="observed rho")
abline(v=rho, lty=3)
lines(density(res2[,1]), col=cols[2])
lines(density(res4[,1]), lty=1, col=cols[3])
lines(density(res8[,1]), lty=1, col=cols[4])
i = 6
modelSE=c(sqrt(median(res[,i])), sqrt(median(res2[,i])),
          sqrt(median(res4[,i])),sqrt(median(res8[,i])))
obserSE=c(sd(res[,1]), sd(res2[,1]),
          sd(res4[,1]),sd(res8[,1]))
legend("topleft", legend=c("", "1 SNP", "2 SNP", "4 SNP", "8 SNP"), bty="n", text.col=c(cols[1], cols))
legend("topright", legend=c("sd:mod, obs", sprintf("%s, %s", round(modelSE,3), round(obserSE,3))), bty="n", text.col=c(cols[1], cols))

plot(density(res[,2]), bty="n", main="eQTL", xlab="observed pi")
abline(v=logiti(b0), lty=3)
lines(density(res2[,2]), col=cols[2])
lines(density(res4[,2]), lty=1, col=cols[3])
lines(density(res8[,2]), lty=1, col=cols[4])

i = 7
modelSE=c(sqrt(median(res[,i])), sqrt(median(res2[,i])),
          sqrt(median(res4[,i])),sqrt(median(res8[,i])))
obserSE=c(sd(res[,2]), sd(res2[,2]),
          sd(res4[,2]),sd(res8[,2]))
legend("topleft", legend=c("", "1 SNP", "2 SNP", "4 SNP", "8 SNP"), bty="n", text.col=c(cols[1], cols))
legend("topright", legend=c("sd:mod, obs", sprintf("%s, %s", round(modelSE,3), round(obserSE,3))), bty="n", text.col=c(cols[1], cols))
dev.off()
 

q("no")
