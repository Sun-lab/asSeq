
# --------------------------------------------------------------------
# check the MLE and its asymtotic distribution when the data 
# were simulated with some clustering structure. Within each 
# cluster, counts following binomial distribution, and between
# clusters, they follow beta-binomial distirbution. 
# --------------------------------------------------------------------

source("../asSeq/_test/aseR.R")

nrep = 500 # number of replicates
nind = 200 # sample size
mu0  = 20  # total number of allele-specific reads
nSNP = 4   # split the allele-specific read counts tin nSNP SNPs.
prob = 0.5 # success prob. of beta-binomial
bbod = 0.5 # over-dispersion of beta-binomial theta = 1/(alpha + beta)
rho0 = 1/(1/bbod + 1) # rho = 1/(1 + alpha + beta), hence theta = 1/(1/rho - 1)
rho0

alphaBB = prob*(1-rho0)/rho0
betaBB  = (1-prob)*(1-rho0)/rho0

mus = theta0s = theta1s = logL0s = logL1s = matrix(nrow=nrep, ncol=2)

for(i in 1:nrep){
  if(i %% 100 == 0){ cat(i, " ", date(), "\n")}
  # asrecT is total ASReC per individaul for all individuals
  # total is the total ASReC per SNP and per individual
  asrecT = rnbinom(n=nind, size=1e6, mu=mu0)
  total1 = round(asrecT/nSNP)
  total2 = asrecT - (nSNP - 1)*total1
  total  = matrix(c(rep(total1, times=nSNP-1), total2), ncol=nSNP)

  while(any(total == 0)){
    asrecT = rnbinom(n=nind, size=1e6, mu=mu0)
    total1 = round(asrecT/nSNP)
    total2 = asrecT - (nSNP - 1)*total1
    total  = matrix(c(rep(total1, times=nSNP-1), total2), ncol=nSNP)
  }
  
  # hap1 is the number of ASReC for haplotype 1
  hap1 = matrix(nrow=nind, ncol=nSNP)
  prb1 = rbeta(nind, alphaBB, betaBB)
  
  for(k in 1:nSNP){
    hap1[,k] = rbinom(nind, size=total[,k], prob=prb1)
  }
  
  nTotal = rowSums(total)
  nA     = rowSums(hap1)
  zeta   = rbinom(nind, 1, 0.5)
  
  fit1 = aseR(nA, nTotal, zeta)
  theta0s[i,1] = fit1$parH0[1]
  theta1s[i,1] = fit1$parH1[1]
  mus[i,1]     = fit1$parH1[2]
  logL0s[i,1]  = fit1$logLikH0
  logL1s[i,1]  = fit1$logLikH1
  
  
  grad1 = gradLogH1(fit1$parH1, nA, nTotal, zeta)
    
    
  nTotal = c(total)
  nA     = c(hap1)
  zeta1  = rep(zeta, each=nSNP)
  
  fit2 = aseR(nA, nTotal, zeta1)
  theta0s[i,2] = fit2$parH0[1]
  theta1s[i,2] = fit2$parH1[1]
  mus[i,2]     = fit2$parH1[2]
  logL0s[i,2]  = fit2$logLikH0
  logL1s[i,2]  = fit2$logLikH1
}

pdf("figures/check_beta_binomial_MLE.pdf", width=9, height=6)
par(mfrow=c(2,3), bty="n", mar=c(5,4,3,1), cex=0.8)

plot(density(mus[,1]), col="darkred", lty=2, main="mu", ylim=c(0,20))
lines(density(mus[,2]), col="darkblue", lty=3)

plot(density(theta0s[,1]), col="darkred", lty=2, main="theta, H0")
lines(density(theta0s[,2]), col="darkblue", lty=3)

plot(density(theta1s[,1]), col="darkred", lty=2, main="theta, H1")
lines(density(theta1s[,2]), col="darkblue", lty=3)

lrt1 = 2*(logL1s[,1] - logL0s[,1])
lrt2 = 2*(logL1s[,2] - logL0s[,2])

plot(density(lrt1), col="darkred", lty=2, main="LRT")
lines(density(lrt2), col="darkblue", lty=3)
x1 = seq(0.01, 8, by=0.01)
lines(x1, dchisq(x1, 1), lwd=1.5, col="grey")

plot(lrt1, lrt2)
abline(0, 1, col="red")
dev.off()
