setwd('/fh/fast/sun_w/licai/_tumor_eQTL/GitHub/asSeq/asSeq2/_test')

Rcpp::sourceCpp('../Rcpp/asSeq_test1.cpp')
source("trecR.R")
load('_test_func.Rdata')

N
dat[1:5,]
bxj
lgy1 = dat$LGX1
zz = ZZ[,1]
table(zz)
zz[which(zz == 2)] = 1
zz[which(zz == 3)] = 2
table(zz)
X = data.matrix(XX)
betas = BETA
phi = PHI 
y = dat$total
b0 = 0

offsets = rep(0, N)
mu1 = rep(0, N)
mu = exp(X%*%betas)

## negtive binomial
Rcpp_logLTReC(bxj, lgy1, y, zz, mu, phi, T, mu1, offsets)
logLTReC(bxj, y, zz, mu, b0, phi, 'negbin')

Rcpp_grad_bxj_trec(bxj, y, zz, mu1, phi, T)
grad.bxj.trec(bxj, y, zz, mu, b0, phi, 'negbin')

Rcpp_reg_LL(y, X, offsets, PARAMS = c(betas, log(phi)), fam_nb = T, lgy1, mu)
Rcpp_reg_grad(y, X, mu, PARAMS = c(betas, log(phi)), fam_nb = T)
Rcpp_reg_Hess(y, X, mu, PARAMS = c(betas, log(phi)), fam_nb = T)
Rcpp_reg(y, X, offsets, rep(0.1, 5),fam_nb = T, lgy1)
Rcpp_reg_BFGS(y, X, offsets, rep(0.1, 5),fam_nb = T, lgy1) 
library(MASS)
g1  = glm.nb(y ~ X+offset(offsets))
g1
log(1/g1$theta)
g1$twologlik/2

optim(b0, fn=logLTReC, gr=grad.bxj.trec, y=y, x=zz, mu=g1$fitted, b0=b0, phi=0.1, 
      fam='negbin', method="BFGS", control=list(fnscale=-1.0, trace=0))


Rcpp_reg_LL(y, X, offsets, PARAMS = c(betas),fam_nb = F,lgy1, mu)
Rcpp_reg_grad(y, X, mu, PARAMS = c(betas),fam_nb = F)
Rcpp_reg_Hess(y, X, mu, PARAMS = c(betas),fam_nb = F)
Rcpp_reg(y,X,offsets,rep(0,4),fam_nb = F, lgy1)

g1 = glm(y~X+offset(offsets), family = 'poisson')
g1
logLik(g1)


