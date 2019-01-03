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
mu = rep(0, N)

#-------------------------------------------------------------
# TReC model
#-------------------------------------------------------------

## negtive binomial
Rcpp_logLTReC(bxj, y, X, zz, betas, phi, T, lgy1, mu) 
logLTReC(bxj, y, zz, exp(X%*%betas), b0, phi, 'negbin')

Rcpp_grad_hess_bxj_trec(bxj, y, zz, mu, phi, T)
grad.bxj.trec(bxj, y, zz, exp(X%*%betas), b0, phi, 'negbin')

Rcpp_reg_LL(y, X, offsets, PARAMS = c(betas, log(phi)), fam_nb = T, lgy1, mu)
Rcpp_reg_grad(y, X, mu, PARAMS = c(betas, log(phi)), fam_nb = T)
Rcpp_reg_Hess(y, X, mu, PARAMS = c(betas, log(phi)), fam_nb = T)
Rcpp_reg(y, X, offsets, rep(0.1, 5),fam_nb = T, lgy1)
Rcpp_reg_BFGS(y, X, offsets, rep(0.1, 5),fam_nb = T, lgy1) 
library(MASS)

g1  = glm.nb(y ~ X+offset(offsets))
BETA = summary(g1)$coef[,1]
BETA
phi = 1/g1$theta
log(phi)
g1$twologlik/2

Rcpp_trec_bxj(y, X, 0, zz, BETA, phi, T, lgy1, max_iter = 4000L, 
              eps = 1e-07, show = TRUE) 
optim(0, fn=logLTReC, gr=grad.bxj.trec, y=y, x=zz, mu=exp(X%*%BETA), b0=b0, phi=phi, 
      fam='negbin', method="BFGS", control=list(fnscale=-1.0, trace=0))

Rcpp_trec(y, X, zz, T, lgy1)
trecR(y, X, zz, 'negbin')


## poisson
Rcpp_logLTReC(bxj, y, X, zz, betas, phi, F, lgy1, mu) 
logLTReC(bxj, y, zz, exp(X%*%betas), b0, phi, 'poisson')

Rcpp_grad_bxj_trec(bxj, y, zz, mu, phi, F)
grad.bxj.trec(bxj, y, zz, exp(X%*%betas), b0, 0, 'poisson')

Rcpp_reg_LL(y, X, offsets, PARAMS = c(betas, log(phi)), fam_nb = T, lgy1, mu)
Rcpp_reg_grad(y, X, mu, PARAMS = c(betas, log(phi)), fam_nb = T)
Rcpp_reg_Hess(y, X, mu, PARAMS = c(betas, log(phi)), fam_nb = T)
Rcpp_reg(y, X, offsets, rep(0.1, 5),fam_nb = T, lgy1)
Rcpp_reg_BFGS(y, X, offsets, rep(0.1, 5),fam_nb = T, lgy1) 
g1 = glm(y~X+offset(offsets), family = 'poisson')
BETA = summary(g1)$coef[,1]
BETA
logLik(g1)

Rcpp_trec_bxj(y, X, 0, zz, BETA, phi, F, lgy1) 
optim(0, fn=logLTReC, gr=grad.bxj.trec, y=y, x=zz, mu=exp(X%*%BETA), b0=b0, phi=phi, 
      fam='poisson', method="BFGS", control=list(fnscale=-1.0, trace=0))

Rcpp_trec(y, X, zz, F, lgy1)
trecR(y, X, zz, 'poisson')

#-------------------------------------------------------------
# ASE model
#-------------------------------------------------------------
mat = NULL
for(pis in c(1:9*0.1)){
  aa=Rcpp_loglikBB(dat$total_phased, dat$hapB, pis, PSI, dat$LBC)
  bb=sum(logBB1(dat$total_phased, dat$hapB, pis, PSI))
  cc=Rcpp_vec_log_BB( dat$hapB,dat$total_phased, pis, PSI, dat$LBC)
  mat = rbind(mat, c(aa,bb,cc))
}


