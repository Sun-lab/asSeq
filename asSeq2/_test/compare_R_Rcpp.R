setwd('/fh/fast/sun_w/licai/_tumor_eQTL/GitHub/asSeq/asSeq2/_test')
library(Rcpp)
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

Cop = Rcpp_trec_bxj(y, X, 0, zz, BETA, phi, T, lgy1, max_iter = 4000L, 
              eps = 1e-07, show = TRUE) 
Rop = optim(0, fn=logLTReC, gr=grad.bxj.trec, y=y, x=zz, mu=exp(X%*%BETA), b0=b0, phi=phi, 
      fam='negbin', method="BFGS", control=list(fnscale=-1.0, trace=0))
Cop
Rop

c(Cop$bxj, Rop$par)

Rcpp_logLTReC(Cop$bxj, y, X, zz, BETA, phi, T, lgy1, mu) 
Rcpp_logLTReC(Rop$par, y, X, zz, BETA, phi, T, lgy1, mu) 
Rcpp_logLTReC((Rop$par+Cop$bxj)/2, y, X, zz, BETA, phi, T, lgy1, mu) 

sapply(seq(Rop$par, Cop$bxj, by = 1e-7), Rcpp_logLTReC, y, X, zz, BETA, phi, T, lgy1, mu)
sapply(seq(Rop$par, Cop$bxj, by = 1e-7), Rcpp_grad_hess_bxj_trec, y, zz, mu, phi, T)


Ctrec = Rcpp_trec(y, X, zz, T, lgy1)
Ctrec
Rtrec = trecR(y, X, zz, 'negbin', yfit = T)
Rtrec

Ctrec$bxj
Rtrec$b

Ctrec$reg_par[1:4]
Rtrec$betas

exp(Ctrec$reg_par[5])
Rtrec$phi

Rcpp_logLTReC(Ctrec$bxj, y, X, zz, Ctrec$reg_par[1:4], exp(Ctrec$reg_par[5]), T, lgy1, mu) 
Rcpp_logLTReC(Rtrec$b, y, X, zz, Rtrec$betas, Rtrec$phi, T, lgy1, mu) 

Rcpp_grad_hess_bxj_trec(Ctrec$bxj, y, zz, mu, exp(Ctrec$reg_par[5]), T)
Rcpp_grad_hess_bxj_trec(Rtrec$b, y, zz, mu, Rtrec$phi, T)


q('no')




## poisson
Rcpp_logLTReC(bxj, y, X, zz, betas, phi, F, lgy1, mu) 
logLTReC(bxj, y, zz, exp(X%*%betas), b0, phi, 'poisson')

Rcpp_grad_hess_bxj_trec(bxj, y, zz, mu, 0, F)
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
source("aseR.R")

theta = PSI
Pi = exp(bxj)/(exp(bxj) + 1)
ind = which(dat$total_phased > 7)
zz_AS = zz[ind]
ni = dat$total_phased[ind]
ni0 = dat$hapB[ind]
lbc = dat$LBC[ind]
zeta = zz_AS==1

mat = NULL

for(pis in c(1:9*0.1)){
  aa  = Rcpp_loglikBB(ni, ni0, pis, log(theta),lbc, zeta)
  bb  = logH1(c(theta, pis), ni0, ni, zeta)
  mat = rbind(mat, c(aa,bb))
}
mat
mat[,2] - mat[,1]

Rcpp_loglikBB(ni, ni0, 0.5, log(theta),lbc, zeta)
logH0(c(theta, pis), ni0, ni)
logH1(c(theta, 0.5), ni0, ni, zeta)

Rcpp_ase_grad(ni, ni0, Pi, log(theta), zeta)
gradLogH1(c(theta, Pi), ni0, ni, zeta)

Rcpp_ase_grad_H0(ni, ni0, 0.5, log(theta), zeta)
gradLogH0(c(theta, 0.5), ni0, ni)

# Rcpp_ase_hess(ni, ni0, 0.51602267, 0.06677156, zeta) 

par0 = 0.1
Rcpp_ase_BFGS(ni, ni0, zeta, 0, log(par0), lbc, max_iter = 4000L, 
              eps = 1e-5, show = TRUE) #failed in search ???
op0  = optim(par0, logH0, gr=gradLogH0, nA=ni0, nTotal=ni, 
             method="L-BFGS-B", lower=1e-16, upper=Inf, 
             control=list(fnscale=-1))
op0

Rcpp_ase_BFGS(ni, ni0, zeta, 1, c(0.5, 0.1), lbc, max_iter = 4000L, 
              eps = 1e-5, show = TRUE)
optim(c(op0$par, 0.5), logH1, gr=gradLogH1, nA=ni0, nTotal=ni, zeta=zeta, 
      method="L-BFGS-B", lower=c(0,0) + 1e-16, upper=c(Inf, 1-1e-16), 
      control=list(fnscale=-1),  hessian = T)


Rcpp_ase(ni, ni0, zeta, lbc,show = F) 
aseR(ni0, ni, zeta)

#-------------------------------------------------------------
# TReCASE
#-------------------------------------------------------------
source("trecaseR.R")
mu = rep(0, length(y))
Rcpp_trecase_LL(bxj, y, X, z, BETA, phi, 1, lgy1, mu, ni, ni0, log(theta), 
                lbc, zeta) 
loglikJoin(bxj, y, z, ni0, ni, zeta, exp(X%*%BETA), 0, phi, theta)

Rcpp_trecase_grad_bxj(bxj, y, X, z, BETA, phi, 1, lgy1, mu, ni, ni0, log(theta), 
                      lbc, zeta)
grad.bxj(bxj, y, z, ni0, ni, zeta, exp(X%*%BETA), 0, phi, theta)

Rcpp_trecase_BFGS(0, y, X, z, BETA, phi, T, lgy1, ni, ni0, log(theta), 
                  lbc, zeta, max_iter = 4000L, eps = 1e-07, show = TRUE) 
op = optim(par=0, fn=loglikJoin, gr=grad.bxj, y=y, x=z, nA=ni0, nTotal=ni, 
           zeta=zeta, mu=exp(X%*%BETA), b0=0, phi=phi, theta=theta, method = "BFGS", 
           control=list(fnscale=-1.0))
op
