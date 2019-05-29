setwd("/fh/fast/sun_w/licai/_tumor_eQTL/GitHub/asSeq/asSeq2/_test/")
source("./tumor_eQTL.R")
sourceCpp("../asSeq2/src/tumor_eQTL.cpp")

load("_test_func_tumorEQTL.Rdata")

N
dat[1:5,]
ETA
GAMMA
KAPPA
lgy1 = dat$LGX1
zz = zz2 = ZZ[,1]
table(zz)
zz[which(zz == 2)] = 1
zz[which(zz == 3)] = 2
table(zz)

X = data.matrix(XX)
betas = BETA
phi = PHI
y = dat$total

theta = PSI
ind = which(dat$total_phased > 7)
zz_AS = zz2[ind]
ni = dat$total_phased[ind]
ni0 = dat$hapB[ind]
lbc = dat$LBC[ind]
tauB = dat$tauB[ind]
tau = dat$tau[ind]
RHO_AS = RHO[ind]

offsets = rep(0, N)
mu = rep(0, N)
expXbeta = rep(0, N)

#-------------------------------------------------------------
# TReC model (negtive binomial)
#-------------------------------------------------------------

Roffset = compute_offset(zz2, RHO, KAPPA, ETA, GAMMA, tau1, tau2)
RcppT_compute_offset(zz2, RHO, KAPPA, ETA, GAMMA, tau1, tau2, offsets)
summary(Roffset - offsets)
RcppT_compute_expXbeta(X, betas, expXbeta)
expXbetaR = exp(X %*% betas)
summary(as.vector(expXbeta - expXbetaR))

# regression
RcppT_reg_LL(y, X, offsets, betas, lgy1, mu) 
RcppT_reg_grad(y, X, mu, betas) 
reg  = RcppT_reg_BFGS(y, X, offsets, rep(0, length(betas)+1), lgy1)
regR = update_beta_phi(y, X , offsets)
reg$PAR
c(regR$betas, log(regR$phi))

# eQTL
para = log(c(KAPPA, ETA, GAMMA))
H0 = 0
loglikNB(para, H0, y, zz2, X, betas, phi, RHO, tau1, tau2)
RcppT_compute_offset(zz2, RHO, KAPPA, ETA, GAMMA, tau1, tau2, offsets)
RcppT_compute_expXbeta(X, betas, expXbeta) 
Grad_NB(para, H0, y, zz2, X, betas, phi, RHO, tau1, tau2)
RcppT_loglikNB_KEG(para, H0, y, zz2, phi, RHO, tau1, tau2, lgy1, 
                    expXbeta, offsets, mu)
RcppT_grad_NB(para, H0, y, zz2, RHO, phi, tau1, tau2, expXbeta, 
              offsets, mu) 
mu = expXbeta*exp(offsets)
RcppT_loglikNB(y, phi, lgy1, mu) 

H0 = 1 
para = log(c(KAPPA, GAMMA))
loglikNB(para, H0, y, zz2, X, betas, phi, RHO, tau1, tau2)
RcppT_compute_offset(zz2, RHO, KAPPA, 1, GAMMA, tau1, tau2, offsets)
RcppT_compute_expXbeta(X, betas, expXbeta) 
Grad_NB(para, H0, y, zz2, X, betas, phi, RHO, tau1, tau2)
RcppT_loglikNB_KEG(para, H0, y, zz2, phi, RHO, tau1, tau2, lgy1, 
                    expXbeta, offsets, mu)
RcppT_grad_NB(para, H0, y, zz2, RHO, phi, tau1, tau2, expXbeta, 
              offsets, mu) 
mu = expXbeta*exp(offsets)
RcppT_loglikNB(y, phi, lgy1, mu) 

H0 = 2
para = log(c(KAPPA, ETA))
loglikNB(para, H0, y,zz2, X, betas, phi, RHO, tau1, tau2)
RcppT_compute_offset(zz2, RHO, KAPPA, ETA, 1, tau1, tau2, offsets)
RcppT_compute_expXbeta(X, betas, expXbeta) 
Grad_NB(para, H0, y, zz2, X, betas, phi, RHO, tau1, tau2)
RcppT_loglikNB_KEG(para, H0, y, zz2, phi, RHO, tau1, tau2, lgy1, 
                    expXbeta, offsets, mu)
RcppT_grad_NB(para, H0, y, zz2, RHO, phi, tau1, tau2, expXbeta, 
              offsets, mu) 
mu = expXbeta*exp(offsets)
RcppT_loglikNB(y, phi, lgy1, mu) 

H0 = 0
para = c(0,0,0)
RcppT_trec_KEG_BFGS(para, H0, y, zz2, RHO, X, betas, phi,
                     tau1, tau2, lgy1)
update_keg_trec(H0, para, y,zz2, X, betas, phi, RHO, tau1, tau2)
sfitR = TReC_sfit(H0, para, y,zz2, X, RHO, tau1, tau2)
sfit  = RcppT_trec_sfit(H0, para, y, zz2, RHO, X, tau1, tau2, lgy1) 
c(sfitR$KAPPA, sfitR$ETA, sfitR$GAMMA)
exp(sfit$PAR)
sfitR$loglik
sfit$LL
c(sfitR$betas, sfitR$phi)
c(sfit$reg_par[1:length(betas)], exp(sfit$reg_par[length(betas)+1]))

H0 = 1
para = c(0,0)
RcppT_trec_KEG_BFGS(para, H0, y, zz2, RHO, X, betas, phi, tau1, tau2, lgy1)
update_keg_trec(H0,para, y,zz2, X,betas, phi, RHO, tau1, tau2)
sfitR = TReC_sfit(H0, para, y,zz2, X, RHO, tau1, tau2)
sfit  = RcppT_trec_sfit(H0, para, y, zz2, RHO, X, tau1, tau2, lgy1) 
c(sfitR$KAPPA, sfitR$ETA, sfitR$GAMMA)
exp(sfit$PAR)
sfitR$loglik
sfit$LL
c(sfitR$betas, sfitR$phi)
c(sfit$reg_par[1:length(betas)], exp(sfit$reg_par[length(betas)+1]))


H0 = 2
RcppT_trec_KEG_BFGS(para, H0, y, zz2, RHO, X, betas, phi, tau1, tau2, lgy1)
update_keg_trec(H0, para, y,zz2, X,betas, phi, RHO, tau1, tau2)
sfitR = TReC_sfit(H0, para, y,zz2, X, RHO, tau1, tau2)
sfit  = RcppT_trec_sfit(H0, para, y, zz2, RHO, X, tau1, tau2, lgy1) 
c(sfitR$KAPPA, sfitR$ETA, sfitR$GAMMA)
exp(sfit$PAR)
sfitR$loglik
sfit$LL
c(sfitR$betas, sfitR$phi)
c(sfit$reg_par[1:length(betas)], exp(sfit$reg_par[length(betas)+1]))


TReCR = TReC_test(y,zz2, X, RHO, tau1, tau2)
TReC  = RcppT_trec(y,zz2, RHO, X, tau1, tau2, lgy1)
c(TReCR$p.eta, TReCR$p.gamma)
c(TReC$p_eta, TReC$p_gamma)
c(TReCR$loglik.full, TReCR$loglik.eta, TReCR$loglik.gamma)
c(TReC$LL, TReC$LL_eta, TReC$LL_gamma)
c(TReCR$KAPPA, TReCR$ETA, TReCR$GAMMA)
exp(TReC$PAR)
c(TReCR$betas, TReCR$phi)
c(TReC$reg_par[1:length(betas)], exp(TReC$reg_par[length(betas)+1]))

#-------------------------------------------------------------
# ASE likelihood (beta binomial)
#-------------------------------------------------------------
pis = rep(0, length(zz_AS))

piR = compute_pi(zz_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau)
RcppT_compite_pi(zz_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau, pis)
summary(piR - pis)

H0 = 0
para = log(c(KAPPA, ETA, GAMMA))
loglikBB(para, H0, ni0, ni, theta, zz_AS, RHO_AS, tauB, tau)
RcppT_loglikBB_THETA(ni, ni0, log(theta), pis, lbc) 
loglikBB_THETA(theta, ni0, ni, pis)

RcppT_grad_BB_THETA(ni, ni0, log(theta), pis)
lASE.dTHETA(theta, ni0, ni, pis)*theta

RcppT_grad_BB_KEG(para, H0, zz_AS, RHO_AS, ni, ni0, log(theta), lbc, tauB, 
                  tau, pis)
Grad_BB_keg(para, H0,ni0,ni,theta, zz_AS, RHO_AS,tauB, tau)
RcppT_loglikBB_KEG(para, H0, zz_AS, RHO_AS, ni, ni0, log(theta), lbc, tauB, 
                   tau, pis) 


