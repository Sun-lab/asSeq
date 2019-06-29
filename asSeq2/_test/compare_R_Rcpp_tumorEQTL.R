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
reg  = RcppT_reg_BFGS(y, X, offsets, rep(0, length(betas)+1), lgy1,
                    max_iter = 400, eps=1e-5, show = F)
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
Rcpp_loglikNB(y, phi, lgy1, mu) 

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
Rcpp_loglikNB(y, phi, lgy1, mu) 

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
Rcpp_loglikNB(y, phi, lgy1, mu) 

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
RcppT_compite_pi(zz_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau, pis)
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

RcppT_compite_pi(zz_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau, pis)
bfgsC = RcppT_ase_theta_BFGS(-4, ni, ni0, pis, lbc) 
bfgsR = update_theta(para, H0, zz_AS,RHO_AS, ni0, ni, 0.1, tauB, tau)
c(bfgsC$LL, bfgsR$value)
c(exp(bfgsC$PAR), bfgsR$par)

RcppT_ase_KEG_BFGS(para, H0, zz_AS, RHO_AS, ni0, ni, log(theta), tauB, 
                              tau, lbc)
update_keg_ase(para, H0, zz_AS, RHO_AS,  ni0, ni,theta, tauB, tau)



H0 = 1
para = log(c(KAPPA, GAMMA))
RcppT_compite_pi(zz_AS, RHO_AS, KAPPA, 1, GAMMA, tauB, tau, pis)
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

bfgsC = RcppT_ase_theta_BFGS(-4, ni, ni0, pis, lbc) 
bfgsR = update_theta(para, H0, zz_AS,RHO_AS, ni0, ni, 0.1, tauB, tau)
c(bfgsC$LL, bfgsR$value)
c(exp(bfgsC$PAR), bfgsR$par)


RcppT_ase_KEG_BFGS(para, H0, zz_AS, RHO_AS, ni0, ni, log(theta), tauB, 
                   tau, lbc)
update_keg_ase(para, H0, zz_AS, RHO_AS,  ni0, ni,theta, tauB, tau)


H0 = 2
para = log(c(KAPPA, ETA))
RcppT_compite_pi(zz_AS, RHO_AS, KAPPA, ETA, 1, tauB, tau, pis)
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

bfgsC = RcppT_ase_theta_BFGS(-4, ni, ni0, pis, lbc) 
bfgsR = update_theta(para, H0, zz_AS,RHO_AS, ni0, ni, 0.1, tauB, tau)
c(bfgsC$LL, bfgsR$value)
c(exp(bfgsC$PAR), bfgsR$par)


RcppT_ase_KEG_BFGS(para, H0, zz_AS, RHO_AS, ni0, ni, log(theta), tauB, 
                   tau, lbc)
update_keg_ase(para, H0, zz_AS, RHO_AS,  ni0, ni,theta, tauB, tau)

aseC = RcppT_ase(zz_AS, RHO_AS, ni0, ni, tauB, tau, lbc)
aseR = ASE_test(zz_AS, RHO_AS, ni0, ni,0.05, tauB, tau)
exp(c(aseC$PAR))
c(aseR$KAPPA, aseR$ETA, aseR$GAMMA)
c(aseC$p_eta, aseC$p_gamma)
c(aseR$p.eta, aseR$p.gamma)
c(aseC$LL, aseC$LL_eta, aseC$LL_gamma)
c(aseR$loglik.full, aseR$loglik.eta, aseR$loglik.gamma)



#-------------------------------------------------------------
# TReCASE
#-------------------------------------------------------------

H0 = 0
para = log(c(KAPPA, ETA, GAMMA))
RcppT_TReCASE_LL_KEG(para, H0, y, z, zz_AS, phi, RHO, RHO_AS, tau1, tau2, 
          ni0, ni, log(theta), tauB, tau, lgy1, expXbeta, lbc, offsets, 
          pis, mu) 
RcppT_TReCASE_LL(y, phi, lgy1, mu, ni0, ni, log(theta), pis, lbc) 
loglik.TReCASE(para, H0, y,z, zz_AS, X, betas, phi, RHO, 
               RHO_AS, tau1, tau2,  ni0, ni, theta, tauB, tau)
RcppT_TReCASE_grad_KEG(para, H0, y, z, zz_AS, phi, RHO, RHO_AS, tau1, tau2, 
                     ni0, ni, log(theta), tauB, tau, lgy1, expXbeta, lbc, offsets, 
                     pis, mu) 
grad.TReCASE.keg(para, H0, y,z, zz_AS, X, betas, phi, RHO,
                             RHO_AS, tau1, tau2, ni0, ni, theta, tauB, tau)

H0 = 1
para = log(c(KAPPA, GAMMA))
RcppT_TReCASE_LL_KEG(para, H0, y, z, zz_AS, phi, RHO, RHO_AS, tau1, tau2, 
                     ni0, ni, log(theta), tauB, tau, lgy1, expXbeta, lbc, offsets, 
                     pis, mu) 
RcppT_TReCASE_LL(y, phi, lgy1, mu, ni0, ni, log(theta), pis, lbc) 
loglik.TReCASE(para, H0, y,z, zz_AS, X, betas, phi, RHO, 
               RHO_AS, tau1, tau2,  ni0, ni, theta, tauB, tau)
RcppT_TReCASE_grad_KEG(para, H0, y, z, zz_AS, phi, RHO, RHO_AS, tau1, tau2, 
                       ni0, ni, log(theta), tauB, tau, lgy1, expXbeta, lbc, offsets, 
                       pis, mu) 
grad.TReCASE.keg(para, H0, y,z, zz_AS, X, betas, phi, RHO,
                 RHO_AS, tau1, tau2, ni0, ni, theta, tauB, tau)


H0 = 2
para = log(c(KAPPA, ETA))
RcppT_TReCASE_LL_KEG(para, H0, y, z, zz_AS, phi, RHO, RHO_AS, tau1, tau2, 
                     ni0, ni, log(theta), tauB, tau, lgy1, expXbeta, lbc, offsets, 
                     pis, mu) 
RcppT_TReCASE_LL(y, phi, lgy1, mu, ni0, ni, log(theta), pis, lbc) 
loglik.TReCASE(para, H0, y,z, zz_AS, X, betas, phi, RHO, 
               RHO_AS, tau1, tau2,  ni0, ni, theta, tauB, tau)
RcppT_TReCASE_grad_KEG(para, H0, y, z, zz_AS, phi, RHO, RHO_AS, tau1, tau2, 
                       ni0, ni, log(theta), tauB, tau, lgy1, expXbeta, lbc, offsets, 
                       pis, mu) 
grad.TReCASE.keg(para, H0, y,z, zz_AS, X, betas, phi, RHO,
                 RHO_AS, tau1, tau2, ni0, ni, theta, tauB, tau)


H0 = 0
para0 = c(0,0,0)
RcppT_trecase_KEG_BFGS(para0, H0, y, z, zz_AS, RHO, RHO_AS, X, BETA, phi, tau1, 
                       tau2, lgy1, ni0, ni, log(theta), tauB, tau, lbc)
update_keg_trecase(para0, H0, y,z, zz_AS, X, betas, phi, RHO,
                   RHO_AS, tau1, tau2, ni0, ni, theta, tauB, tau)
TReCASE_sfitC = RcppT_trecase_sfit(H0, para0, y, z, zz_AS, RHO, RHO_AS, X,
                                   tau1, tau2, lgy1, ni0, ni, tauB, tau, lbc) 
TReCASE_sfitR = TReCASE_sfit(H0, para0, y, z, zz_AS, X, RHO, RHO_AS, 
             tau1, tau2, ni0, ni, tauB, tau, maxiter =500, tol = 1e-7)
c(exp(TReCASE_sfitC$PAR))
c(TReCASE_sfitR$KAPPA, TReCASE_sfitR$ETA, TReCASE_sfitR$GAMMA)
c(TReCASE_sfitC$reg_par)
c(TReCASE_sfitR$betas, log(TReCASE_sfitR$phi))
c(TReCASE_sfitC$LL, TReCASE_sfitR$loglik)
c(TReCASE_sfitC$log_theta, log(TReCASE_sfitR$THETA))


H0 = 1
para0 = c(0,0)
RcppT_trecase_KEG_BFGS(para0, H0, y, z, zz_AS, RHO, RHO_AS, X, BETA, phi, tau1, 
                       tau2, lgy1, ni0, ni, log(theta), tauB, tau, lbc)
update_keg_trecase(para0, H0, y,z, zz_AS, X, betas, phi, RHO,
                   RHO_AS, tau1, tau2, ni0, ni, theta, tauB, tau)
TReCASE_sfitC = RcppT_trecase_sfit(H0, para0, y, z, zz_AS, RHO, RHO_AS, X,
                                   tau1, tau2, lgy1, ni0, ni, tauB, tau, lbc) 
TReCASE_sfitR = TReCASE_sfit(H0, para0, y, z, zz_AS, X, RHO, RHO_AS, 
                             tau1, tau2, ni0, ni, tauB, tau, maxiter =500, tol = 1e-7)
c(exp(TReCASE_sfitC$PAR))
c(TReCASE_sfitR$KAPPA, TReCASE_sfitR$ETA, TReCASE_sfitR$GAMMA)
c(TReCASE_sfitC$reg_par)
c(TReCASE_sfitR$betas, log(TReCASE_sfitR$phi))
c(TReCASE_sfitC$LL, TReCASE_sfitR$loglik)
c(TReCASE_sfitC$log_theta, log(TReCASE_sfitR$THETA))

H0 = 2
para0 = c(0,0)
RcppT_trecase_KEG_BFGS(para0, H0, y, z, zz_AS, RHO, RHO_AS, X, BETA, phi, tau1, 
                       tau2, lgy1, ni0, ni, log(theta), tauB, tau, lbc)
update_keg_trecase(para0, H0, y,z, zz_AS, X, betas, phi, RHO,
                   RHO_AS, tau1, tau2, ni0, ni, theta, tauB, tau)
TReCASE_sfitC = RcppT_trecase_sfit(H0, para0, y, z, zz_AS, RHO, RHO_AS, X,
                                   tau1, tau2, lgy1, ni0, ni, tauB, tau, lbc) 
TReCASE_sfitR = TReCASE_sfit(H0, para0, y, z, zz_AS, X, RHO, RHO_AS, 
                             tau1, tau2, ni0, ni, tauB, tau, maxiter =500, tol = 1e-7)
c(exp(TReCASE_sfitC$PAR))
c(TReCASE_sfitR$KAPPA, TReCASE_sfitR$ETA, TReCASE_sfitR$GAMMA)
c(TReCASE_sfitC$reg_par)
c(TReCASE_sfitR$betas, log(TReCASE_sfitR$phi))
c(TReCASE_sfitC$LL, TReCASE_sfitR$loglik)
c(TReCASE_sfitC$log_theta, log(TReCASE_sfitR$THETA))

trecaseC = RcppT_trecase(y, z, zz_AS, RHO, RHO_AS, X, tau1, tau2, lgy1, ni0, 
              ni, tauB, tau, lbc)
trecaseR = TReCASE_test(y,z, zz_AS, X, RHO, RHO_AS, tau1, tau2, ni0, ni, 
             tauB, tau)
exp(c(trecaseC$PAR))
c(trecaseR$KAPPA, trecaseR$ETA, trecaseR$GAMMA)
(c(trecaseC$reg_par))
c(trecaseR$betas, log(trecaseR$phi))
c(trecaseC$p_eta, trecaseC$p_gamma)
c(trecaseR$p.eta, trecaseR$p.gamma)
c(trecaseC$LL, trecaseC$LL_eta, trecaseC$LL_gamma)
c(trecaseR$loglik.full, trecaseR$loglik.eta, trecaseR$loglik.gamma)


# ----------------------------------------------------------------------
# TReC + ASE
# ----------------------------------------------------------------------
KEG_EaseGase =  log(c(KAPPA, ETA, GAMMA, 1.5, 1.8))
RcppT_TReC_ASE_LL_KEG(KEG_EaseGase, y, z, zz_AS, phi, RHO, RHO_AS, tau1, tau2, 
                      ni0, ni, log(theta), tauB, tau, lgy1, expXbeta, lbc, offsets, 
                      pis, mu) 
loglik.TReC_ASE_sep(KEG_EaseGase, y,z, zz_AS, X, betas, phi, RHO, 
                    RHO_AS, tau1, tau2, ni0, ni,  theta, tauB, tau)
RcppT_TReCASE_LL(y, phi, lgy1, mu, ni0, ni, log(theta), pis, lbc) 


RcppT_grad_NB(KEG_EaseGase[c(1:3)], H0, y, z, RHO, phi, tau1, tau2, expXbeta, offsets, 
              mu)
RcppT_grad_BB_KEG(KEG_EaseGase[c(1,4:5)], H0, zz_AS, RHO_AS, ni, ni0, 
                  log(theta), lbc, tauB, tau, pis) 
RcppT_TReC_ASE_grad_KEG(KEG_EaseGase, y, z, zz_AS, phi, RHO, RHO_AS, tau1, tau2, 
                         ni0, ni, log(theta), tauB, tau, lgy1, expXbeta, lbc, offsets, 
                         pis, mu) 
grad.TRe_CASE.keg_sep(KEG_EaseGase, y,z, zz_AS, X, betas, phi, RHO, 
                      RHO_AS, tau1, tau2, ni0, ni,  theta, tauB, tau)

RcppT_trec_ase_KEG_BFGS(KEG_EaseGase, y, z, zz_AS, RHO, RHO_AS, X, BETA, phi, tau1, 
                        tau2, lgy1, ni0, ni, log(theta), tauB, tau, lbc)
update_keg_trecase_sep(KEG_EaseGase, 0, y,z, zz_AS, X, betas, phi, RHO,
                       RHO_AS, tau1, tau2,ni0, ni, theta, tauB, tau)


trec_aseC = RcppT_trec_ase(c(0,0,0,0,0), y, z, zz_AS, RHO, RHO_AS, X, tau1, tau2, lgy1, 
               ni0, ni, tauB, tau, lbc, max_iter = 400)

trec_aseR = TReCASE_sep_sfit(KEG_EaseGase=rep(0,5), y,z,zz_AS, X, RHO, RHO_AS, 
                         tau1, tau2, ni0, ni, tauB, tau)
exp(c(trec_aseC$PAR))
c(trec_aseR$KEG)
(c(trec_aseC$reg_par))
c(trec_aseR$betas, log(trec_aseR$phi))
exp(trec_aseC$log_theta)
c(trec_aseR$THETA)
c(trec_aseC$LL)
c(trec_aseR$loglik)

# ----------------------------------------------------------------------
# Cis-Trans score test
# ----------------------------------------------------------------------
para = log(c(KAPPA, ETA, GAMMA))
source("./tumor_eQTL.R")
CisTrans_ScoreObs(para, y,z,zz_AS, X, BETA, phi,
                  RHO, RHO_AS, tau1, tau2, ni0, ni, tauB, tau, theta,
                  Power = F)

RcppT_CisTrans_ScoreObs(para, y, z, zz_AS, RHO, RHO_AS, X, BETA, phi, tau1, 
                  tau2, lgy1, ni0, ni, log(theta), tauB, tau, lbc) 

Expvec = c(0,0,0,0)
RcppT_ASE_ExpFunc(14, 0.5, 1/0.1, Expvec)
Expvec
ASE_ExpFunc(14, 0.5, 1/0.1)

CisTrans_Score(para, y,z,zz_AS, X, BETA, phi,
                  RHO, RHO_AS, tau1, tau2, ni0, ni, tauB, tau, theta,
                  Power = F)

RcppT_CisTrans_Score(para, y, z, zz_AS, RHO, RHO_AS, X, BETA, phi, tau1, 
                    tau2, lgy1, ni0, ni, log(theta), tauB, tau, lbc) 

