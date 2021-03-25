args = commandArgs(TRUE)
ETA = as.numeric(args[1])
GAMMA = as.numeric(args[2])
ETA
GAMMA

source("./tumor_eQTL.R")
library(doParallel)
library(MASS)
## ---------------------------------------------------------------------
## Simulation Function
## ---------------------------------------------------------------------

Data_generate <- function(rep, N, betas, X, phi, KAPPA, ETA, GAMMA, THETA, 
                          prob_phased=0.05, MAF=0.2)
{
  # simulated data 
  set.seed(2019+rep)
  RHO         = runif(N)
  geno_probs  = c((1-MAF)^2,MAF*(1-MAF),MAF*(1-MAF),MAF^2)
  dat         = data.frame(X) 
  dat$y       = NA
  prob_phased = 0.05
  dat$z       = sample(0:3,N,replace = TRUE,prob = geno_probs)
  prob_tau    = c(0.14, 0.68, 0.17, 0.01)
  tau1        = sample(0:3,N,replace = TRUE,prob_tau)
  tau2        = sample(0:3,N,replace = TRUE,prob_tau)
  tauB        = tau2
  tau         = tau1 + tau2
  
  offsets = compute_offset(dat$z, RHO, KAPPA, ETA, GAMMA, tau1, tau2)
  mus     = exp(X%*%betas+offsets)
  pis     = compute_pi(dat$z, RHO, KAPPA, ETA, GAMMA, tauB, tau)
  for(ii in seq(N)){
    dat$y[ii]  = rnbinom(n = 1,mu =mus[ii],size = 1/phi)
    dat$total_phased[ii] = rbinom(n = 1,size = dat$y[ii],prob = prob_phased)
    dat$y1[ii] = emdbook::rbetabinom(n = 1, prob = pis[ii],
                                     size = dat$total_phased[ii],theta = 1/THETA)
    dat$y2[ii] = dat$total_phased[ii] - dat$y1[ii]
  }
  
  z = dat$z 
  y = dat$y
  
  ind = which(dat$y1 + dat$y2 > 7)

  A = dat$y1[ind]
  D = A + dat$y2[ind]
  tauB = tau2[ind]
  tau = (tau1+tau2)[ind]
  z_AS = dat$z[ind]
  RHO_AS=RHO[ind]
  
  # trec = TReC_test(y,z, X, RHO, tau1, tau2, maxiter =200, tol = 1e-5)
  trecase = TReCASE_test(y,z,z_AS, X, RHO, RHO_AS, tau1, tau2, A, D, tauB, tau,
                        maxiter =200, tol = 1e-5)
  para = log(c(trecase$KAPPA, trecase$ETA, trecase$GAMMA))
  phi = trecase$phi
  THETA = trecase$THETA
  betas = trecase$betas

  score_obs = CisTrans_ScoreObs(para, y,z,z_AS, X, betas, phi,
                                RHO, RHO_AS, tau1, tau2, A, D, tauB, tau, THETA)
  
  #return(list(trec, trecase))
  return(c(unlist(trecase), unlist(score_obs)))
  # return(trec)
  
}


## ---------------------------------------------------------------------
## Simulation
## ---------------------------------------------------------------------
nsim = 500
N = 1000
set.seed(1)
betas=c(5, 0.15, -0.5, 0.2)
p <- length(betas)
X =data.matrix(cbind(1,matrix(runif(N*(p-1)),N,(p-1) )))

cl = makeCluster(4, outfile="")
registerDoParallel(cl)
nw <- getDoParWorkers() 
nw
time1 = Sys.time()
res_eg <- foreach(rep = 1:nsim, .combine = rbind) %dopar% {
  library(MASS)
  Data_generate(rep, N=1000, betas=betas,X=X,
                phi=0.2, KAPPA=1.5,ETA=ETA, GAMMA=GAMMA, THETA=0.2,
                prob_phased=0.05, MAF=0.2 )
}
time2 = Sys.time()
time2 - time1

mean(res_eg[,'p.score']<0.05)
save(ETA, GAMMA,res_eg, file = sprintf('trecase_ETA%s_GAMMA%s.Rdata', ETA, GAMMA))
#save(ETA, GAMMA,res_eg, file = sprintf('trec_ETA%s_GAMMA%s.Rdata', ETA, GAMMA))
q('no')
## ---------------------------------------------------------------------
## Summary and plot
## ---------------------------------------------------------------------

mean(unlist(res_eg[,'KAPPA']))
mean(unlist(res_eg[,'GAMMA']))
mean(unlist(res_eg[,'ETA']))
mean(unlist(res_eg[,'phi']))
mean(unlist(res_eg[,'THETA']))
mean(unlist(res_eg[,'p.eta']) <0.05)
mean(unlist(res_eg[,'p.gamma']) <0.05)
betas_est = sapply(res_eg[,'betas'], '[')
dim(betas_est)
apply(betas_est, 1, mean)

q(save='no')

E = 1
for(G in c(0.5,0.8,1,1.2,1.5))
  cat(sprintf("sbatch R CMD BATCH '--args %s %s' Simulated_data.R Simulated_data_E%s_G%s.Rout \n", E, G,E,G))

G = 1
for(E in c(0.5,0.8,1.2,1.5))
  cat(sprintf("sbatch R CMD BATCH '--args %s %s' Simulated_data.R Simulated_data_E%s_G%s.Rout \n", E, G,E,G))


library(MASS)
options(digits = 7)
load("/Users/lhuang/Dropbox/_tumor_eQTL/pTReCASE/data/calibrationData.RData")
head(calibrationData)
tau1 = tau2 = sample(1, 500, replace = T)
table(calibrationData$z)
calibrationData$z[calibrationData$z == 1] = 2
calibrationData$z[calibrationData$z == 0] = 1
table(calibrationData$z, useNA = 'ifany')

z=calibrationData$z;RHO=calibrationData$rhos;
X = data.matrix(calibrationData$Xs)
y = calibrationData$y

ind = which(calibrationData$y1 + calibrationData$y2 > 7)

A = calibrationData$y1[ind]
D = A + calibrationData$y2[ind]
tauB = tau1[ind]
tau = (tau1+tau2)[ind]
z_AS=calibrationData$z[ind]
A[which(z_AS ==2)] = (D - A)[which(z_AS ==2)]
RHO_AS=calibrationData$rhos[ind]

trec = TReC_test(y,z, X, RHO, tau1, tau2, maxiter =200, tol = 1e-5)
trecase = TReCASE_test(y,z, z_AS, X, RHO, RHO_AS, tau1, tau2, A, D, tauB, tau, maxiter =200, tol = 1e-5)
ase = ASE_test(z_AS, RHO_AS, A, D, THETA, tauB, tau)
CisTrans_lrt(trecase, trec, ase)
  
para = log(c(1.175004,	0.9976487, 1.63688))
phi =  0.2150183
betas =3.956203
THETA = 0.182807
H0 =0
if(H0 == 0){
  KAPPA = exp(para[1])
  ETA = exp(para[2])
  GAMMA = exp(para[3])
} else if(H0 == 1){
  
  ETA = 1.0;
  KAPPA = exp(para[1])
  GAMMA = exp(para[2])
  
} else if(H0 == 2){
  GAMMA = 1.0
  
  KAPPA = exp(para[1])
  ETA = exp(para[2])
}	

loglikNB(para, H0, y,z, X, betas, phi, RHO, tau1, tau2)
Grad_NB(para, H0, y,z, X, betas, phi, RHO, tau1, tau2)
loglikBB(para,H0, A,D,THETA, z_AS, RHO_AS, tauB, tau)
Grad_BB_keg(para,H0, A,D,THETA, z_AS, RHO_AS, tauB, tau)
loglikNB(para, H0, y,z, X, betas, phi, RHO, tau1, tau2) +
  loglikBB(para,H0, A,D,THETA, z_AS, RHO_AS, tauB, tau)

Grad_NB(para, H0, y,z, X, betas, phi, RHO, tau1, tau2) +
  Grad_BB_keg(para,H0, A,D,THETA, z_AS, RHO_AS, tauB, tau)


KAPPA = exp(para[1])
ETA = exp(para[2])
GAMMA = exp(para[3])
pis <- compute_pi(z_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau)
loglikBB_THETA(THETA, A,D, pis)
lASE.dTHETA(THETA, A,D, pis)

CisTrans_ScoreObs(para, y,z,z_AS, X, betas, phi,
                RHO, RHO_AS, tau1, tau2, A, D, tauB, tau, THETA)

CisTrans_Score(para, y,z,z_AS, X, betas, phi,
               RHO, RHO_AS, tau1, tau2, A, D, tauB, tau, THETA)




offsets <- compute_offset(z, RHO, KAPPA, ETA, GAMMA, tau1, tau2)
exp((X%*%betas) + offsets)
