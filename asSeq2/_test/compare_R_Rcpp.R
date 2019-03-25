setwd('/fh/fast/sun_w/licai/_tumor_eQTL/GitHub/asSeq/asSeq2/_test')
library(Rcpp)
Rcpp::sourceCpp('../asSeq2/src/asSeq2.cpp')
source("trecR.R")
source("aseR.R")
source("trecaseR.R")

# load('_test_func.Rdata')
load('../_test/_test_func.Rdata')
options(digits=10)

N
dat[1:5,]
bxj
lgy1 = dat$LGX1
zz =zz2 = ZZ[,1]
table(zz)
zz[which(zz == 2)] = 1
zz[which(zz == 3)] = 2
table(zz)
X = data.matrix(XX)
betas = BETA
phi = PHI
y = dat$total
b0 = 0

theta = PSI
Pi = exp(bxj)/(exp(bxj) + 1)
ind = which(dat$total_phased > 7)
zz_AS = zz[ind]
ni = dat$total_phased[ind]
ni0 = dat$hapB[ind]
lbc = dat$LBC[ind]
zeta = zz_AS==1

offsets = rep(0, N)
mu = rep(0, N)

#-------------------------------------------------------------
# TReC model
#-------------------------------------------------------------

## negtive binomial
Rcpp_logLTReC(bxj, y, X, zz, betas, phi, T, lgy1, mu)
logLTReC(bxj, y, zz, exp(X%*%betas), b0, phi, 'negbin')

Rcpp_trec_grad_bxj(bxj, y, zz, mu, phi, T)
grad.bxj.trec(bxj, y, zz, exp(X%*%betas), b0, phi, 'negbin')

Rcpp_reg_LL(y, X, offsets, PARAMS = c(betas, log(phi)), fam_nb = T, lgy1, mu)
Rcpp_reg_grad(y, X, mu, PARAMS = c(betas, log(phi)), fam_nb = T)
Rcpp_reg_Hess(y, X, mu, PARAMS = c(betas, log(phi)), fam_nb = T)
Rcpp_reg(y, X, offsets, rep(0.1, 5),fam_nb = T, lgy1, max_iter = 4000L,
         eps = 1e-5, show = TRUE)
reg0 = Rcpp_reg_BFGS(y, X, offsets, rep(0, 5),fam_nb = T, lgy1, max_iter = 4000L,
              eps = 1e-5, show = TRUE)
reg0
library(MASS)

g1  = glm.nb(y ~ X+offset(offsets))
BETA = summary(g1)$coef[,1]
BETA
phi = 1/g1$theta
log(phi)
g1$twologlik/2


Cop = Rcpp_trec_bxj_BFGS(0, y, X, zz, BETA, phi, T, lgy1, max_iter = 4000L,
              eps = 1e-5, show = TRUE)
Rop = optim(0, fn=logLTReC, gr=grad.bxj.trec, y=y, x=zz, mu=exp(X%*%BETA), b0=b0, phi=phi,
      fam='negbin', method="BFGS", control=list(fnscale=-1.0, trace=0))
Cop
Rop

library(microbenchmark)
microbenchmark(Rcpp_trec_bxj_BFGS(0, y, X, zz, BETA, phi, T, lgy1, max_iter = 4000L,
                               eps = 1e-5, show = F),
               optim(0, fn=logLTReC, gr=grad.bxj.trec, y=y, x=zz, mu=exp(X%*%BETA), b0=b0, phi=phi,
                     fam='negbin', method="BFGS", control=list(fnscale=-1.0, trace=0)))

c(Cop$PAR, Rop$par)

Rcpp_logLTReC(Cop$PAR, y, X, zz, BETA, phi, T, lgy1, mu)
Rcpp_logLTReC(Rop$par, y, X, zz, BETA, phi, T, lgy1, mu)
Rcpp_logLTReC((Rop$par+Cop$PAR)/2, y, X, zz, BETA, phi, T, lgy1, mu)

lls = sapply(seq(Rop$par, Cop$PAR, by = 1e-6), Rcpp_logLTReC, y, X, zz,
             BETA, phi, T, lgy1, mu)
lls
grs = sapply(seq(Rop$par, Cop$PAR, by = 1e-6), Rcpp_trec_grad_bxj, y,
             zz, mu, phi, T)
grs

time1 = Sys.time()
Ctrec3 = Rcpp_trec2(y, X, zz, T, lgy1, max_iter = 4000L, eps = 1e-07, show = T)
Ctrec = Rcpp_trec2(y, X, zz, T, lgy1,ini_bxj = 0.0, LL_null=reg0$LL,
                   ini_reg_par = reg0$PAR,
                   max_iter = 4000L, eps = 1e-07, show = T)
time2 = Sys.time()
Rtrec = trecR(y, X, zz, 'negbin', yfit = T)
time3 = Sys.time()
c(time2-time1, time3-time2)

Ctrec
Rtrec
Ctrec$bxj
Rtrec$b
Ctrec$bxj - Rtrec$b

Ctrec$reg_par[1:4]
Rtrec$betas

exp(Ctrec$reg_par[5])
Rtrec$phi

# q('no')




## poisson
Rcpp_logLTReC(bxj, y, X, zz, betas, phi, F, lgy1, mu)
logLTReC(bxj, y, zz, exp(X%*%betas), b0, phi, 'poisson')

Rcpp_trec_grad_bxj(bxj, y, zz, mu, 0, F)
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

Rcpp_trec_bxj_BFGS(0, y, X,  zz, BETA, phi, F, lgy1)
optim(0, fn=logLTReC, gr=grad.bxj.trec, y=y, x=zz, mu=exp(X%*%BETA), b0=b0, phi=phi,
      fam='poisson', method="BFGS", control=list(fnscale=-1.0, trace=0))

time1 = Sys.time()
Rcpp_trec(y, X, zz, F, lgy1, max_iter = 4000L, eps = 1e-04,
          show = T)
time2 = Sys.time()
trecR(y, X, zz, 'poisson')
time3 = Sys.time()
c(time2-time1, time3-time2)

#-------------------------------------------------------------
# ASE model
#-------------------------------------------------------------
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
Rcpp_ase_grad_H0(ni, ni0, Pi, log(theta), zeta)
gradLogH1(c(theta, Pi), ni0, ni, zeta)
Rcpp_ase_grad_Pi(ni, ni0, Pi, log(theta), zeta)

Rcpp_ase_grad_H0(ni, ni0, 0.5, log(theta), zeta)
gradLogH0(c(theta, 0.5), ni0, ni)

# Rcpp_ase_hess(ni, ni0, 0.51602267, 0.06677156, zeta)

par0 = 0.1
Rcpp_ase_BFGS(ni, ni0, zeta, log(par0), lbc, max_iter = 4000L,
              eps = 1e-5, show = TRUE)
Rcpp_ase_theta_BFGS(ni, ni0, zeta, 0.5, log(par0), lbc, max_iter = 4000L,
              eps = 1e-5, show = TRUE)
op0  = optim(par0, logH0, gr=gradLogH0, nA=ni0, nTotal=ni,
             method="L-BFGS-B", lower=1e-16, upper=Inf,
             control=list(fnscale=-1))
op0

Rcpp_ase_BFGS(ni, ni0, zeta, c(0.5, -2.705298513+0.2), lbc, max_iter = 4000L,
              eps = 5e-5, show = TRUE)
optim(c(op0$par, 0.5), logH1, gr=gradLogH1, nA=ni0, nTotal=ni, zeta=zeta,
      method="L-BFGS-B", lower=c(0,0) + 1e-16, upper=c(Inf, 1-1e-16),
      control=list(fnscale=-1),  hessian = T)

microbenchmark(
  Rcpp_ase_theta_BFGS(ni, ni0, zeta, 0.5, -1, lbc, max_iter = 4000L,
                    eps = 5e-04, show = FALSE),
  Rcpp_ase_theta_BFGS2(ni, ni0, zeta, 0.5, -1, lbc, max_iter = 4000L,
                    eps = 5e-04, show = FALSE)

)
optimize(loglikTheta, interval=c(0, 1000), pi=0.5, nA=ni0,
         nTotal=ni, zeta=zeta, maximum=TRUE)

time1 = Sys.time()
Rcpp_ase(ni, ni0, zeta, lbc, max_iter = 4000L,
         eps = 5e-05, show = FALSE)
time2 = Sys.time()
aseR(ni0, ni, zeta)
time3 = Sys.time()
c(time2-time1, time3-time2)

#-------------------------------------------------------------
# TReCASE
#-------------------------------------------------------------
source("trecaseR.R")
mu = rep(0, length(y))
Rcpp_trecase_LL(bxj, y, X, zz, BETA, phi, 1, lgy1, mu, ni, ni0, log(theta),
                lbc, zeta)
loglikJoin(bxj, y, zz, ni0, ni, zeta, exp(X%*%BETA), 0, phi, theta)

Rcpp_trecase_grad_bxj(bxj, y, X, zz, BETA, phi, 1, lgy1, mu, ni, ni0, log(theta),
                      lbc, zeta)
grad.bxj(bxj, y, zz, ni0, ni, zeta, exp(X%*%BETA), 0, phi, theta)

time1 = Sys.time()
Rcpp_trecase_BFGS(0, y, X, zz, BETA, phi, T, lgy1, ni, ni0, log(theta),
                  lbc, zeta, max_iter = 4000L, eps = 1e-07, show = TRUE)
time2 = Sys.time()
optim(par=0, fn=loglikJoin, gr=grad.bxj, y=y, x=zz, nA=ni0, nTotal=ni,
           zeta=zeta, mu=exp(X%*%BETA), b0=0, phi=phi, theta=theta, method = "BFGS",
           control=list(fnscale=-1.0))
time3 = Sys.time()
c(time2-time1, time3-time2)

time1 = Sys.time()
Ctrecase = Rcpp_trecase(y, X, zz, 1, lgy1, ni, ni0, zeta, lbc, max_iter = 4000L,
             eps = 5e-05, show = T)
time2 = Sys.time()
Rtrecase = trecaseR(y, ni0, ni, X, zz, plotIt=FALSE, traceIt=FALSE)
time3 = Sys.time()
c(time2-time1, time3-time2)

c(Ctrecase$bxj, Rtrecase$b)
c(exp(Ctrecase$lg_theta), Rtrecase$theta)
c(exp(Ctrecase$reg_par[5]), Rtrecase$phi)
c((Ctrecase$lrt), Rtrecase$lrt)

#-------------------------------------------------------------
# Rcpp wrapper
#-------------------------------------------------------------
geneloc = data.frame(gene = paste0('gene', 1:2), chr = 1:2, start = c(1, 1e8),
                     end = c(1, 1e8)+1000, stringsAsFactors = F)
geneloc

SNPloc = data.frame(snp = paste0("SNP", 1:100), chr = c(rep(1, 40), rep(2, 60)),
                    pos = c(1:30*100, 31:100*100 + 1e8), stringsAsFactors = F)
head(SNPloc)

source("../Rcpp/asSeq_interface.R")

load('../_test/_test_func.Rdata')
Y = cbind(dat$total, dat$total)
y1= dat$hapA
y1[ZZ[,1]==2] = dat$hapB[ZZ[,1]==2]
y2 = dat$total_phased - y1
Y1 = cbind(y1,y1)
Y2 = cbind(y2,y2)

## R version
time2 = Sys.time()
res2 = trecase2(Y, Y1, Y2, ZZ,
                XX, SNPloc, geneloc, 1,
              cis_window = 1e5, useASE = F,
              min_ASE_total=8, min_nASE=10, eps=1e-4, show = F)
time3 = Sys.time()
c( time3-time2)

length(res)
res2[[1]][[1]][1:5,]
res2[[2]]

str(res)


Y = data.matrix(Y)
Y1 =data.matrix(Y1)
Y2 =data.matrix(Y2)
ZZ =data.matrix(ZZ)
XX = data.matrix(XX)

# new version
time1 = Sys.time()
Rcpp_trecase_mtest(Y, Y1, Y2, ZZ, XX, SNPloc[,3], SNPloc[,2], T,
                   geneloc[,3], geneloc[,4], geneloc[,2],"trec.txt",
                   "trecase.txt",
                    eps = 5e-05, useASE = 1)
time2 = Sys.time()
time2-time1

# asSeq package
library(asSeq)
Z = ZZ
Z[ZZ > 1] = Z[ZZ>1] + 1
time1 = Sys.time()
trecase(Y, Y1, Y2, XX[,-1], Z, "output_trecase",1, min.AS.reads=8, min.AS.sample=10,
                eChr=geneloc[,2], ePos=geneloc[,4], mChr=SNPloc[,2], mPos=SNPloc[,3])
time2 = Sys.time()
time2-time1

res= read.table("trec.txt", sep = "\t", header = T)
res2= read.table("output_eqtl.txt", sep = "\t", header = T)
dim(res)
dim(res2)

res[,"TReC_Pvalue"]
res2[,"TReC_Pvalue"]
summary(res[,"TReC_Pvalue"]- res2[,"TReC_Pvalue"])
plot(res[,"TReC_Pvalue"], res2[,"TReC_Pvalue"])
abline(0,1)
dim(Rres)


time1 = Sys.time()
Rcpp_trecase_mtest(Y, Y1, Y2, ZZ, XX, as.numeric(SNPloc[,3]), as.numeric(SNPloc[,2]), T,
                   as.numeric(geneloc[,3]),as.numeric(geneloc[,4]),
                   as.numeric(geneloc[,2]), "trec2.txt",
                   "trecase2.txt",
                   eps = 5e-05, useASE = 1,
                   show =F)
time2 = Sys.time()
time2-time1

library(microbenchmark)
microbenchmark(
  trec(Y, XX[,-1], ZZ3, "output_trec", 1, local.distance = 1e+05,
       eChr=geneloc[,2], ePos=geneloc[,4], mChr=SNPloc[,2], mPos=SNPloc[,3]),
  Rcpp_trecase_mtest(Y, Y1, Y2, ZZ, XX[, 1, drop = F], SNPloc[,3], SNPloc[,2], T,
                     geneloc[,3], geneloc[,4], geneloc[,2],"trec.txt",
                     "trecase.txt",
                     eps = 5e-05, useASE = 0),
  times = 30L
)
#-------------------------------------------------------------
# Rcpp wrapper
#-------------------------------------------------------------
geneloc = data.frame(gene = paste0('gene', 1:2), chr = 1:2, start = c(1, 1e8),
                     end = c(1, 1e8)+1000, stringsAsFactors = F)
geneloc

SNPloc = data.frame(snp = paste0("SNP", 1:100), chr = c(rep(1, 40), rep(2, 60)),
                    pos = c(1:30*100, 31:100*100 + 1e8), stringsAsFactors = F)
head(SNPloc)

load('./_test_func.Rdata')
Y = cbind(dat$total, dat$total)
y1= dat$hapA
y1[ZZ[,1]==2] = dat$hapB[ZZ[,1]==2]
y2 = dat$total_phased - y1
Y1 = cbind(y1,y1)
Y2 = cbind(y2,y2)



library(asSeq2)
time1 = Sys.time()
trecase(Y, Y1, Y2, ZZ, XX, SNPloc, geneloc, fam_nb = T,eps=1e-8, show = F)
time2 = Sys.time()
time2-time1

time1 = Sys.time()
ase(Y1, Y2, ZZ, geneloc,SNPloc, "ase.txt")
time2 = Sys.time()
time2-time1

ZZ3 = ZZ
ZZ3[ZZ==2] = 1
ZZ3[ZZ==3] = 2
time1 = Sys.time()
trec(Y, XX[,-1], ZZ3, "output_trec", 1, local.distance = 1e+05,
     eChr=geneloc[,2], ePos=geneloc[,4], mChr=SNPloc[,2], mPos=SNPloc[,3])
time2 = Sys.time()
time2-time1

TrecFastTest:::trec_fast(Y, XX[,-1], ZZ3, "output_trec", 1, local.distance = 1e+05,
                         eChr=geneloc[,2], ePos=geneloc[,4], mChr=SNPloc[,2], mPos=SNPloc[,3])