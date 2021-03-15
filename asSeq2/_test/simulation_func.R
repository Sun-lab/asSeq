# setwd("/fh/fast/sun_w/licai/_tumor_eQTL/_simulation_asSeq2/")
library(MASS)

args = commandArgs(TRUE)
# args = c("1.5", "1")
GAMMA = as.numeric(args[1])
ETA = as.numeric(args[2])
GAMMA
ETA
KAPPA = 1.5
useASE = 1

# source("simulation_function_tumor_eqtl.R")

set.seed(2020)
NN = 500
numSNPs = 1

BETA    = c(5,0.15,-0.5,0.25)
PHI     = 0.1
MAF     = 0.2
XX      = data.frame(cbind(1,matrix(runif(NN*3),NN,3)), stringsAsFactors = F) 
names(XX) = paste0("x",1:4)
RHO     = runif(NN, 0.3, 0.9)

PSI = THETA = 100 # over-dispersion = 0.01

prob_phased = 0.1
corrSNPs    = 0
eQTL_index  = 1
prob_tau    = c(0.14, 0.68, 0.17, 0.01)

nSims  = 22

## ---------------------------------------------------------------------
## Simulation
## ---------------------------------------------------------------------

compute_offset2 <- function(z, RHO, KAPPA, ETA, GAMMA, tau1, tau2){
  n = length(z)
  offsets = vector('numeric', n)
  indAA = which(z == 0)
  offsets[indAA] = log(2*(1-RHO[indAA]) +
                         (tau1[indAA]+tau2[indAA])*RHO[indAA]*KAPPA)
  
  indAB = which(z == 1)
  offsets[indAB] = log((1-RHO[indAB]) + tau1[indAB]*RHO[indAB]*KAPPA +
                         (1-RHO[indAB])*ETA + tau2[indAB]*RHO[indAB]*KAPPA*GAMMA)
  
  indBA = which(z == 2)
  offsets[indBA] = log((1-RHO[indBA]) + tau2[indBA]*RHO[indBA]*KAPPA +
                         (1-RHO[indBA])*ETA + tau1[indBA]*RHO[indBA]*KAPPA*GAMMA)
  
  indBB = which(z == 3)
  offsets[indBB] = log(2*(1-RHO[indBB])*ETA +
                         (tau1[indBB]+tau2[indBB])*RHO[indBB]*KAPPA*GAMMA)
  
  return(offsets)
}

compute_pi <- function(z_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau){
  pis = vector('numeric', length(z_AS))
  indhomo = which(z_AS %in% c(0,3))
  pis[indhomo] = (RHO_AS[indhomo]*tauB[indhomo]*KAPPA + 1 -
                    RHO_AS[indhomo])/(RHO_AS[indhomo]*tau[indhomo]*KAPPA
                                      + 2*(1 - RHO_AS[indhomo]))
  indhet = which(z_AS %in% c(1,2))
  tmp1 = RHO_AS[indhet]*tauB[indhet]*GAMMA*KAPPA + (1 - RHO_AS[indhet])*ETA
  pis[indhet] = tmp1/(RHO_AS[indhet]*(tau[indhet]-tauB[indhet])*KAPPA +
                        (1 - RHO_AS[indhet]) + tmp1)
  return(pis)
}


Data_generate <- function(rep, N, betas, XX, phi, KAPPA, ETA, GAMMA, THETA, RHO,
                          prob_phased=0.1, MAF=0.2)
{
  # simulated data 
  set.seed(2020+rep)
  # RHO         = runif(N, 0.3, 0.9)
  geno_probs  = c((1-MAF)^2,MAF*(1-MAF),MAF*(1-MAF),MAF^2)
  dat         = data.frame(RHO)
  ZZ          = sample(0:3,N,replace = TRUE,prob = geno_probs)
  prob_tau    = c(0.14, 0.68, 0.17, 0.01)
  tau1        = sample(0:3,N,replace = TRUE,prob_tau)
  tau2        = sample(0:3,N,replace = TRUE,prob_tau)
  tauB        = tau2
  ind0 = which(ZZ %in% c(2,3))
  tauB[ind0]  = tau1[ind0] 
  tau         = tau1 + tau2
  # cbind(tauB, tau1, tau2, dat$z)[1:30,]
  
  
  offsets = compute_offset2(ZZ, RHO, KAPPA, ETA, GAMMA, tau1, tau2)
  mus     = exp(as.matrix(XX)%*%betas+offsets)
  pis     = compute_pi(ZZ, RHO, KAPPA, ETA, GAMMA, tauB, tau)
  for(ii in seq(N)){
    dat$total[ii]  = rnbinom(n = 1,mu =mus[ii],size = 1/phi)
    dat$total_phased[ii] = rbinom(n = 1,size = dat$total[ii],prob = prob_phased)
    dat$hapB[ii] = emdbook::rbetabinom(n = 1, prob = pis[ii],
                                       size = dat$total_phased[ii],theta = 1/THETA)
    dat$hapA[ii] = dat$total_phased[ii] - dat$hapB[ii]
  }
  
  # final data
  dat$hap1 = dat$hapA
  ind1 = which(ZZ %in% c(2,3))
  dat$hap1[ind1] = dat$total_phased[ind1] - dat$hapA[ind1]
  dat$hap2 = dat$total_phased -  dat$hap1
  
  dat$tau1 = tau1
  dat$tau2 = tau2
  
  return(list(dat,ZZ))
  
}

func1 <- function(repi){
  set.seed(2020 + repi)
  gene_name = paste0("G", repi) 
  # gen_gene_RC_SNPs_cnv(gene_name,XX,numSNPs,RHO, BETA,PHI,ETA,KAPPA,GAMMA,
  #                      PSI,MAF,tau1,tau2,prob_phased,corrSNPs,eQTL_index)
  Data_generate(repi, NN, BETA, XX, PHI, KAPPA, ETA, GAMMA, THETA, RHO,
                prob_phased=0.05, MAF=0.2)
}

dataSim_raw = lapply(1:nSims, func1) 
dataGene    = sapply(dataSim_raw, "[", 1)

cnm = colnames(dataGene[[1]])
cnm

Y  = data.matrix(sapply(dataGene, "[[", which(cnm == "total")))
Y1 = data.matrix(sapply(dataGene, "[[", which(cnm == "hap1")))
Y2 = data.matrix(sapply(dataGene, "[[", which(cnm == "hap2")))
Z  = data.matrix(sapply(dataSim_raw, "[[", 2))
XX = data.matrix(XX)
CNV1 = data.matrix(sapply(dataGene, "[[", which(cnm == "tau1")))
CNV2 = data.matrix(sapply(dataGene, "[[", which(cnm == "tau2")))
# RHO = data.matrix(sapply(dataGene, "[[", which(cnm == "RHO")))

geneloc = data.frame(gene = paste0('gene', 1:nSims), chr = 1:nSims, start = c(1),
                     end = c(1e3)+1000, stringsAsFactors = F)
head(geneloc)

SNPloc = data.frame(snp = paste0("SNP", 1:nSims), chr = 1:nSims,
                    pos = 1:nSims, stringsAsFactors = F)
head(SNPloc)


res_all = NULL
library(asSeq2, lib = "/home/lhuang2/R/x86_64-pc-linux-gnu-library/4.0")

file_trecase_all = sprintf("simu1_trecase_GAMMA%s_ETA%s.txt",GAMMA,ETA)
file_trec_all    = sprintf("simu1_trec_GAMMA%s_ETA%s.txt",GAMMA,ETA)


if(useASE == 1){
  file_name = file_trecase_all
}else{
  file_name = file_trec_all
}

time1 = Sys.time()
trecaseT(Y, Y1, Y2, Z, XX, RHO, CNV1, CNV2,
         SNPloc, geneloc, GeneSnpList = split(1:nSims, seq(nSims)),
         file_trec = file_trec_all, file_trecase = file_name,
         useLRT = T, transTestP = 0, cis_window = 100000, useASE = useASE,
         min_ASE_total = 8L, min_nASE = 5L, min_nASE_het = 5L, eps = 5e-5,
         max_iter = 200L, show = FALSE)

time2 = Sys.time()
time2 - time1

res_all = read.table(file_name, sep = "\t", header = T, stringsAsFactors = FALSE)


summary(res_all)
mean(res_all$TReCASE_pGamma<0.05)
mean(res_all$TReCASE_pEta<0.05)
mean(res_all$TReC_pGamma<0.05)
mean(res_all$TReC_pEta<0.05)
mean(res_all$CisTrans_Pvalue<0.05, na.rm = T)

table(res_all$Converge)

q("no")


Eta = 1
Gamma = 1.5
cat(sprintf("sbatch R CMD BATCH '--args %s %s' simulation_func.R simulation_func_GAMMA%s_ETA%s.Rout \n", Gamma,Eta, Gamma, Eta))

