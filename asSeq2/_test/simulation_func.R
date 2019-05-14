
smart_df = function(...){
  data.frame(...,stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------
## Simulation Function
## ---------------------------------------------------------------------

gen_hap_v2<-function(SNPcor, MAF, ZZ){
  ## SNPcor: correlation between adjecent SNPs
  ## MAF: minor allele frequency 
  ## ZZ: adjecent haplotype
  library(rootSolve)
  N = length(ZZ)
  fn <- function(x, SNPcor, MAF){
    EG1 = x[1]*MAF +x[2]*(1-MAF)
    c(F1 = (EG1 - MAF ) ,
      F2 = (x[1]*MAF - EG1*MAF)/sqrt((EG1*(1-EG1)*MAF*(1-MAF))) - SNPcor)
  }
  delta.p= multiroot(fn, c(-1e-3,1e-3), SNPcor=SNPcor, MAF=MAF)$root
  if(any(delta.p < 0) | any(delta.p > 1)){
    stop('delta.p out of bound')
  }
  hap = vector('numeric', N)
  ind1 = which(ZZ==1)
  hap[ind1] = rbinom(length(ind1), 1, prob = delta.p[1])
  hap[-ind1] = rbinom(N- length(ind1), 1, prob = delta.p[2])
  # print(cor(hap, ZZ))
  return(hap)
}


gen_gene_RC_SNPs = function(gene_name,XX,numSNPs,BETA,PHI,bxj,
                            PSI,MAF,prob_phased,corrSNPs,eQTL_index){
  if(FALSE){
    set.seed(1)
    gene_name = "G1"
    NN = 300
    XX = smart_df(cbind(1,matrix(runif(NN*3),NN,3))); names(XX) = paste0("x",1:4)
    numSNPs = 100
    BETA = c(5,0.15,-0.5,0.25)
    PHI = 0.1
    bxj = log(1.5)
    PSI = 0.05
    MAF = 0.2
    prob_phased = 0.05
    corrSNPs = 0
    eQTL_index = 1
  }
  
  N = nrow(XX)
  
  # Generate SNPs
  # MAF = 0.2; N = 500; numSNPs = 10; corrSNPs = 0
  ZZ1 = ZZ2 = matrix(NA, nrow = N, ncol = numSNPs )
  ZZ1[,1] = rbinom(N, 1, prob = MAF)
  ZZ2[,1] = rbinom(N, 1, prob = MAF)
  if(corrSNPs == 0){
    for(coli in 2:numSNPs){
      ZZ1[,coli] = rbinom(N, 1, prob = MAF)
      ZZ2[,coli] = rbinom(N, 1, prob = MAF)
    }
  }else{
    for(coli in 2:numSNPs){
      ZZ1[,coli] = gen_hap_v2(corrSNPs, MAF, ZZ1[, coli-1]) 
      ZZ2[,coli] = gen_hap_v2(corrSNPs, MAF, ZZ2[, coli-1]) 
  }
  
  }

  ZZ = ZZ1+2*ZZ2

  colnames(ZZ) = paste0("Z",seq(numSNPs))
  # round(cor(ZZ),3); ZZ[1:5,]
  
  dat = smart_df(total = rep(NA,N))
  dat$total_phased = NA
  dat$hapA = 0
  dat$hapB = 0
  dat$LBC = 1
  dat$phased = NA
  dat$Z = ZZ[,paste0("Z",eQTL_index)]
  dat$pp = NA
  # First snp corresponds to cts eQTL, other snps are not associated
  
  vec_log_mu_AA = as.numeric(as.matrix(XX) %*% BETA)
  for(ii in seq(N)){
    # ii = 1
    if( dat$Z[ii] == 0 ){ # AA
      log_mu = vec_log_mu_AA[ii]
    } else if( dat$Z[ii] %in% c(1,2) ){ # AB or BA
      log_mu = vec_log_mu_AA[ii] + log((1+exp(bxj))/2)
    } else if( dat$Z[ii] == 3 ){ # BB
      log_mu = vec_log_mu_AA[ii] + bxj
    }
    dat$total[ii] = rnbinom(n = 1,mu = exp(log_mu),size = 1/PHI)
    dat$total_phased[ii] = rbinom(n = 1,size = dat$total[ii],prob = prob_phased)
    dat$phased[ii] = ifelse(dat$total_phased[ii] > 0,1,0)
    if( dat$phased[ii] == 1 ){
      if( dat$Z[ii] == 1 ){ # AB
        dat$pp[ii] = exp(bxj) / (1 + exp(bxj))
      } else if( dat$Z[ii] == 2 ){ # BA
        dat$pp[ii] = exp(bxj) / (1 + exp(bxj))
      } else { # AA or BB
        dat$pp[ii] = 0.5
      }
    }
  }
  for(Z in c(0,1,2,3)){
    # Z = 0
    tmp_index = which(dat$phased == 1 & dat$Z == Z)
    len_index = length(tmp_index)
    if(len_index > 0){
      if(Z == 1){
        dat$hapB[tmp_index] = emdbook::rbetabinom(n = len_index, 
                                                  prob = dat$pp[tmp_index],
                                                  size = dat$total_phased[tmp_index],
                                                  theta = 1/PSI)
        dat$hapA[tmp_index] = dat$total_phased[tmp_index] - dat$hapB[tmp_index]
      } else if(Z == 2){
        dat$hapA[tmp_index] = emdbook::rbetabinom(n = len_index,
                                                  prob = dat$pp[tmp_index],
                                                  size = dat$total_phased[tmp_index],
                                                  theta = 1/PSI)
        dat$hapB[tmp_index] = dat$total_phased[tmp_index] - dat$hapA[tmp_index]
      } else {
        dat$hapA[tmp_index] = emdbook::rbetabinom(n = len_index,
                                                  prob = dat$pp[tmp_index],
                                                  size = dat$total_phased[tmp_index],
                                                  theta = 1/PSI)
        dat$hapB[tmp_index] = dat$total_phased[tmp_index] - dat$hapA[tmp_index]
      }
    }
  }
  
  tmp_index = which(dat$phased == 1)
  dat$LBC[tmp_index] = as.numeric(apply(dat[tmp_index,c("total_phased","hapB")],1,function(x)
    lchoose(x[1],x[2])))
  dat$LGX1 = sapply(dat$total,function(x) lgamma(x+1))
  dat = dat[,which(!(names(dat) %in% c("pp","Z")))]
  
  list(dat,ZZ)
}

## ---------------------------------------------------------------------
## Simulation Function (Tumor eQTL)
## ---------------------------------------------------------------------

compute_offset <- function(z, RHO, KAPPA, ETA, GAMMA, tau1, tau2){
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

gen_gene_RC_SNPs_cnv = function(gene_name,XX,numSNPs,RHO, BETA,PHI,ETA,KAPPA,GAMMA,
                                PSI,MAF,tau1,tau2,prob_phased,corrSNPs,eQTL_index){
  if(FALSE){
    set.seed(1)
    gene_name = "G1"
    NN = 300
    XX = smart_df(cbind(1,matrix(runif(NN*3),NN,3))); names(XX) = paste0("x",1:4)
    numSNPs = 100
    RHO  = runif(NN)
    BETA = c(5,0.15,-0.5,0.25)
    PHI = 0.1
    KAPPA = 1.5
    ETA = 0.8
    GAMMA = 1.2
    PSI = 0.05
    MAF = 0.2
    prob_phased = 0.05
    corrSNPs = 0
    eQTL_index = 1
    prob_tau    = c(0.14, 0.68, 0.17, 0.01)
    tau1        = sample(0:3,NN,replace = TRUE, prob_tau)
    tau2        = sample(0:3,NN,replace = TRUE, prob_tau)
  }
  
  
  N = nrow(XX)
  
  # Generate SNPs
  # MAF = 0.2; N = 500; numSNPs = 10; corrSNPs = 0
  ZZ1 = ZZ2 = matrix(NA, nrow = N, ncol = numSNPs )
  ZZ1[,1] = rbinom(N, 1, prob = MAF)
  ZZ2[,1] = rbinom(N, 1, prob = MAF)
  if(corrSNPs == 0){
    for(coli in 2:numSNPs){
      ZZ1[,coli] = rbinom(N, 1, prob = MAF)
      ZZ2[,coli] = rbinom(N, 1, prob = MAF)
    }
  }else{
    for(coli in 2:numSNPs){
      ZZ1[,coli] = gen_hap_v2(corrSNPs, MAF, ZZ1[, coli-1]) 
      ZZ2[,coli] = gen_hap_v2(corrSNPs, MAF, ZZ2[, coli-1]) 
    }
    
  }
  
  ZZ = ZZ1+2*ZZ2
  
  colnames(ZZ) = paste0("Z",seq(numSNPs))
  # round(cor(ZZ),3); ZZ[1:5,]
  
  dat = smart_df(total = rep(NA,N))
  dat$total_phased = NA
  dat$hapA = 0
  dat$hapB = 0
  dat$LBC = 1
  dat$phased = NA
  dat$Z = ZZ[,paste0("Z",eQTL_index)]
  dat$tau = tau1 + tau2
  dat$tau1 = tau1
  dat$tau2 = tau2
  dat$tauB = tau2
  dat$tauB[which(dat$Z==2)] = tau1[which(dat$Z==2)]
  # First snp corresponds to cts eQTL, other snps are not associated
  offsets = compute_offset(dat$Z, RHO, KAPPA, ETA, GAMMA, tau1, tau2)
  mus     = exp(as.matrix(XX)%*%BETA+offsets)
  dat$pis =  compute_pi(dat$Z, RHO, KAPPA, ETA, GAMMA, dat$tauB, dat$tau)

  for(ii in seq(N)){

    dat$total[ii] = rnbinom(n = 1,mu = mus[ii], size = 1/PHI)
    dat$total_phased[ii] = rbinom(n = 1,size = dat$total[ii],prob = prob_phased)
    dat$phased[ii] = ifelse(dat$total_phased[ii] > 0,1,0)
    dat$hapB[ii] = emdbook::rbetabinom(n = 1, 
                                   prob = dat$pis[ii],
                                   size = dat$total_phased[ii],
                                   theta = 1/PSI)
    dat$hapA[ii] = dat$total_phased[ii] - dat$hapB[ii]
  }
  
  tmp_index = which(dat$phased == 1)
  dat$LBC[tmp_index] = as.numeric(apply(dat[tmp_index,c("total_phased","hapB")],
                                        1,function(x) lchoose(x[1],x[2])))
  dat$LGX1 = sapply(dat$total,function(x) lgamma(x+1))
  dat = dat[,which(!(names(dat) %in% c("pp","Z")))]
  
  return(list(dat,ZZ))
  # save.image("_test_func_tumorEQTL.Rdata")
}


## ---------------------------------------------------------------------
## Simulation 
## ---------------------------------------------------------------------


set.seed(1)
gene_name = "G1"
NN = 300
XX = smart_df(cbind(1,matrix(runif(NN*3),NN,3))); names(XX) = paste0("x",1:4)
numSNPs = 100
BETA = c(5,0.15,-0.5,0.25)
PHI = 0.1
PSI = 0.05
MAF = 0.2
prob_phased = 0.05
corrSNPss = c(0, 0.3, 0.5)
eQTL_index = 1



simDat = list()
for(exp_bxj in 1:2){
  for(corrSNPs in corrSNPss){
    simDat[[paste0("ratioPi_",exp_bxj,"_Cor_", corrSNPs)]] =  
      gen_gene_RC_SNPs(paste0(bxj,"_", corrSNPs),XX,numSNPs,BETA,PHI,log(exp_bxj),PSI,
                       MAF,prob_phased,corrSNPs,eQTL_index) 
  }
}


