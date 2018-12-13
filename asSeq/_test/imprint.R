imprint <- function(nA, nTotal, maxIt=500, verbose=0) {
    
  # --------------------------------------------------------------
  # logsumexp(v) = log(sum(exp(v)))
  # --------------------------------------------------------------

  logsumexp <- function(v){
    if(any(v==Inf)){
      stop("positive inifinite value in v\n")
    }else{
      w = which(v == -Inf)
      if(length(w)>0){
        v = v[-w]
      }
    }
    if(length(v)==0){
      lse = -Inf
    }else if(length(v)==1){
      lse = v[1]
    }else{
      wv  = which.max(v)
      mv  = max(v)
      res = sum(exp(v[-wv] - mv))
      lse = mv + log(1+res)
    }
    lse
  }

  
  # --------------------------------------------------------------    
  # log density of beta-binomial, or log likelihood for one sample
  # --------------------------------------------------------------
  logBB1 <- function(ni, ni0, pi, theta){
    alpha = pi/theta
    beta  = (1-pi)/theta
    
    lchoose(ni, ni0) + lbeta(ni0+alpha, ni-ni0+beta) - lbeta(alpha, beta)
  }
  
  logBB <- function(ni, ni0, pi, theta){
    
    tmp0 = lchoose(ni, ni0)
    
    if(ni0 > 0){
      seq1 = seq(0, ni0-1, by=1)
      for(k in seq1){
        tmp0 = tmp0 + log(pi + k*theta)
      }
    }
    
    if(ni0 < ni){
      seq2 = seq(0, ni-ni0-1, by=1)
      for(k in seq2){
        tmp0 = tmp0 + log(1 - pi + k*theta)
      }
    }
    
    seq3 = 1:(ni-1)
    for(k in seq3){
      tmp0 = tmp0 - log(1 + k*theta)
    }
    
    tmp0
  }
  
  # --------------------------------------------------------------    
  # Gradient of log beta-binomial
  # --------------------------------------------------------------
  gradBB <- function(ni, ni0, pi, theta){
    
    tmp0 = tmp1 = 0
    
    if(ni0 > 0){
      seq1 = seq(0, ni0-1, by=1)
      for(k in seq1){
        tmp0 = tmp0 + 1/(pi + k*theta)
        tmp1 = tmp1 + k/(pi + k*theta)
      }
    }
    
    if(ni0 < ni){
      seq2 = seq(0, ni-ni0-1, by=1)
      for(k in seq2){
        tmp0 = tmp0 - 1/(1 - pi + k*theta)
        tmp1 = tmp1 + k/(1 - pi + k*theta)
      }
    }
    
    seq3 = 1:(ni-1)
    for(k in seq3){
      tmp1 = tmp1 - k/(1 + k*theta)
    }
    
    c(tmp0, tmp1)
    
  }
  
  # --------------------------------------------------------------    
  # log likelihood: H1
  # --------------------------------------------------------------
  
  logH1 <- function(par){
    pi    = par[1]
    theta = par[2]
    
    sumL = 0
    
    for(i in 1:length(nA)){
      ni   = nTotal[i]
      ni0  = nA[i]

      sumL = sumL + logBB(ni, ni0, pi, theta)
    }
    
    sumL
  }
  
  # --------------------------------------------------------------    
  # Gradient of log likelihood: H1
  # --------------------------------------------------------------
  
  gradLogH1 <- function(par){
    pi    = par[1]
    theta = par[2]
    
    gradPi = gradTheta = 0
    
    for(i in 1:length(nA)){
      ni   = nTotal[i]
      ni0  = nA[i]
      tmp  = gradBB(ni, ni0, pi, theta)
      
      gradPi = gradPi + tmp[1]
      gradTheta = gradTheta + tmp[2]
    }
    
    c(gradPi, gradTheta)
  }
  
  # --------------------------------------------------------------    
  # log likelihood: H3
  # --------------------------------------------------------------
  
  logH3 <- function(par){
    pi0    = par[1]
    theta0 = par[2]
    pi1    = par[3]
    theta1 = par[4]
    
    sumL = 0
    
    for(i in 1:length(nA)){
      ni   = nTotal[i]
      ni0  = nA[i]
      
      sumL = sumL + tau0[i]*logBB(ni, ni0, pi0, theta0)
      sumL = sumL + (1-tau0[i])*logBB(ni, ni0, pi1, theta1)
    }
    
    sumL
  }
  
  # --------------------------------------------------------------    
  # Gradient of log likelihood: H3
  # --------------------------------------------------------------
  
  gradLogH3 <- function(par){
    pi0    = par[1]
    theta0 = par[2]
    pi1    = par[3]
    theta1 = par[4]
    
    gradPi0 = gradTheta0 = 0
    
    for(i in 1:length(nA)){
      ni   = nTotal[i]
      ni0  = nA[i]
      tmp  = gradBB(ni, ni0, pi0, theta0)
      
      gradPi0 = gradPi0 + tau0[i]*tmp[1]
      gradTheta0 = gradTheta0 + tau0[i]*tmp[2]
    }
    
    gradPi1 = gradTheta1 = 0
    
    for(i in 1:length(nA)){
      ni   = nTotal[i]
      ni0  = nA[i]
      tmp  = gradBB(ni, ni0, pi1, theta1)
      
      gradPi1 = gradPi1 + (1 - tau0[i])*tmp[1]
      gradTheta1 = gradTheta1 + (1-tau0[i])*tmp[2]
    }
    
    c(gradPi0, gradTheta0, gradPi1, gradTheta1)
  }
  
  # --------------------------------------------------------------    
  # check input
  # --------------------------------------------------------------
  
  if(!is.numeric(nA)){
    stop("nA must be a numerical vector\n")
  }
  
  if(!is.numeric(nTotal)){
    stop("nTotal must be a numerical vector\n")
  }
  
  if(length(nA) != length(nTotal)){
    stop("nA and nTotal have different lengths\n")  
  }
  
  if(any(nA > nTotal)){
    stop("nA must be smaller than nTotal\n")  
  }

  # -------------------------------------------------------------
  # N is the number of samples
  # -------------------------------------------------------------
  
	N = length(nA)

  # -------------------------------------------------------------
  # First, find MLE under H1
  # -------------------------------------------------------------
  
  par1 = c(0.5, 0.1)
  op1  = optim(par1, logH1, gr=gradLogH1, method="L-BFGS-B", 
               lower=c(0,0) + 1e-16, upper=c(1-1e-16, Inf), 
               control=list(fnscale=-1))

  # -------------------------------------------------------------
  # Secondly, find MLE under H3
  # -------------------------------------------------------------
  
  # -------------------------------------------------------------
  # Initialization:
  # tau0 is the posterior that gene i is imprinted
  # pi0  is the prior that gene i is imprinted
  # -------------------------------------------------------------
  
  pi0  = rep(0.5, N)
  tau0 = rep(NA,  N)
  
  # -------------------------------------------------------------
  # Initialize alpha and beta
  #
  # Tried to divdide the data into two parts based on the ranks 
  # and then use Method of moments to esimate alpha and beta
  # for each group. Well, it does not work well.
  # -------------------------------------------------------------
  
  pi0    = 0.2
  theta0 = 0.1
  pi1    = 0.8
  theta1 = 0.1
  
  parOld = parNew = c(pi0, theta0, pi1, theta1)
  
  oldLik = 1e16
  logLik = Q = NULL
  converged = 0

  # -------------------------------------------------------------
  # iteration
  # -------------------------------------------------------------
  for(it in 1:maxIt){
    
    if(verbose > 0 && it %% verbose == 0){
      message("\niteration ", it, " ", date())
    }
    
    # -----------------------------------------------------------
    # E-step, estimate tau0
    # -----------------------------------------------------------
    
    pi0    = parOld[1]
    theta0 = parOld[2]
    pi1    = parOld[3]
    theta1 = parOld[4]
    
    for(i in 1:N){
      ni  = nTotal[i]
      ni0 = nA[i]
      tmp1 = logBB(ni, ni0, pi1, theta1)
      tmp0 = logBB(ni, ni0, pi0, theta0)
      tau0[i] = 1/(1 + exp(tmp1 - tmp0))
    }
    
    # -----------------------------------------------------------
    # M-step, estimate alpha, beta
    # -----------------------------------------------------------
        
    op3    = optim(parOld, logH3, gr=gradLogH3, method="L-BFGS-B", 
                   lower=c(0,0,0,0) + 1e-16, upper=c(1-1e-16, Inf),
                   control=list(fnscale=-1))
    parNew = op3$par

    # -----------------------------------------------------------
    # calculate the likelihood
    # (1) to check the covergenece of EM algorithm
    # (2) to compare the likelihoods under other situations
    # -----------------------------------------------------------
    
    newLik = 0
    for(i in 1:N){
      pi0    = parOld[1]
      theta0 = parOld[2]
      pi1    = parOld[3]
      theta1 = parOld[4]
      
      ni     = nTotal[i]
      ni0    = nA[i]
      
      tmp0   = log(0.5) + logBB(ni, ni0, pi0, theta0)
      tmp1   = log(0.5) + logBB(ni, ni0, pi1, theta1)
      newLik = newLik + logsumexp(c(tmp0, tmp1))
    }
    
    cond2 = max(abs(parNew - parOld)/abs(parOld), na.rm =TRUE) < 1e-3
    
    if(cond2){
      converged = 1
      break  
    }
    
    logLik = c(logLik, newLik)
    
    parOld = parNew
    oldLik = newLik
  }
  
  if(converged == 0){
    warning("fail to coverge\n")
  }else{
    message("converged after ", it, " iterations\n")
  }
      
  list(parH1=op1$par, logLikH1=op1$value, parH3=op3$par, 
      postPH3=tau0, logLikH3=logLik)
}


