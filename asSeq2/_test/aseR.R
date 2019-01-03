
# --------------------------------------------------------------    
# log density of beta-binomial, or log likelihood for one sample
# logBB1 is the concise version using log beta function
# logBB  is another version avoiding beta function
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
  
  grTh = grPi = 0
  
  if(ni0 > 0){
    seq1 = seq(0, ni0-1, by=1)
    for(k in seq1){
      xx = 1/(pi + k*theta)
      grPi = grPi + xx
      grTh = grTh + k*xx
    }
  }
  
  if(ni0 < ni){
    seq2 = seq(0, ni-ni0-1, by=1)
    for(k in seq2){
      xx = 1/(1 - pi + k*theta)
      grPi = grPi - xx
      grTh = grTh + k*xx
    }
  }
  
  seq3 = 1:(ni-1)
  for(k in seq3){
    grTh = grTh - k/(1 + k*theta)
  }
  
  c(grTh, grPi)
}

# --------------------------------------------------------------    
# log likelihood: H0
# --------------------------------------------------------------

logH0 <- function(par, nA, nTotal){
  theta = par[1]
  pi0   = 0.5
  
  sumL = 0
  
  for(i in 1:length(nA)){
    ni   = nTotal[i]
    ni0  = nA[i]
    sumL = sumL + logBB(ni, ni0, pi0, theta)
  }
  
  sumL
}

# --------------------------------------------------------------    
# Gradient of log likelihood: H0
# --------------------------------------------------------------

gradLogH0 <- function(par, nA, nTotal){
  theta = par[1]
  pi0   = 0.5
  
  gradTh = 0
  
  for(i in 1:length(nA)){
    ni   = nTotal[i]
    ni0  = nA[i]
    gri  = gradBB(ni, ni0, pi0, theta)
    
    gradTh = gradTh + gri[1]
  }
  
  gradTh
}

# --------------------------------------------------------------    
# log likelihood: H1
# --------------------------------------------------------------

logH1 <- function(par, nA, nTotal, zeta){
  theta = par[1]
  pi1   = par[2]
  
  sumL  = 0
  
  for(i in 1:length(nA)){
    ni  = nTotal[i]
    ni0 = nA[i]
    
    if(zeta[i]){
      sumL = sumL + logBB(ni, ni0, pi1, theta)
    }else{
      sumL = sumL + logBB(ni, ni0, 0.5, theta)
    }
  }
  
  sumL
}

# --------------------------------------------------------------    
# Gradient of log likelihood: H1
# --------------------------------------------------------------

gradLogH1 <- function(par, nA, nTotal, zeta){
  theta = par[1]
  pi1   = par[2]
  
  gradPi1 = gradTh = 0
  
  for(i in 1:length(nA)){
    ni   = nTotal[i]
    ni0  = nA[i]
    
    if(zeta[i]){
      gri     = gradBB(ni, ni0, pi1, theta)
      gradTh  = gradTh  + gri[1]
      gradPi1 = gradPi1 + gri[2]
    }else{
      gri     = gradBB(ni, ni0, 0.5, theta)
      gradTh  = gradTh + gri[1]
    }
    
  }
  
  c(gradTh, gradPi1)
}

# --------------------------------------------------------------    
# main function
# --------------------------------------------------------------

aseR <- function(nA, nTotal, zeta, maxIt=50, trace=0) {    
  
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
  # First, find MLE under H0
  # -------------------------------------------------------------
  
  par0 = 0.1
  op0  = optim(par0, logH0, gr=gradLogH0, nA=nA, nTotal=nTotal, 
               method="L-BFGS-B", lower=1e-16, upper=Inf, 
               control=list(fnscale=-1))
  
  theta0 = op0$par
  
  # -------------------------------------------------------------
  # Next, find MLE under H1
  # -------------------------------------------------------------
  
  par1 = c(theta0, 0.5)
  op1  = optim(par1, logH1, gr=gradLogH1, nA=nA, nTotal=nTotal, zeta=zeta, 
               method="L-BFGS-B", lower=c(0,0) + 1e-16, upper=c(Inf, 1-1e-16), 
               control=list(fnscale=-1))
  
  list(parH0=theta0,  logLikH0=op0$value, 
       parH1=op1$par, logLikH1=op1$value)
}


