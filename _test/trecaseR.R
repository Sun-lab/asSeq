
loglikJoin <- function(bxj, y, x, nA, nTotal, zeta, mu, b0, phi, theta){
  
  mu1 = mu
  ww2 = which(x==2)
  ww1 = which(x==1)
  
  mu1[ww2] = exp(log(mu[ww2]) + bxj - b0)
  mu1[ww1] = exp(log(mu[ww1]) + log(1 + exp(bxj)) - log(1 + exp(b0)))
  
  pi1 = exp(bxj)/(1 + exp(bxj))
  par = c(theta, pi1)
  
  logTReC = loglikNB(phi, mu1, y)
  
  logASE  = logH1(par, nA, nTotal, zeta)
  
  logTReC + logASE
}

grad.bxj <- function(bxj, y, x, nA, nTotal, zeta, mu, b0, phi, theta)
{
  mu1 = mu
  ww2 = which(x==2)
  ww1 = which(x==1)
  
  mu1[ww2] = exp(log(mu[ww2]) + bxj - b0)
  mu1[ww1] = exp(log(mu[ww1]) + log((1 + exp(bxj))/2) - log((1 + exp(b0))/2))
  
  dg.dmu = y/mu1 - (1 + phi*y)/(1 + phi*mu1)
  
  
  dmu.db = rep(0, length(mu1))
  dmu.db[ww1] = mu1[ww1]*exp(bxj)/(1 + exp(bxj))
  dmu.db[ww2] = mu1[ww2]
  
  w2use  = which(zeta > 0)
  
  pi = exp(bxj)/(1 + exp(bxj))
  
  dh.dpi = 0.0
  for(w1 in w2use){
    if(nA[w1] > 0){
      ks = 0:(nA[w1] - 1)
      dh.dpi = dh.dpi + sum(1/(pi + ks*theta)) 
    }
    
    if(nA[w1] < nTotal[w1]){
      ks = 0:(nTotal[w1] - nA[w1] - 1)
      dh.dpi = dh.dpi - sum(1/(1 - pi + ks*theta))
    }
  }
  
  dpi.db = exp(bxj)/((1 + exp(bxj))^2)
  
  grad = sum(dg.dmu*dmu.db) + dh.dpi*dpi.db
  grad
}

loglikTheta <- function(theta, pi, nA, nTotal, zeta){
  
  sumL = 0
  
  for(i in 1:length(nA)){
    ni   = nTotal[i]
    ni0  = nA[i]
    
    if(zeta[i]){
      sumL = sumL + logBB(ni, ni0, pi, theta)
    }else{
      sumL = sumL + logBB(ni, ni0, 0.5, theta)
    }
  }
  
  sumL
}

loglikNB <- function(phi, mu, y){
  logL = 0.0
  vphi = 1/phi
  lik0 = vphi*log(vphi) - lgamma(vphi)
  
  for(i in 1:length(y)){
    yi   = y[i]
    mui  = mu[i]
    if(yi==0){
      logL = logL + vphi*log(vphi) - vphi*log(vphi + mui)
    }else{
      logL = logL + lgamma(yi + vphi) - lgamma(yi + 1.0) + yi*log(mui) 
      logL = logL - (vphi+yi)*log(vphi+mui) + lik0
    }
  }
  
  logL
}


grad.theta <- function(theta, nA, nTotal, zeta, pi)
{
  grad = 0.0
  
  for(i in 1:length(nA)){
    if(nA[i] > 0){
      ks = 0:(nA[i] - 1)
      
      if(zeta[i] > 0){
        grad = grad + sum(ks/(pi  + ks*theta)) 
      }else{
        grad = grad + sum(ks/(0.5 + ks*theta)) 
      }
    }
    
    if(nA[i] < nTotal[i]){
      ks = 0:(nTotal[i] - nA[i] - 1)
      
      if(zeta[i] > 0){
        grad = grad + sum(ks/(1 - pi + ks*theta))
      }else{
        grad = grad + sum(ks/(0.5 + ks*theta))
      }
    }
    
    ks = 1:(nTotal[i] - 1)
    grad = grad - sum(ks/(1 + ks*theta))
    
  }
  
  grad
}

#-------------------------------------------------------------
# joint model
#-------------------------------------------------------------

trecaseR <- function(y, y1, y2, X, z1, z2, plotIt=FALSE, traceIt=FALSE){
  
  #----------------------------------------------------------
  # initial model fitting
  #---------------------------------------------------------- 

  g0 = glm.nb(y ~ X)
  g0

  nTotal = as.numeric(y1 + y2)
  nA     = as.numeric(y2)
  nA[z2==3] = y1[z2==3]

  wkp    = which(nTotal >=5)
  if(length(wkp) < 5){
    stop("no enough allele specific reads\n")  
  }
  
  nTotal = nTotal[wkp]
  nA     = nA[wkp]

  wkp1   = which(abs(z2[wkp] - 2) < 1.5)      
  zeta   = rep(0, length(nTotal))
  zeta[wkp1] = 1

  a1 = aseR(nA, nTotal, zeta)
  a1

  pi0    = 0.5
  theta0 = a1$parH0
  phi0   = 1/g0$theta
  b0     = 0.0
  mu0    = g0$fitted
  b0_ase = log(pi0/(1-pi0))
  g1     = g0
  
  logLik = NULL
  
  if(traceIt){
    message(sprintf("Initial: b0=%e, pi0=%e, b0_ase=%e, theta0=%e, phi0=%e\n", 
                  b0, pi0, b0_ase, theta0, phi0))
  }
  
  par0       = c(theta0, pi0)
  logLikNULL = 2*logH1(par0, nA, nTotal, zeta) + 2*loglikNB(phi0, mu0, y)

  for(g in 1:100){
    parDiff = 0.0

    #-----------------------------------------------------------
    # estimate bxj
    #-----------------------------------------------------------

    if(plotIt){
      bs = seq(-1.0, 1.0, by=0.1)
      ll = rep(NA, length(bs))

      gs = rep(NA, length(bs))

      for(i in 1:length(bs)){
        bxj   = bs[i]
        gs[i] = grad.bxj(bxj, y, x=z1, nA, nTotal, zeta, mu0, b0, phi0, theta0)
        ll[i] = loglikJoin(bxj, y, x=z1, nA, nTotal, zeta, mu0, b0, phi0, theta0)
      }
      
      quartz()
      par(mfrow=c(2,1))
      plot(bs, gs, type="l")
      abline(h=0)
      plot(bs, ll, type="l")
    }
    
    par0    = c(theta0, pi0)
    logLik0 = 2*logH1(par0, nA, nTotal, zeta) + g1$twologlik
        
    op = optim(par=b0, fn=loglikJoin, gr=grad.bxj, y=y, x=z1, nA=nA, nTotal=nTotal, 
          zeta=zeta, mu=mu0, b0=b0, phi=phi0, theta=theta0, method = "BFGS", 
          control=list(fnscale=-1.0))
    
    b1  = op$par

    pi1 = exp(b1)/(1 + exp(b1))
    
    par1 = c(theta0, pi1)
    mu1  = mu0
    ww2  = which(z1==2)
    ww1  = which(z1==1)
    
    mu1[ww2] = mu0[ww2] * exp(b1 - b0)
    mu1[ww1] = mu0[ww1] * (1 + exp(b1))/(1 + exp(b0))
    
    logLik1 = 2*logH1(par1, nA, nTotal, zeta) + 2*loglikNB(phi0, mu1, y)

    if(traceIt){
      message(sprintf("\ng=%d", g))
      message(sprintf("  logLik_TReC=(%e, %e)", g1$twologlik, 2*loglikNB(phi0, mu1, y)))
      tmp0 = 2*logH1(par0, nA, nTotal, zeta)
      tmp1 = 2*logH1(par1, nA, nTotal, zeta)
      message(sprintf("  logLik_ASE =(%e, %e)", tmp0, tmp1))
    }
  
    if((logLik1 - logLik0)/abs(logLik0) < -0.01){stop("liklihood decreases for bxj\n")}
    
    if(parDiff < abs(b1 - b0)) parDiff = abs(b1 - b0)
    
    b0  = b1
    pi0 = pi1
    mu0 = mu1
    
    #-----------------------------------------------------------
    # estimate theta
    #-----------------------------------------------------------
    if(plotIt){
      ths = seq(0,0.1,by=0.001)
      gs  = rep(NA, length(ths))
      ll  = rep(NA, length(ths))
      
      for(i in 1:length(ths)){
        gs[i] = grad.theta(ths[i], nA, nTotal, zeta, pi0)
        ll[i] = logH1(par=c(ths[i], pi0), nA, nTotal, zeta)
      }
      
      quartz()
      par(mfrow=c(2,1))
      plot(ths, gs, type="l")
      abline(h=0)
      plot(ths, ll, type="l")
    }
    
    # uk = uniroot(grad.theta, c(0, 1), nA=nA, nTotal=nTotal, 
    #               zeta=zeta, pi=pi0, tol=1e-10)
    # theta1 = uk$root
    op = optimize(loglikTheta, interval=c(0, 1000), pi=pi0, nA=nA, 
             nTotal=nTotal, zeta=zeta, maximum=TRUE)
    
    theta1 = op$maximum
      
    llase0 = 2*loglikTheta(theta0, pi0, nA, nTotal, zeta)
    llase1 = 2*loglikTheta(theta1, pi0, nA, nTotal, zeta)
    
    if(traceIt){
      message(sprintf("\ng=%d, estimate theta", g))
      message(sprintf("  (llase0, llase1)=(%e, %e)", llase0, llase1))
    }
    
    if((llase1 - llase0)/abs(llase0) < -0.01){stop("liklihood decreases for theta\n")}

    if(parDiff < abs(theta1 - theta0)) parDiff = abs(theta1 - theta0)
    theta0 = theta1
    
    #-----------------------------------------------------------
    # estimate mu and phi
    #-----------------------------------------------------------
    
    bx1  = rep(0, length(y))
    ww1  = which(z1==1)
    ww2  = which(z1==2)
    bx1[ww1] = log((1 + exp(b0))/2)
    bx1[ww2] = b0
    g1   = glm.nb(y ~ X + offset(bx1))
    mu1  = g1$fitted
    phi1 = 1/g1$theta
    
    lltrec0 = 2*loglikNB(phi0, mu0, y)
    lltrec1 = 2*loglikNB(phi1, mu1, y)
    
    if(traceIt){
      message(sprintf("\ng=%d, estimate phi", g))
      message(sprintf("  (lltrec0, lltrec1)=(%e, %e)", lltrec0, lltrec1))
    }
    
    if((lltrec1 - lltrec0)/abs(lltrec0) < -0.01){stop("liklihood decreases for phi")}
    
    if(parDiff < abs(phi1 - phi0)) parDiff = abs(phi1 - phi0)
    
    phi0 = phi1
    mu0  = mu1
    
    #-----------------------------------------------------------
    # check convergence
    #-----------------------------------------------------------

    par0 = c(theta0, pi0)
    logLik = c(logLik, 2*logH1(par0, nA, nTotal, zeta) + g1$twologlik)
              
    if(traceIt){
      message(sprintf("g=%d, parDiff=%.e, b0=%e, theta0=%e, phi0=%e", 
              g, parDiff, b0, theta0, phi0))
    }
              
    if(parDiff < 1e-8) break
    
  }
  
  list(b=b0, theta=theta0, phi=phi0, logLik=logLik,
       lrt=logLik[length(logLik)] - logLikNULL, 
       lrtASE=2.0*(a1$logLikH1 - a1$logLikH0))
}
