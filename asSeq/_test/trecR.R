
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

logLTReC <- function(bxj, y, x, mu, b0, phi, fam){
  mu1 = mu
  ww2 = which(x==2)
  ww1 = which(x==1)
  
  mu1[ww2] = mu[ww2] * exp(bxj - b0)
  mu1[ww1] = mu[ww1] * (1 + exp(bxj))/(1 + exp(b0))
  
  if(fam=="negbin"){
    logL = loglikNB(phi, mu1, y)
  }else{
    logL = 0.0
    for(i in 1:length(y)){
      logL = logL + dpois(y[i], mu1[i], log = TRUE)  
    }
  }
  logL
}

grad.bxj.trec <- function(bxj, y, x, mu, b0, phi, fam)
{
  mu1 = mu
  ww2 = which(x==2)
  ww1 = which(x==1)

  mu1[ww2] = mu[ww2] * exp(bxj - b0)
  mu1[ww1] = mu[ww1] * (1 + exp(bxj))/(1 + exp(b0))
  
  dg.dmu = y/mu1 - (1 + phi*y)/(1 + phi*mu1)
  
  dmu.db = rep(0, length(mu1))
  dmu.db[ww1] = mu1[ww1]*exp(bxj)/(1 + exp(bxj))
  dmu.db[ww2] = mu1[ww2]
  
  grad = sum(dg.dmu*dmu.db)
  grad
}

#-------------------------------------------------------------
# TReC model
#-------------------------------------------------------------

trecR <- function(y, X, z1, fam, nIter=100, plotIt=FALSE, trace=FALSE, yfit=FALSE){
  
  #----------------------------------------------------------
  # initial model fitting
  #----------------------------------------------------------
  if(fam=="negbin"){
    g1      = glm.nb(y ~ X)
    phi0    = 1/g1$theta
    logLik0 = g1$twologlik
  }else{
    g1      = glm(y ~ X, family=poisson())
    phi0    = 0.0
    logLik0 = 2*logLik(g1)
  }
  
  b0   = 0.0
  mu0  = g1$fitted
    
  logLik = logLikOld =logLik0
  
  if(trace)
    message(sprintf("Initial: b0=%e, phi0=%e\n", b0, phi0))

  for(g in 1:nIter){
    
    parDiff = 0.0

    #--------------------------------------------------------
    # estimate bxj
    #--------------------------------------------------------

    if(plotIt){
      bs = seq(-0.3, 0.3, by=0.01)
      ll = rep(NA, length(bs))
      gs = rep(NA, length(bs))

      for(i in 1:length(bs)){
        bxj   = bs[i]
        gs[i] = grad.bxj.trec(bxj, y, x=z1, mu0, b0, phi0, fam)
        ll[i] = logLTReC(bxj, y, x=z1, mu0, b0, phi0, fam)
      }
      
      quartz()
      par(mfrow=c(2,1))
      plot(bs, gs, type="l")
      abline(h=0)
      plot(bs, ll, type="l")
    }
    
    o1 = optim(b0, fn=logLTReC, gr=grad.bxj.trec, y=y, x=z1, mu=mu0, b0=b0, phi=phi0, 
          fam=fam, method="BFGS", control=list(fnscale=-1.0, trace=0))
    

    if(o1$convergence != 0){ 
      warning("g =", g, " fail BFGS\n")
      return(NULL)
    }
    
    if(plotIt){
      abline(v=o1$par)
    }
    
    b1 = o1$par
    grad.bxj.trec(b1, y, z1, mu0, b0, phi0, fam)

    mu1  = mu0
    ww2  = which(z1==2)
    ww1  = which(z1==1)
    
    mu1[ww2] = mu0[ww2] * exp(b1 - b0)
    mu1[ww1] = mu0[ww1] * (1 + exp(b1))/(1 + exp(b0))
    
    logLik1 = 2.0*o1$value

    if(trace){
      message(sprintf("\ng=%d", g))
      message(sprintf("  Estimate bxj: logLik_TReC=(%e, %e)", logLikOld, logLik1))
    }
    
    if((logLik1 - logLikOld)/abs(logLikOld) < -1e-4){stop("liklihood decreases for bxj\n")}
    
    if(parDiff < abs(b1 - b0)) parDiff = abs(b1 - b0)
    
    b0  = b1
    mu0 = mu1
    logLikOld = logLik1
    
    #-----------------------------------------------------------
    # estimate mu and phi
    #-----------------------------------------------------------
    
    bx1  = rep(0, length(y))
    ww1  = which(z1==1)
    ww2  = which(z1==2)
    bx1[ww1] = log((1 + exp(b0))/2)
    bx1[ww2] = b0
    
    if(fam=="negbin"){
      g1   = glm.nb(y ~ X + offset(bx1))
      phi1 = 1/g1$theta
      logLik1 = g1$twologlik
    }else{
      g1 = glm(y ~ X, family=poisson(), offset=bx1)
      phi1 = 0.0
      logLik1 = 2*logLik(g1)
    }
    
    mu1  = g1$fitted
    
    if((logLik1 - logLikOld)/abs(logLikOld) < -1e-5){stop("liklihood decreases for phi\n")}
    
    if(parDiff < abs(phi1 - phi0)) parDiff = abs(phi1 - phi0)
    
    if(trace){
      message(sprintf("  Estimate mu & phi: phi=%.2f, logLik_TReC=(%e, %e)", phi1, logLikOld, logLik1))
    }
    
    phi0 = phi1
    mu0  = mu1
    logLikOld = logLik1

    #-----------------------------------------------------------
    # check convergence
    #-----------------------------------------------------------

    logLik = c(logLik, logLik1)
    
    if(trace){
      message(sprintf("g=%d, parDiff=%.e, b0=%e, phi0=%e", 
                    g, parDiff, b0, phi0))
    }
    
    if(parDiff < 5e-5) break
    
  }
  
  if(yfit){
    l1 = list(b=b0, phi=phi0, logLik=logLik, lrt=logLik1 - logLik0, fitted=g1$fitted)
  }else{
    l1 = list(b=b0, phi=phi0, logLik=logLik, lrt=logLik1 - logLik0)
  }
  
  l1
}
