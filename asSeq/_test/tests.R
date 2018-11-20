
# --------------------------------------------------------------
# Modififed from DeanB is from R/DCluster version 0.2-2
# Copyright 2004 Virgilio Gómez Rubio and Roger Bivand
#
# Score Test for overdispersion of glm model
# Dean 1992 JASA, vol 87, No, 418, p451-457
#
# There are three test statistics in Dean 1992
# Pa, Pb and Pc (Table 1 of Dean 1992)
# Pb is the one appropriate for Poisson versus Negative Binomial
# or preceisely, NB2 where variance is a quadratic function of 
# mean: Yang Z, Hardin JW, Addy CL, Vuong QH. et al. 
# Testing approaches for overdispersion in poisson regression 
# versus the generalized poisson model. Biom J. 2007
# Aug;49(4):565-84. PubMed PMID: 17638291.
# --------------------------------------------------------------

DeanB <- function(x.glm, alternative="greater") {
	alternative = match.arg(alternative, c("less", "greater", "two.sided"))
  
	if (!(inherits(x.glm, "glm"))) stop("not a glm object")
  
	y  = model.response(model.frame(x.glm))
	mu = fitted(x.glm)
	Pb = sum((y-mu)^2-y)/sqrt(2*sum(mu^2))
	names(Pb) = "P_B"
  
	Pv = NA
	
  if (is.finite(Pb)) {
    if (alternative == "two.sided"){
      Pv = 2 * pnorm(abs(Pb), lower.tail=FALSE)
    }else if (alternative == "greater"){
      Pv = pnorm(Pb, lower.tail=FALSE)
    }else{
      Pv = pnorm(Pb)
    }
	}
	
  res = list(statistic=Pb, p.value=Pv, alternative=alternative, 
              method="Dean's P_B test for overdispersion",
              data.name=deparse(substitute(x.glm)))
  res
}

# --------------------------------------------------------------
# A Score Test for zero-inflation in a Poisson Distribution
# Jan van den Broek, Biometrics 51, 738-743
# --------------------------------------------------------------

ZIPtest <- function(pglm, alternative="greater") {
	alternative = match.arg(alternative, c("less", "greater", "two.sided"))
  
  if (!(inherits(pglm, "glm"))) stop("not a glm object")
  
	if (family(pglm)$family != "poisson") stop("not a glm Poisson object")
  
	y  = model.response(model.frame(pglm))
  X  = pglm$x
  
  n  = length(y)
  w0 = which(y==0)
  n0 = length(w0)
  
  S1   = NA
  pval = NA
  
  if(n0 > 0){
    
    mu = fitted(pglm)
    S1 = sum(exp(mu[w0])) - n
    MM = solve(t(X) %*% diag(mu) %*% X)
    MM = X %*% MM %*% t(X)
    SV = sum(exp(mu)) - n - t(mu) %*% MM %*% mu 
    S1 = S1/sqrt(SV)
      
    if (is.finite(S1)) {
      if (alternative == "two.sided"){
        pval = 2 * pnorm(abs(S1), lower.tail=FALSE)
      }else if (alternative == "greater"){
        pval = pnorm(S1, lower.tail=FALSE)
      }else{
        pval = pnorm(S1)
      }
    }
  
  }
  
	res = list(statistic=S1, p.value=pval, alternative=alternative)
	res
}

# --------------------------------------------------------------
# A Score Test for zero-inflation in a Negative Binomail GLM
# Deng, D. and Paul, S.R. Statistica Sinica Vol. 15, 257-276
# Score tests for zero-inflation and over-dispersion in 
# generalized linear models
# --------------------------------------------------------------

ZINBtest <- function(nbglm) {
  
  if (!(inherits(nbglm, "glm"))) stop("not a glm object")
  
	if (!grepl("Negative Binomial", family(nbglm)$family)){
    stop("not a glm negative binomial object")
  }
  
	y  = nbglm$y
  X  = nbglm$x
  
  n  = length(y)
  w0 = which(y==0)
  n0 = length(w0)
  
  Z4   = NA
  pval = NA
  
  if(n0 > 0){
    
    mu  = fitted(nbglm)
    # phi is the invsers of dispersion parameter, 
    # which is denoted by c in Deng and Paul
    phi = 1/nbglm$theta 
    mu0 = mu[w0]
    
    Z4  = sum((1+phi*mu0)^(1/phi)) - n
    
    W3  = -mu/(1+phi*mu)
    W1  = -diag(W3)

    Irr = sum((1+phi*mu)^(1/phi)) - n
    Irc = sum(log(1+phi*mu)/(phi*phi) - mu/(phi*(1+phi*mu)))
    
    maxI = rep(NA, length(mu))
    for(i in 1:length(mu)){
      maxI[i]= qnbinom(0.99, size=1/phi, mu=mu[i])
    }
    maxY= max(maxI)
    
    l   = 1:maxY
    vv1 = (l-1)/(1+(l-1)*phi)
    vv1 = cumsum(vv1*vv1)
    Ev1 = rep(NA, length(mu))
    
    for(i in 1:length(mu)){
      probI  = dnbinom(1:maxI[i], size=1/phi, mu=mu[i])
      Ev1[i] = sum(probI*vv1[1:maxI[i]])
    }
    
    Icc = sum(Ev1 - 2*log(1 + phi*mu)/(phi^3) - (2*mu + phi*mu*mu)/(phi*phi*(1+phi*mu)))
    
    MM  = X %*% solve(t(X) %*% W1 %*% X) %*% t(X)
    Vr  = Irr - t(W3) %*% MM %*% W3 - Irc*Irc/Icc
    
    if(Vr <=0){ stop("Vr <= 0 \n") }
    
    Z4 = drop(Z4*Z4/Vr)
    
    if (is.finite(Z4)) {
      pval = pchisq(Z4, df=1, lower.tail = FALSE)
    }
    
  }
  
	res = list(statistic=Z4, p.value=pval)
	res
}

