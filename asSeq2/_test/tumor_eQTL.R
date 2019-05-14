# ----------------------------------------------------------------------
# Negative Binomial (TREC)
# ----------------------------------------------------------------------

compute_offset <- function(z, RHO, KAPPA, ETA, GAMMA, tau1, tau2){
  n = length(z)
  offsets = vector('numeric', n)
  indAA = which(z == 1)
  offsets[indAA] = log(2*(1-RHO[indAA]) +
                         (tau1[indAA]+tau2[indAA])*RHO[indAA]*KAPPA)

  indAB = which(z == 2)
  offsets[indAB] = log((1-RHO[indAB]) + tau1[indAB]*RHO[indAB]*KAPPA +
                         (1-RHO[indAB])*ETA + tau2[indAB]*RHO[indAB]*KAPPA*GAMMA)

  indBA = which(z == 3)
  offsets[indBA] = log((1-RHO[indBA]) + tau2[indBA]*RHO[indBA]*KAPPA +
                         (1-RHO[indBA])*ETA + tau1[indBA]*RHO[indBA]*KAPPA*GAMMA)

  indBB = which(z == 4)
  offsets[indBB] = log(2*(1-RHO[indBB])*ETA +
                         (tau1[indBB]+tau2[indBB])*RHO[indBB]*KAPPA*GAMMA)

  return(offsets)
}

loglikNB <- function(para, H0, y,z, X, betas, phi, RHO, tau1, tau2){
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

  if(is.null(dim(X)))
    X = data.matrix(X)

  offsets <- compute_offset(z=z, RHO=RHO ,KAPPA= KAPPA, ETA=ETA,
                            GAMMA=GAMMA, tau1=tau1, tau2=tau2)
  expx.t.betas = exp(X%*%betas)
  mu <- expx.t.betas * exp(offsets)

  vphi = 1/phi
  loglik = lgamma(y + vphi) - lgamma(vphi) - lgamma(y + 1) +
    vphi * log(vphi) - (vphi + y)*log(vphi + mu) + y * log(mu)
  return(-sum(loglik))
}

Grad_NB <- function(para, H0, y,z, X, betas, phi, RHO, tau1, tau2){

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

  n = length(RHO)
  if(is.null(dim(X)))
    X = data.matrix(X)

  offsets <- compute_offset(z=z, RHO=RHO ,KAPPA= KAPPA, ETA=ETA,
                            GAMMA=GAMMA, tau1=tau1, tau2=tau2)
  expx.t.betas = exp(X%*%betas)
  mu <- expx.t.betas * exp(offsets)

  dlTRec.dmu = y/mu - (1 + phi * y)/(1 + phi * mu)
  RHO_1 = 1 - RHO
  tau = tau1 + tau2

  dmu.dKAPPA =dmu.dETA =dmu.dGAMMA = rep(0, n)

  # z = AA
  indAA = which(z == 1)
  dmu.dKAPPA[indAA] = tau[indAA]

  # z = AB
  indAB = which(z == 2)
  dmu.dKAPPA[indAB] =  tau1[indAB]+tau2[indAB]*GAMMA
  dmu.dETA[indAB] = 1
  dmu.dGAMMA[indAB] = tau2[indAB]

  # z = BA
  indBA = which(z == 3)
  dmu.dKAPPA[indBA] = tau2[indBA] + tau1[indBA]*GAMMA
  dmu.dETA[indBA] = 1
  dmu.dGAMMA[indBA] = tau1[indBA]

  # z = BB
  indBB = which(z == 4)
  dmu.dKAPPA[indBB] = tau[indBB]*GAMMA
  dmu.dETA[indBB] = 2
  dmu.dGAMMA[indBB] = tau[indBB]

  dmu.dKAPPA = expx.t.betas*RHO*dmu.dKAPPA
  dmu.dETA = expx.t.betas*RHO_1*dmu.dETA
  dmu.dGAMMA = expx.t.betas*RHO*KAPPA*dmu.dGAMMA

  dlTRec.dlKAPPA = sum(dlTRec.dmu * dmu.dKAPPA)*KAPPA
  dlTRec.dlETA = sum(dlTRec.dmu * dmu.dETA)*ETA
  dlTRec.dlGAMMA = sum(dlTRec.dmu * dmu.dGAMMA)*GAMMA


  if(H0 == 0){
    return(c(-dlTRec.dlKAPPA, -dlTRec.dlETA, -dlTRec.dlGAMMA))
  } else if(H0 == 1){
    return(c(-dlTRec.dlKAPPA, -dlTRec.dlGAMMA))
  } else if(H0 == 2){
    return(c(-dlTRec.dlKAPPA, -dlTRec.dlETA))
  }


}
# ----------------------------------------------------------------------
# TREC
# ----------------------------------------------------------------------

# update betas and phi
update_beta_phi <- function(y, X , offsets){
  if(dim(X)[2] == 1){
    nb_out =glm.nb(y~offset(offsets))
  }else{
    X = X[,-1]
    nb_out =glm.nb(y~X+offset(offsets))
  }
  return(list(betas = as.numeric(nb_out$coefficients),
              phi = 1/nb_out$theta))
}

# updtade KAPPA, ETA, GAMMA
# L-BFGS-B

update_keg_trec <- function(H0, para, y,z, X,betas, b, phi, RHO, tau1, tau2){
  # gr = Grad_NB
  if(H0 == 0){
    lbfgs = optim(par = para , fn = loglikNB,gr = Grad_NB,
                  H0 =0,y =y, z = z, X=X, betas = betas, phi = phi,
                  RHO = RHO, tau1= tau1, tau2 = tau2,
                  method = "L-BFGS-B" , lower = -50, upper = 10)
  }
  if(H0 == 1){
    lbfgs = optim(par = para , fn = loglikNB,gr = Grad_NB,
                  H0 =1, y =y, z = z, X=X, betas = betas, phi = phi,
                  RHO = RHO, tau1= tau1, tau2 = tau2,
                  method = "L-BFGS-B" , lower = -50, upper = 10)
  }
  if(H0 == 2){
    lbfgs = optim(par = para , fn = loglikNB,gr = Grad_NB,
                  H0 =2, y =y, z = z, X=X, betas = betas, phi = phi,
                  RHO = RHO, tau1= tau1, tau2 = tau2,
                  method = "L-BFGS-B" , lower = -50, upper = 10)
  }
  return(list(para_cur = lbfgs$par, loglik = lbfgs$value,
              convergence = lbfgs$convergence))
}

TReC_sfit = function(H0, para, y,z, X, RHO, tau1, tau2, maxiter =500,
                     tol = 1e-7)
{
  iter = 0
  betas = phi = 0
  if(H0 == 0){
    KAPPA = exp(para[1])
    ETA = exp(para[2])
    GAMMA = exp(para[3])
  } else if(H0 == 1){
    KAPPA = exp(para[1])
    ETA = 1.0;
    GAMMA = exp(para[2])
  } else if(H0 == 2){
    KAPPA = exp(para[1])
    ETA = exp(para[2])
    GAMMA = 1.0
  }

  offsets <- compute_offset(z=z, RHO=RHO ,KAPPA= KAPPA, ETA=ETA,
                            GAMMA=GAMMA, tau1=tau1, tau2=tau2)
  if(is.null(dim(X)))
    X = data.matrix(X)
  #1 initialize phi and betas
  bp = update_beta_phi(y,X, offsets = offsets )
  beta_cur = bp$betas
  phi_cur =  bp$phi

  #2 initial value for KAPPA, ETA, GAMMA
  if(H0 == 0){
    para = c(0,0,0)
  }else{
    para = c(0,0)
  }

  keg_lbfgs = update_keg_trec(H0=H0, para = para , y =y ,z=z, X=X,
                              betas = beta_cur, phi=phi_cur, RHO=RHO,
                              tau1=tau1, tau2=tau2)

  if(keg_lbfgs$convergence != 0)
    message(sprintf('lbfgs error code: %s' ,keg_lbfgs$convergence))

  para_cur = keg_lbfgs$para_cur

  # update parameters
  while(iter < maxiter & (min(abs(para_cur- para)) > tol |
                          min(abs(beta_cur - betas)) > tol |
                          abs(phi_cur - phi) >tol )){
    para = para_cur
    betas = beta_cur
    phi = phi_cur
    if(H0 == 0){
      KAPPA = exp(para_cur[1])
      ETA = exp(para_cur[2])
      GAMMA = exp(para_cur[3])
    } else if(H0 == 1){
      KAPPA = exp(para_cur[1])
      ETA = 1.0;
      GAMMA = exp(para_cur[2])
    } else if(H0 == 2){
      KAPPA = exp(para_cur[1])
      ETA = exp(para_cur[2])
      GAMMA = 1.0
    }

    offsets <- compute_offset(z=z, RHO=RHO ,KAPPA= KAPPA, ETA=ETA,
                              GAMMA=GAMMA, tau1=tau1, tau2=tau2)
    bp = update_beta_phi(y,X, offsets = offsets )
    beta_cur = bp$betas
    phi_cur =  bp$phi


    keg_lbfgs = update_keg_trec(H0=H0, para = para_cur , y =y ,z=z, X=X,
                                betas = beta_cur, phi=phi_cur, RHO=RHO,
                                tau1=tau1, tau2=tau2)

    para_cur = keg_lbfgs$para_cur

    if(keg_lbfgs$convergence != 0)
      message(sprintf('lbfgs error code: %s' ,keg_lbfgs$convergence))

    iter = iter + 1
  }

  if(H0 == 0){
    KAPPA = exp(para_cur[1])
    ETA = exp(para_cur[2])
    GAMMA = exp(para_cur[3])
  } else if(H0 == 1){
    KAPPA = exp(para_cur[1])
    ETA = 1.0;
    GAMMA = exp(para_cur[2])
  } else if(H0 == 2){
    KAPPA = exp(para_cur[1])
    ETA = exp(para_cur[2])
    GAMMA = 1.0
  }
  return(list(KAPPA = KAPPA, ETA = ETA, GAMMA = GAMMA, phi = phi_cur,
              betas = beta_cur, loglik =keg_lbfgs$loglik, iters = iter ))
}

TReC_test = function(y,z, X, RHO, tau1, tau2, maxiter =500, tol = 1e-7){
  H0 =0
  para = c(0,0,0)
  res1 = TReC_sfit(H0, para, y,z, X, RHO, tau1, tau2, maxiter=maxiter, tol=tol)

  H0 =1
  para = c(0,0)
  res1.eta1 = TReC_sfit(H0, para, y,z, X, RHO, tau1, tau2, maxiter=maxiter,
                        tol=tol)

  H0=2
  para = c(0,0)
  res1.gamma1 = TReC_sfit(H0, para, y,z, X, RHO, tau1, tau2, maxiter=maxiter,
                          tol=tol)

  H0=0
  para = c(res1.gamma1$KAPPA, res1.gamma1$ETA, res1.eta1$GAMMA)
  res1b = TReC_sfit(H0, para, y,z, X, RHO, tau1, tau2, maxiter=maxiter, tol=tol)

  if(res1b$loglik < res1$loglik)
    res1 = res1b
  p.eta = pchisq(2*( res1.eta1$loglik- res1$loglik ), df =1, lower.tail = F)
  p.gamma = pchisq(2*( res1.gamma1$loglik- res1$loglik ), df =1, lower.tail = F)

  return(list(p.eta = p.eta, p.gamma = p.gamma, loglik.full = - res1$loglik,
              loglik.eta = - res1.eta1$loglik,
              loglik.gamma = - res1.gamma1$loglik,
              KAPPA = res1$KAPPA, ETA = res1$ETA, GAMMA = res1$GAMMA,
              phi = res1$phi, betas = res1$betas
              #,iters = c(res1$iters, res1.eta1$iters, res1.gamma1$iters)
  ))
}

# ----------------------------------------------------------------------
# betas-Binomial
# ----------------------------------------------------------------------
# tauB is the ASCN corresponding to allele B and is one of tau1, tau2

compute_pi <- function(z_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau){
  pis = vector('numeric', length(z_AS))
  indhomo = which(z_AS %in% c(1,4))
  pis[indhomo] = (RHO_AS[indhomo]*tauB[indhomo]*KAPPA + 1 -
                    RHO_AS[indhomo])/(RHO_AS[indhomo]*tau[indhomo]*KAPPA
                                      + 2*(1 - RHO_AS[indhomo]))
  indhet = which(z_AS %in% c(2,3))
  tmp1 = RHO_AS[indhet]*tauB[indhet]*GAMMA*KAPPA + (1 - RHO_AS[indhet])*ETA
  pis[indhet] = tmp1/(RHO_AS[indhet]*(tau[indhet]-tauB[indhet])*KAPPA +
                        (1 - RHO_AS[indhet]) + tmp1)
  return(pis)
}

loglikBB <- function(para,H0, A,D,THETA, z_AS, RHO_AS, tauB, tau){

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

  pis <- compute_pi(z_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau)

  lASE = mapply(lchoose,D,A) # can be saved as input
  lASE = lASE + lgamma(A+pis/THETA)+lgamma(D-A+(1-pis)/THETA)+lgamma(1.0/THETA)-
    lgamma(D+1.0/THETA)-lgamma(pis/THETA)-lgamma((1-pis)/THETA);
  return(-sum(lASE))
}

Grad_BB_keg<-function(para, H0,A,D,THETA, z_AS, RHO_AS,tauB, tau){
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

  pis <- compute_pi(z_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau)

  n = length(RHO_AS)
  pis_1 = 1 - pis
  dlASE.dpi = (1/THETA)*(digamma(A+pis/THETA) + digamma(pis_1/THETA) -
                           digamma(D-A+pis_1/THETA) - digamma(pis/THETA))

  dpi.dKAPPA = dpi.dETA = dpi.dGAMMA = rep(0, n)
  indhomo = which((z_AS %in% c(1,4)))
  tmp.homo = (RHO_AS[indhomo]*tau[indhomo]*KAPPA + 2*(1-RHO_AS[indhomo]))
  dpi.dKAPPA[indhomo] = RHO_AS[indhomo]*tauB[indhomo]/tmp.homo -
    RHO_AS[indhomo]*tau[indhomo]*
    (RHO_AS[indhomo]*tauB[indhomo]*KAPPA+(1-RHO_AS[indhomo]))/tmp.homo^2


  indhet = which(z_AS %in% c(2,3))
  tmp1 = RHO_AS[indhet]*tauB[indhet]*GAMMA*KAPPA+(1-RHO_AS[indhet])*ETA
  tmp2 = RHO_AS[indhet]*(tau[indhet]-tauB[indhet])*KAPPA + (1 - RHO_AS[indhet]) +
    tmp1
  tmp3 = tmp1/tmp2^2

  dpi.dKAPPA[indhet] = RHO_AS[indhet]*tauB[indhet]*GAMMA/tmp2 -
    (RHO_AS[indhet]*(tau[indhet]-tauB[indhet])+
       RHO_AS[indhet]*tauB[indhet]*GAMMA)*tmp3
  dpi.dETA[indhet] = (1-RHO_AS[indhet])/tmp2 - (1-RHO_AS[indhet])*tmp3
  dpi.dGAMMA[indhet] =RHO_AS[indhet]*tauB[indhet]*KAPPA/tmp2 -
    (RHO_AS[indhet]*tauB[indhet]*KAPPA)*tmp3


  dlASE.dlKAPPA = sum(dlASE.dpi * dpi.dKAPPA * KAPPA)
  dlASE.dlETA = sum(dlASE.dpi * dpi.dETA * ETA)
  dlASE.dlGAMMA = sum(dlASE.dpi * dpi.dGAMMA * GAMMA)

  if(H0 == 0){
    return(c(-dlASE.dlKAPPA, -dlASE.dlETA, -dlASE.dlGAMMA))
  } else if(H0 == 1){

    return(c(-dlASE.dlKAPPA, -dlASE.dlGAMMA))

  } else if(H0 == 2){
    return(c(-dlASE.dlKAPPA, -dlASE.dlETA))
  }

}




loglikBB_THETA <- function(THETA,A,D,pis){
  lASE = mapply(lchoose,D,A)
  lASE = lASE + lgamma(A+pis/THETA)+lgamma(D-A+(1-pis)/THETA)+lgamma(1.0/THETA)-
    lgamma(D+1.0/THETA)-lgamma(pis/THETA)-lgamma((1-pis)/THETA)
  return(-sum(lASE))
}

lASE.dTHETA <- function(THETA,A,D,pis){
  lASE.THETA = (1.0/(THETA*THETA))*(digamma(D+1.0/THETA)+digamma(pis/THETA)*pis+
                                      digamma((1.0-pis)/THETA)*(1.0-pis)-
                                      digamma(A+pis/THETA)*pis-digamma(D-A+(1.0-pis)/THETA)*
                                      (1.0-pis)-digamma(1.0/THETA))
  return(-sum(lASE.THETA))
}


# ----------------------------------------------------------------------
# TRECASE
# ----------------------------------------------------------------------

loglik.TReCASE <- function(para, H0, y,z, z_AS, X, betas, phi, RHO,
                           RHO_AS, tau1, tau2, A,D, THETA, tauB, tau){
  lNB <- loglikNB(para=para, H0=H0, y=y, z=z, X=X, betas=betas, phi=phi,
                  RHO=RHO, tau1=tau1, tau2=tau2)
  lASE <- loglikBB(para=para,H0=H0,A=A,D=D,THETA=THETA, z_AS=z_AS,
                   RHO=RHO_AS, tauB=tauB, tau = tau)
  return(lNB + lASE)

}
grad.TReCASE.keg <- function(para, H0, y,z, z_AS, X, betas, phi, RHO,
                             RHO_AS, tau1, tau2, A,D, THETA, tauB, tau){
  return(Grad_NB(para=para, H0=H0, y=y, z=z, X=X, betas=betas, phi=phi,
                 RHO=RHO, tau1=tau1, tau2=tau2) +
           Grad_BB_keg(para=para,H0,A=A,D=D,THETA=THETA, z_AS=z_AS,
                       RHO=RHO_AS, tauB=tauB, tau = tau))
}

update_keg_trecase <- function(para, H0, y,z, z_AS, X, betas, phi, RHO,
                               RHO_AS, tau1, tau2, A,D, THETA, tauB, tau){
  if(H0 == 0){
    lbfgs = optim(par = para, fn = loglik.TReCASE, gr = grad.TReCASE.keg, H0 =0,
                  y =y, z = z, z_AS=z_AS, X=X, betas = betas, phi = phi,
                  RHO = RHO, RHO_AS=RHO_AS,  tau1= tau1,
                  tau2 = tau2,A=A,D=D,THETA=THETA,tauB=tauB, tau=tau,
                  lower = -50, upper = 10,method = "L-BFGS-B")
  }
  if(H0 == 1){
    lbfgs = optim(par = para, fn = loglik.TReCASE, gr = grad.TReCASE.keg, H0 =1,
                  y =y, z = z, z_AS=z_AS, X=X, betas = betas, phi = phi,
                  RHO = RHO, RHO_AS=RHO_AS,  tau1= tau1,
                  tau2 = tau2,A=A,D=D,THETA=THETA, tauB=tauB, tau=tau,
                  lower = -50, upper = 10,method = "L-BFGS-B")
  }
  if(H0 == 2){
    lbfgs = optim(par = para, fn = loglik.TReCASE,gr = grad.TReCASE.keg, H0 =2,
                  y =y, z = z, z_AS=z_AS, X=X, betas = betas, phi = phi,
                  RHO = RHO, RHO_AS=RHO_AS,  tau1= tau1,
                  tau2 = tau2,A=A,D=D,THETA=THETA,tauB=tauB, tau=tau,
                  lower = -50, upper = 10,method = "L-BFGS-B")
  }
  return(lbfgs)
}

#have to figure out why adding gradient will fail the function: fixed
update_theta <- function(para, H0, z_AS, RHO_AS, A,D, THETA, tauB, tau){
  if(H0 == 0){
    KAPPA = exp(para[1])
    ETA = exp(para[2])
    GAMMA = exp(para[3])
  } else if(H0 == 1){
    KAPPA = exp(para[1])
    ETA = 1.0;
    GAMMA = exp(para[2])
  } else if(H0 == 2){
    KAPPA = exp(para[1])
    ETA = exp(para[2])
    GAMMA = 1.0
  }
  pis <- compute_pi(z_AS=z_AS, RHO_AS=RHO_AS,
                    KAPPA=KAPPA, ETA=ETA, GAMMA=GAMMA,
                    tauB=tauB, tau=tau)
  lbfgs = optim(par = THETA , fn = loglikBB_THETA, gr = lASE.dTHETA,
                A=A,D=D, pis = pis,
                lower = 1e-8, upper = 50,method = "L-BFGS-B")
}

# tried uniroot
# uniroot(lASE.dTHETA,c(1e-9, 2), A=A,D=D, pis = pis )

TReCASE_sfit = function(H0, para, y,z,z_AS, X, RHO, RHO_AS,
                        tau1, tau2, A, D, tauB, tau, maxiter =500, tol = 1e-7){
  iter = 0
  #1 initialize all the parameters
  betas =  0
  THETA = phi = 0.05
  if(H0 == 0){
    KAPPA = exp(para[1])
    ETA = exp(para[2])
    GAMMA = exp(para[3])
  } else if(H0 == 1){
    KAPPA = exp(para[1])
    ETA = 1.0;
    GAMMA = exp(para[2])
  } else if(H0 == 2){
    KAPPA = exp(para[1])
    ETA = exp(para[2])
    GAMMA = 1.0
  }

  offsets <- compute_offset(z=z, RHO=RHO ,KAPPA= KAPPA, ETA=ETA,
                            GAMMA=GAMMA, tau1=tau1, tau2=tau2)
  if(is.null(dim(X)))
    X = data.matrix(X)

  # update beta phi
  bp = update_beta_phi(y,X, offsets = offsets )
  beta_cur = bp$betas
  phi_cur =  bp$phi

  # update theta
  theta_lbfgs = update_theta(para=para, H0=H0, z_AS=z_AS, RHO_AS=RHO_AS,
                             A=A, D=D, THETA=THETA, tauB=tauB, tau=tau)

  THETA_cur = theta_lbfgs$par

  # update keg
  keg_lbfgs = update_keg_trecase(H0=H0, para=para, y =y, z=z, z_AS=z_AS, X=X,
                                 betas=beta_cur, phi=phi_cur, RHO=RHO, RHO_AS=RHO_AS,
                                 tau1=tau1, tau2=tau2, A=A,D=D,
                                 THETA=THETA_cur, tauB=tauB, tau=tau)
  para_cur = keg_lbfgs$par

  if(any(c(keg_lbfgs$convergence,theta_lbfgs$convergence ) != 0))
    message(sprintf('keg and theta error code: %s, %s',
                    keg_lbfgs$convergence, theta_lbfgs$convergence))

  # update parameters
  while(iter < maxiter & (min(abs(para_cur- para)) > tol |
                          min(abs(beta_cur - betas)) > tol |
                          abs(phi_cur - phi) >tol |
                          abs(THETA_cur - THETA) > tol)){
    para = para_cur
    betas = beta_cur
    phi = phi_cur
    THETA = THETA_cur

    if(H0 == 0){
      KAPPA = exp(para_cur[1])
      ETA = exp(para_cur[2])
      GAMMA = exp(para_cur[3])
    } else if(H0 == 1){
      KAPPA = exp(para_cur[1])
      ETA = 1.0;
      GAMMA = exp(para_cur[2])
    } else if(H0 == 2){
      KAPPA = exp(para_cur[1])
      ETA = exp(para_cur[2])
      GAMMA = 1.0
    }

    offsets <- compute_offset(z=z, RHO=RHO ,KAPPA= KAPPA, ETA=ETA,
                              GAMMA=GAMMA, tau1=tau1, tau2=tau2)
    bp = update_beta_phi(y,X, offsets = offsets )
    beta_cur = bp$betas
    phi_cur =  bp$phi

    theta_lbfgs = update_theta(para=para, H0=H0, z_AS=z_AS,
                               RHO_AS=RHO_AS, A=A,D=D, THETA=THETA,
                               tauB=tauB, tau=tau)

    THETA_cur = theta_lbfgs$par

    # update keg
    keg_lbfgs = update_keg_trecase(H0=H0, para=para, y =y, z=z, z_AS=z_AS, X=X,
                                   betas=beta_cur, phi=phi_cur, RHO=RHO,
                                   RHO_AS=RHO_AS, tau1=tau1, tau2=tau2,
                                   A=A,D=D, THETA=THETA_cur, tauB=tauB, tau=tau)
    para_cur = keg_lbfgs$par

    # if(any(c(keg_lbfgs$convergence,theta_lbfgs$convergence ) != 0))
    # message(sprintf('keg and theta error code: %s, %s',
    #                 keg_lbfgs$convergence, theta_lbfgs$convergence))

    iter = iter + 1
  }

  if(H0 == 0){
    KAPPA = exp(para_cur[1])
    ETA = exp(para_cur[2])
    GAMMA = exp(para_cur[3])
  } else if(H0 == 1){
    KAPPA = exp(para_cur[1])
    ETA = 1.0;
    GAMMA = exp(para_cur[2])
  } else if(H0 == 2){
    KAPPA = exp(para_cur[1])
    ETA = exp(para_cur[2])
    GAMMA = 1.0
  }
  return(list(KAPPA = KAPPA, ETA = ETA, GAMMA = GAMMA, phi = phi_cur,
              betas = beta_cur, THETA = THETA_cur,
              loglik =keg_lbfgs$value, iters = iter,
              conv = c(keg_lbfgs$convergence,theta_lbfgs$convergence)))
}

TReCASE_test = function(y,z, z_AS, X, RHO, RHO_AS, tau1, tau2, A, D,
                        tauB, tau, maxiter =500, tol = 1e-7){
  H0 =0
  para = c(0,0,0)
  res1 = TReCASE_sfit(H0, para, y,z, z_AS, X, RHO, RHO_AS, tau1, tau2, A, D,
                      tauB, tau,  maxiter=maxiter, tol=tol)

  H0 =1
  para = c(0,0)
  res1.eta1 = TReCASE_sfit(H0, para, y,z, z_AS, X, RHO, RHO_AS, tau1, tau2,
                           A, D, tauB, tau,  maxiter=maxiter, tol=tol)

  H0=2
  res1.gamma1 = TReCASE_sfit(H0, para, y,z, z_AS, X, RHO, RHO_AS, tau1, tau2,
                             A, D, tauB, tau,  maxiter=maxiter, tol=tol)

  H0=0
  para = c(res1.eta1$KAPPA, res1.gamma1$ETA, res1.eta1$GAMMA)
  res1b = TReCASE_sfit(H0, para, y,z, z_AS, X, RHO, RHO_AS, tau1, tau2,
                       A, D, tauB, tau,  maxiter=maxiter, tol=tol)

  if(res1b$loglik < res1$loglik)
    res1 = res1b
  p.eta = pchisq(2*( res1.eta1$loglik- res1$loglik ), df =1, lower.tail = F)
  p.gamma = pchisq(2*( res1.gamma1$loglik- res1$loglik ), df =1, lower.tail = F)

  return(list(p.eta = p.eta, p.gamma = p.gamma, loglik.full = - res1$loglik,
              loglik.eta = - res1.eta1$loglik,
              loglik.gamma = - res1.gamma1$loglik,
              KAPPA = res1$KAPPA, ETA = res1$ETA, GAMMA = res1$GAMMA,
              betas = res1$betas, phi = res1$phi, THETA = res1$THETA
              #, iters = c(res1$iters, res1.eta1$iters, res1.gamma1$iters)
  ))
}



# ----------------------------------------------------------------------
# Cis-Trans Score test
# ----------------------------------------------------------------------

CisTrans_ScoreObs <- function(para, y,z,z_AS, X, betas, phi,
                              RHO, RHO_AS, tau1, tau2, A, D, tauB, tau, THETA,
                              Power = F){
  KAPPA = exp(para[1])
  ETA   = exp(para[2])
  GAMMA = exp(para[3])

  n = length(RHO)

  # ---------------------------------
  # Information matrix (TReC)
  # ---------------------------------
  offsets = compute_offset(z=z, RHO=RHO ,KAPPA= KAPPA, ETA=ETA,
                           GAMMA=GAMMA, tau1=tau1, tau2=tau2)
  if(is.null(dim(X)))
    X = data.matrix(X)

  expx.t.betas = exp(X%*%betas)
  mu = exp(offsets)*expx.t.betas
  phi_mu_1     = (1 + phi * mu)

  ## Iee
  dlTRec.dmu  = y/mu - (1 + phi * y)/phi_mu_1               #Expectation
  d2lTReC.dmu = - y/mu^2 + phi*(1 + phi * y)/phi_mu_1^2     #Expectation
  RHO_1       = 1 - RHO

  dmu.dKAPPA = dmu.dETA =dmu.dGAMMA = rep(0, n)

  ## z = AA
  indAA = which(z == 1)
  dmu.dKAPPA[indAA] = (tau2[indAA] + tau1[indAA])

  ## z = AB
  indAB = which(z == 2)
  dmu.dKAPPA[indAB] =  tau1[indAB]+tau2[indAB]*GAMMA
  dmu.dETA[indAB]   = 1
  dmu.dGAMMA[indAB] = tau2[indAB]

  ## z = BA
  indBA = which(z == 3)
  dmu.dKAPPA[indBA] = tau2[indBA] + tau1[indBA]*GAMMA
  dmu.dETA[indBA]   = 1
  dmu.dGAMMA[indBA] = tau1[indBA]

  ## z = BB
  indBB = which(z == 4)
  dmu.dKAPPA[indBB] = (tau2[indBB] + tau1[indBB])*GAMMA
  dmu.dETA[indBB]   = 2
  dmu.dGAMMA[indBB] = (tau2[indBB] + tau1[indBB])

  dmu.dKAPPA = expx.t.betas*RHO*dmu.dKAPPA
  dmu.dETA   = expx.t.betas*RHO_1*dmu.dETA
  dmu.dGAMMA = expx.t.betas*RHO*KAPPA*dmu.dGAMMA

  dmu = data.matrix(cbind(dmu.dKAPPA, dmu.dETA, dmu.dGAMMA))

  d2lTReC.dKAPPA = sum(d2lTReC.dmu * dmu.dKAPPA^2)
  d2lTReC.dETA   = sum(d2lTReC.dmu * dmu.dETA^2)
  d2lTReC.dGAMMA = sum(d2lTReC.dmu * dmu.dGAMMA^2)

  d2mu.dKAPPA.dGAMMA = rep(0, n)

  d2mu.dKAPPA.dGAMMA[indAB] = tau2[indAB]
  d2mu.dKAPPA.dGAMMA[indBA] = tau1[indBA]
  d2mu.dKAPPA.dGAMMA[indBB] = (tau1[indBB] + tau2[indBB])
  d2mu.dKAPPA.dGAMMA        = RHO*d2mu.dKAPPA.dGAMMA*expx.t.betas

  Iee = diag(c(d2lTReC.dKAPPA, d2lTReC.dETA, d2lTReC.dGAMMA))
  Iee[1,2] = Iee[2,1] = sum(d2lTReC.dmu*dmu.dKAPPA*dmu.dETA)
  Iee[1,3] = Iee[3,1] = sum(d2lTReC.dmu*dmu.dKAPPA*dmu.dGAMMA +
                              dlTRec.dmu*d2mu.dKAPPA.dGAMMA)
  Iee[2,3] = Iee[3,2] = sum(d2lTReC.dmu*dmu.dETA*dmu.dGAMMA)

  ## check eigen(Iee) is negtive definite
  ## check grad is the correct

  ##Ipp
  vphi = 1/phi
  Ipp_tmp1 = digamma(y+vphi)-digamma(vphi)-log(phi_mu_1)
  Ipp  = sum(2*vphi^3*Ipp_tmp1 + vphi^4*(trigamma(y+vphi)-trigamma(vphi)) +
               2*vphi^2*mu/phi_mu_1 + (vphi+y)*mu^2/phi_mu_1^2 -y*vphi^2)


  delta1 = diag(as.numeric(1/phi_mu_1))
  delta2 = diag(as.numeric(mu/phi_mu_1))
  delta3 = diag(as.numeric(mu*(y-mu)/phi_mu_1^2))
  delta4 = diag(as.numeric((y-mu)/phi_mu_1^2))
  # delta5 = diag(as.numeric(1/(mu*phi_mu_1)))
  # delta6 = diag(as.numeric((y-mu)*(1+2*phi*mu)/(mu^2*phi_mu_1^2 )) )

  ## Ibb
  Ibb = -(t(X)%*%delta2%*%X + phi*t(X)%*%delta3%*%X)

  ##Ibe
  Ibe = - (t(X)%*%delta1%*%dmu + phi*t(X)%*%delta4%*%dmu)

  ##Ibp
  Ibp = - (t(X) %*% ((y-mu)*mu/(phi_mu_1)^2))

  ##Iep
  Iep = - (t(dmu) %*% ((y -mu)/phi_mu_1^2))

  # ---------------------------------
  # Score and Information (ASE)
  # RHO_AS, z_AS
  # ---------------------------------
  indhet  = which(z_AS %in% c(2,3))
  indhomo = which(z_AS %in% c(1,4))
  vTHETA  = 1/THETA
  n_AS    = length(RHO_AS)

  if(Power){
    ## set ETA = GAMMA = 1 in pis to test the power of score test
    pis = compute_pi(z_AS, RHO_AS, KAPPA, 1, 1, tauB, tau)
  }else{
    pis = compute_pi(z_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau)
  }

  # score alpha_GAMMA = alpha_ETA
  ## dlASE.dpi
  pis_1 = 1 - pis
  dlASE.dpi = vTHETA*(digamma(A+pis*vTHETA) + digamma(pis_1*vTHETA) -
                        digamma(D-A+pis_1*vTHETA) - digamma(pis*vTHETA)) ## check

  ## AA,BB
  dpi.dKAPPA = dpi.dETA = dpi.dGAMMA = rep(0, n_AS)
  tmp.homo = (RHO_AS[indhomo]*tau[indhomo]*KAPPA + 2*(1-RHO_AS[indhomo]))
  dpi.dKAPPA[indhomo] = RHO_AS[indhomo]*tauB[indhomo]/tmp.homo -
    RHO_AS[indhomo]*tau[indhomo]*(RHO_AS[indhomo]*tauB[indhomo]*
                                    KAPPA+(1-RHO_AS[indhomo]))/tmp.homo^2

  ## AB, BA
  tmp1 = RHO_AS[indhet]*tauB[indhet]*GAMMA*KAPPA+(1-RHO_AS[indhet])*ETA
  tmp2 = RHO_AS[indhet]*(tau[indhet]-tauB[indhet])*KAPPA +
    (1 - RHO_AS[indhet]) + tmp1
  tmp3 = tmp1/tmp2^2

  dpi.dKAPPA[indhet] = RHO_AS[indhet]*tauB[indhet]*GAMMA/tmp2 -
    (RHO_AS[indhet]*(tau[indhet]-tauB[indhet])+
       RHO_AS[indhet]*tauB[indhet]*GAMMA)*tmp3
  dpi.dETA[indhet] = (1-RHO_AS[indhet])/tmp2 - (1-RHO_AS[indhet])*tmp3
  dpi.dGAMMA[indhet] =RHO_AS[indhet]*tauB[indhet]*KAPPA/tmp2 -
    (RHO_AS[indhet]*tauB[indhet]*KAPPA)*tmp3

  dlASE.dKAPPA = sum(dlASE.dpi * dpi.dKAPPA)
  dlASE.dETA = sum(dlASE.dpi * dpi.dETA)
  dlASE.dGAMMA = sum(dlASE.dpi * dpi.dGAMMA )
  ## consistant with gradBB

  score_alpha = c(dlASE.dETA, dlASE.dGAMMA)

  ## Iee
  d2lASE.dpi = vTHETA^2*(trigamma(vTHETA*pis+A) + trigamma(vTHETA*pis_1+D-A)-
                           trigamma(vTHETA*pis) - trigamma(vTHETA*pis_1)) ##check

  ## d2pi.dlambda
  # AB or BA
  # tmp1 = RHO_AS[indhet]*tauB[indhet]*GAMMA*KAPPA+(1-RHO_AS[indhet])*ETA
  # tmp2 = RHO_AS[indhet]*(tau[indhet]-tauB[indhet])*KAPPA + (1-RHO_AS[indhet])+
  #   tmp1

  d2pi.dKAPPA=d2pi.dGAMMA=d2pi.dETA=rep(0,n_AS)
  tmp.indhet.Kappa =(RHO_AS[indhet]*(tau[indhet]-tauB[indhet])+
                       RHO_AS[indhet]*tauB[indhet]*GAMMA)

  d2pi.dKAPPA[indhet] = -2*GAMMA*RHO_AS[indhet]*tauB[indhet]*
    tmp.indhet.Kappa/tmp2^2 + 2*tmp.indhet.Kappa^2*tmp1/tmp2^3 ## done check
  d2pi.dETA[indhet] = - 2*(1-RHO[indhet])^2/tmp2^2 +
    2*(1-RHO[indhet])^2*tmp1/tmp2^3 ## done check
  d2pi.dGAMMA[indhet] = -2*(RHO[indhet]*tauB[indhet]*KAPPA)^2/tmp2^2 +
    2*(RHO[indhet]*tauB[indhet]*KAPPA)^2*tmp1/tmp2^3 ## done check

  ## AA or BB
  d2pi.dKAPPA[indhomo] = -(2*RHO_AS[indhomo]^2*(1-RHO_AS[indhomo])*
                             (2*tauB[indhomo]-tau[indhomo])*tau[indhomo])/tmp.homo^3
  # d2pi.dKAPPA[indhomo] = -(2*RHO_AS[indhomo]^2*tauB[indhomo]*tau[indhomo])/
  #   tmp.homo^2 + 2*RHO_AS[indhomo]^2*tau[indhomo]^2*
  #   (1-RHO_AS[indhomo]+KAPPA*RHO_AS[indhomo]*tauB[indhomo])/tmp.homo^3 ## done check

  d2lASE.dKAPPA = sum(d2lASE.dpi*dpi.dKAPPA^2 + dlASE.dpi*d2pi.dKAPPA)
  d2lASE.dETA = sum(d2lASE.dpi*dpi.dETA^2 + dlASE.dpi*d2pi.dETA)
  d2lASE.dGAMMA = sum(d2lASE.dpi*dpi.dGAMMA^2 + dlASE.dpi*d2pi.dGAMMA)

  ## d2pi.dlambda1.dlambda2
  d2pi.dKAPPA.dETA = d2pi.dKAPPA.dGAMMA=d2pi.dETA.dGAMMA= rep(0, n_AS)
  d2pi.dKAPPA.dETA[indhet]   = -RHO_AS[indhet]*tauB[indhet]*GAMMA*(1-RHO[indhet])/tmp2^2-
    tmp.indhet.Kappa*(1-RHO[indhet])/tmp2^2 +
    2*(1-RHO[indhet])*tmp1*tmp.indhet.Kappa/tmp2^3 ## done check
  d2pi.dKAPPA.dGAMMA[indhet] = RHO_AS[indhet]*tauB[indhet]/tmp2 -
    (RHO_AS[indhet]*tauB[indhet])^2*KAPPA*GAMMA/tmp2^2 -
    tmp.indhet.Kappa*RHO_AS[indhet]*tauB[indhet]*KAPPA/tmp2^2 -
    RHO_AS[indhet]*tauB[indhet]*tmp1/tmp2^2 +
    2*tmp.indhet.Kappa*RHO_AS[indhet]*tauB[indhet]*KAPPA*tmp1/tmp2^3 ## done check
  d2pi.dETA.dGAMMA[indhet]   = -2*KAPPA*(1-RHO_AS[indhet])*RHO_AS[indhet]*tauB[indhet]/tmp2^2 +
    2*KAPPA*(1-RHO_AS[indhet])*RHO_AS[indhet]*tauB[indhet]*tmp1/tmp2^3 ## done check

  d2lASE.dKAPPA.dETA = sum(d2lASE.dpi*dpi.dKAPPA*dpi.dETA +
                             dlASE.dpi*d2pi.dKAPPA.dETA)
  d2lASE.dKAPPA.dGAMMA = sum(d2lASE.dpi*dpi.dKAPPA*dpi.dGAMMA +
                               dlASE.dpi*d2pi.dKAPPA.dGAMMA)
  d2lASE.dETA.dGAMMA =sum(d2lASE.dpi*dpi.dETA*dpi.dGAMMA +
                            dlASE.dpi*d2pi.dETA.dGAMMA )

  ## Iee2
  Iee[1,1] = Iee[1,1] + d2lASE.dKAPPA
  Iee[2,2] = Iee[2,2] + d2lASE.dETA
  Iee[3,3] = Iee[3,3] + d2lASE.dGAMMA
  Iee[1,2] = Iee[2,1] = Iee[2,1] + d2lASE.dKAPPA.dETA
  Iee[1,3] = Iee[3,1] = Iee[3,1] + d2lASE.dKAPPA.dGAMMA
  Iee[2,3] = Iee[3,2] = Iee[2,3] + d2lASE.dETA.dGAMMA

  ## Iaa alpha, alpha
  Iaa = matrix(0, 2,2)
  Iaa[1,1] = d2lASE.dETA
  Iaa[2,2] = d2lASE.dGAMMA
  Iaa[1,2] = Iaa[2,1] = d2lASE.dETA.dGAMMA

  ##Iae alpha_ETA, alpha_GAMMA * KAPPA, ETA, GAMMA
  Iae = matrix(0, 2,3)
  Iae[1,1] = d2lASE.dKAPPA.dETA
  Iae[1,2] = d2lASE.dETA
  Iae[1,3] = Iae[2,2] = d2lASE.dETA.dGAMMA
  Iae[2,1] = d2lASE.dKAPPA.dGAMMA
  Iae[2,3] = d2lASE.dGAMMA

  ##Itt theta
  W = -(digamma(D+vTHETA)+digamma(pis*vTHETA)*pis+digamma(pis_1*vTHETA)*pis_1 -
          digamma(A+pis*vTHETA)*pis-digamma(D-A+pis_1*vTHETA)*pis_1-
          digamma(vTHETA)) ##check
  dW.dTHETA = -vTHETA^2*(pis^2*(trigamma(vTHETA*pis+A) - trigamma(pis*vTHETA)) +
                           pis_1^2*(trigamma(vTHETA*pis_1 + D-A) -
                                      trigamma(vTHETA*pis_1)) +
                           trigamma(vTHETA) - trigamma(vTHETA+D))
  dlASE.dTHETA = -sum(vTHETA^2*W)
  Itt = sum(2*vTHETA^3*W - vTHETA^2*dW.dTHETA)

  ##Iat  alpha, theta (2,1) same as ETA,GAMMA x Theta
  d2lASE.dTHETA.dpi = -vTHETA^2*(digamma(vTHETA*pis+A) - digamma(vTHETA*pis) -
                                   digamma(vTHETA*pis_1+D-A) +
                                   digamma(vTHETA*pis_1)) -
    vTHETA^3*pis*(trigamma(vTHETA*pis+A) - trigamma(vTHETA*pis)) +
    vTHETA^3*pis_1*(trigamma(vTHETA*pis_1+D-A) - trigamma(vTHETA*pis_1)) ##mathmetica

  d2lASE.dTHETA.dKAPPA = sum(d2lASE.dTHETA.dpi*dpi.dKAPPA)
  d2lASE.dTHETA.dETA = sum(d2lASE.dTHETA.dpi*dpi.dETA)
  d2lASE.dTHETA.dGAMMA = sum(d2lASE.dTHETA.dpi*dpi.dGAMMA)

  Iat = rbind(d2lASE.dTHETA.dETA, d2lASE.dTHETA.dGAMMA)
  ##Iet KEG, theta
  Iet = rbind(d2lASE.dTHETA.dKAPPA, Iat)

  # ---------------------------------
  # Final Information matrix
  # ---------------------------------
  np = length(betas)
  M1 = matrix(0, np+5, np+5)
  M1[1:np, 1:np] = Ibb
  M1[1:np, (np+1):(np+3)] = Ibe
  M1[(np+1):(np+3),1:np] = t(Ibe)
  M1[(np+1):(np+3),(np+1):(np+3)] = Iee
  M1[(np+1):(np+3),(np+4)] = Iep
  M1[(np+4),(np+1):(np+3)] = t(Iep)
  M1[1:np, np+4] = Ibp
  M1[np+4, 1:np] = t(Ibp)
  M1[(np+4),(np+4)] = Ipp
  M1[(np+5),(np+5)] = Itt
  M1[(np+5),(np+1):(np+3)] = Iet
  M1[(np+1):(np+3),(np+5)] = t(Iet)
  M2 = matrix(0, np+5,2)
  M2[1:3+np,]  = t(Iae)
  M2[5+np, ] = t(Iat)
  M1 = -M1
  M2 = -M2
  Iaa = -Iaa

  OImat = matrix(0, np+7, np+7)
  OImat[1:(np+5),1:(np+5)] = M1
  OImat[(np+6):(np+7),(np+6):(np+7)] = Iaa
  OImat[1:(np+5),(np+6):(np+7)] = M2
  OImat[(np+6):(np+7),1:(np+5)] = t(M2)

  M2M1 = t(M2)%*%solve(M1)%*%M2
  score_alpha = data.matrix(score_alpha)
  Score = t(score_alpha) %*% solve(Iaa - M2M1) %*% (score_alpha)
  p.score = pchisq(Score, 2, lower.tail = F)

  if(any(eigen(OImat)$values<0)){
    return(c(Score = -1, p.score = NA)) # Info matrix is not positive definite
  }

  return(c(Score = Score, p.score = p.score))
}


# outvec[1] = Expected value of digamma(vTHETA*pi+A)
# outvec[2] = Expected value of digamma(vTHETA*(1.0-pi)+D-A)
# outvec[3] = Expected value of trigamma(vTHETA*pi+A)
# outvec[4] = Expected value of trigamma(vTHETA*(1.0-pi)+D-A)


ASE_ExpFunc<-function(Di, pis.i, vTHETA){
  THETA = 1/vTHETA
  lpmf = sapply(0:Di, loglikBB_THETA, THETA=THETA, D = Di, pis=pis.i)
  pmfvec = exp(-lpmf) ##
  outvec = vector('numeric', 4)
  outvec[1] = sum(digamma(vTHETA*pis.i + (0:Di)) * pmfvec)
  outvec[2] = sum(digamma(vTHETA*(1-pis.i) + (Di:0)) * pmfvec)
  outvec[3] = sum(trigamma(vTHETA*pis.i+(0:Di)) * pmfvec)
  outvec[4] = sum(trigamma(vTHETA*(1-pis.i)+(Di:0)) * pmfvec)
  return(outvec)
}

CisTrans_Score <- function(para, y,z,z_AS, X, betas, phi,
                           RHO, RHO_AS, tau1, tau2, A, D, tauB, tau, THETA,
                           Power = F){
  KAPPA = exp(para[1])
  ETA   = exp(para[2])
  GAMMA = exp(para[3])

  n = length(RHO)

  # ---------------------------------
  # Information matrix (TReC)
  # ---------------------------------
  offsets = compute_offset(z=z, RHO=RHO ,KAPPA= KAPPA, ETA=ETA,
                           GAMMA=GAMMA, tau1=tau1, tau2=tau2)
  if(is.null(dim(X)))
    X = data.matrix(X)

  expx.t.betas = exp(X%*%betas)
  mu = exp(offsets)*expx.t.betas
  phi_mu_1     = (1 + phi * mu)

  ## Iee
  # dlTRec.dmu  = 0
  d2lTReC.dmu = - 1/(mu*phi_mu_1)
  RHO_1       = 1 - RHO

  dmu.dKAPPA = dmu.dETA =dmu.dGAMMA = rep(0, n)

  ## z = AA
  indAA = which(z == 1)
  dmu.dKAPPA[indAA] = (tau2[indAA] + tau1[indAA])

  ## z = AB
  indAB = which(z == 2)
  dmu.dKAPPA[indAB] =  tau1[indAB]+tau2[indAB]*GAMMA
  dmu.dETA[indAB]   = 1
  dmu.dGAMMA[indAB] = tau2[indAB]

  ## z = BA
  indBA = which(z == 3)
  dmu.dKAPPA[indBA] = tau2[indBA] + tau1[indBA]*GAMMA
  dmu.dETA[indBA]   = 1
  dmu.dGAMMA[indBA] = tau1[indBA]

  ## z = BB
  indBB = which(z == 4)
  dmu.dKAPPA[indBB] = (tau2[indBB] + tau1[indBB])*GAMMA
  dmu.dETA[indBB]   = 2
  dmu.dGAMMA[indBB] = (tau2[indBB] + tau1[indBB])

  dmu.dKAPPA = expx.t.betas*RHO*dmu.dKAPPA
  dmu.dETA   = expx.t.betas*RHO_1*dmu.dETA
  dmu.dGAMMA = expx.t.betas*RHO*KAPPA*dmu.dGAMMA

  dmu = data.matrix(cbind(dmu.dKAPPA, dmu.dETA, dmu.dGAMMA))

  d2lTReC.dKAPPA = sum(d2lTReC.dmu * dmu.dKAPPA^2)
  d2lTReC.dETA   = sum(d2lTReC.dmu * dmu.dETA^2)
  d2lTReC.dGAMMA = sum(d2lTReC.dmu * dmu.dGAMMA^2)

  d2mu.dKAPPA.dGAMMA = rep(0, n)

  d2mu.dKAPPA.dGAMMA[indAB] = tau2[indAB]
  d2mu.dKAPPA.dGAMMA[indBA] = tau1[indBA]
  d2mu.dKAPPA.dGAMMA[indBB] = (tau1[indBB] + tau2[indBB])
  d2mu.dKAPPA.dGAMMA        = RHO*d2mu.dKAPPA.dGAMMA*expx.t.betas

  Iee = diag(c(d2lTReC.dKAPPA, d2lTReC.dETA, d2lTReC.dGAMMA))
  Iee[1,2] = Iee[2,1] = sum(d2lTReC.dmu*dmu.dKAPPA*dmu.dETA)
  Iee[1,3] = Iee[3,1] = sum(d2lTReC.dmu*dmu.dKAPPA*dmu.dGAMMA)
  Iee[2,3] = Iee[3,2] = sum(d2lTReC.dmu*dmu.dETA*dmu.dGAMMA)

  ## check eigen(Iee) is negtive definite
  ## Iee similar to Obs

  delta1 = diag(as.numeric(1/phi_mu_1))
  delta2 = diag(as.numeric(mu/phi_mu_1))

  ## Ibb
  Ibb = -(t(X)%*%delta2%*%X)

  ##Ibe
  Ibe = - (t(X)%*%delta1%*%dmu)

  ## Ibb, Ibe similar to Obs
  # ---------------------------------
  # Score and Information (ASE)
  # RHO_AS, z_AS
  # ---------------------------------
  indhet  = which(z_AS %in% c(2,3))
  indhomo = which(z_AS %in% c(1,4))
  vTHETA  = 1/THETA
  n_AS    = length(RHO_AS)

  if(Power){
    ## set ETA = GAMMA = 1 in pis to test the power of score test
    pis = compute_pi(z_AS, RHO_AS, KAPPA, 1, 1, tauB, tau)
  }else{
    pis = compute_pi(z_AS, RHO_AS, KAPPA, ETA, GAMMA, tauB, tau)
  }

  # score alpha_GAMMA = alpha_ETA
  Expmat = mapply(ASE_ExpFunc, D, pis, vTHETA)
  # outvec[1] = Expected value of digamma(vTHETA*pi+A)
  # outvec[2] = Expected value of digamma(vTHETA*(1.0-pi)+D-A)
  # outvec[3] = Expected value of trigamma(vTHETA*pi+A)
  # outvec[4] = Expected value of trigamma(vTHETA*(1.0-pi)+D-A)

  ## dlASE.dpi
  pis_1 = 1 - pis
  dlASE.dpi = vTHETA*(digamma(A+pis*vTHETA) - digamma(D-A+pis_1*vTHETA) +
                        digamma(pis_1*vTHETA) - digamma(pis*vTHETA))

  ## AA,BB
  dpi.dKAPPA = dpi.dETA = dpi.dGAMMA = rep(0, n_AS)
  tmp.homo = (RHO_AS[indhomo]*tau[indhomo]*KAPPA + 2*(1-RHO_AS[indhomo]))
  dpi.dKAPPA[indhomo] = RHO_AS[indhomo]*tauB[indhomo]/tmp.homo -
    RHO_AS[indhomo]*tau[indhomo]*(RHO_AS[indhomo]*tauB[indhomo]*
                                    KAPPA+(1-RHO_AS[indhomo]))/tmp.homo^2

  ## AB, BA
  tmp1 = RHO_AS[indhet]*tauB[indhet]*GAMMA*KAPPA+(1-RHO_AS[indhet])*ETA
  tmp2 = RHO_AS[indhet]*(tau[indhet]-tauB[indhet])*KAPPA +
    (1 - RHO_AS[indhet]) + tmp1
  tmp3 = tmp1/tmp2^2

  dpi.dKAPPA[indhet] = RHO_AS[indhet]*tauB[indhet]*GAMMA/tmp2 -
    (RHO_AS[indhet]*(tau[indhet]-tauB[indhet])+
       RHO_AS[indhet]*tauB[indhet]*GAMMA)*tmp3
  dpi.dETA[indhet] = (1-RHO_AS[indhet])/tmp2 - (1-RHO_AS[indhet])*tmp3
  dpi.dGAMMA[indhet] =RHO_AS[indhet]*tauB[indhet]*KAPPA/tmp2 -
    (RHO_AS[indhet]*tauB[indhet]*KAPPA)*tmp3

  dlASE.dKAPPA = sum(dlASE.dpi * dpi.dKAPPA)
  dlASE.dETA = sum(dlASE.dpi * dpi.dETA)
  dlASE.dGAMMA = sum(dlASE.dpi * dpi.dGAMMA )


  score_alpha = c(dlASE.dETA, dlASE.dGAMMA)

  ## Iee
  d2lASE.dpi = vTHETA^2*(Expmat[3,] + Expmat[4,]-
                           trigamma(vTHETA*pis) - trigamma(vTHETA*pis_1))

  ## d2pi.dlambda
  # AB or BA
  # tmp1 = RHO_AS[indhet]*tauB[indhet]*GAMMA*KAPPA+(1-RHO_AS[indhet])*ETA
  # tmp2 = RHO_AS[indhet]*(tau[indhet]-tauB[indhet])*KAPPA + (1-RHO_AS[indhet])+
  #   tmp1

  d2pi.dKAPPA=d2pi.dGAMMA=d2pi.dETA=rep(0,n_AS)
  tmp.indhet.Kappa =(RHO_AS[indhet]*(tau[indhet]-tauB[indhet])+
                       RHO_AS[indhet]*tauB[indhet]*GAMMA)

  d2pi.dKAPPA[indhet] = -2*GAMMA*RHO_AS[indhet]*tauB[indhet]*
    tmp.indhet.Kappa/tmp2^2 + 2*tmp.indhet.Kappa^2*tmp1/tmp2^3 ## done check
  d2pi.dETA[indhet] = - 2*(1-RHO[indhet])^2/tmp2^2 +
    2*(1-RHO[indhet])^2*tmp1/tmp2^3 ## done check
  d2pi.dGAMMA[indhet] = -2*(RHO[indhet]*tauB[indhet]*KAPPA)^2/tmp2^2 +
    2*(RHO[indhet]*tauB[indhet]*KAPPA)^2*tmp1/tmp2^3 ## done check

  ## AA or BB
  d2pi.dKAPPA[indhomo] = -(2*RHO_AS[indhomo]^2*(1-RHO_AS[indhomo])*
                             (2*tauB[indhomo]-tau[indhomo])*tau[indhomo])/tmp.homo^3
  # d2pi.dKAPPA[indhomo] = -(2*RHO_AS[indhomo]^2*tauB[indhomo]*tau[indhomo])/
  #   tmp.homo^2 + 2*RHO_AS[indhomo]^2*tau[indhomo]^2*
  #   (1-RHO_AS[indhomo]+KAPPA*RHO_AS[indhomo]*tauB[indhomo])/tmp.homo^3 ## done check

  d2lASE.dKAPPA = sum(d2lASE.dpi*dpi.dKAPPA^2 + dlASE.dpi*d2pi.dKAPPA)
  d2lASE.dETA = sum(d2lASE.dpi*dpi.dETA^2 + dlASE.dpi*d2pi.dETA)
  d2lASE.dGAMMA = sum(d2lASE.dpi*dpi.dGAMMA^2 + dlASE.dpi*d2pi.dGAMMA)

  ## d2pi.dlambda1.dlambda2
  d2pi.dKAPPA.dETA = d2pi.dKAPPA.dGAMMA=d2pi.dETA.dGAMMA= rep(0, n_AS)
  d2pi.dKAPPA.dETA[indhet]   = -RHO_AS[indhet]*tauB[indhet]*GAMMA*(1-RHO[indhet])/tmp2^2-
    tmp.indhet.Kappa*(1-RHO[indhet])/tmp2^2 +
    2*(1-RHO[indhet])*tmp1*tmp.indhet.Kappa/tmp2^3 ## done check
  d2pi.dKAPPA.dGAMMA[indhet] = RHO_AS[indhet]*tauB[indhet]/tmp2 -
    (RHO_AS[indhet]*tauB[indhet])^2*KAPPA*GAMMA/tmp2^2 -
    tmp.indhet.Kappa*RHO_AS[indhet]*tauB[indhet]*KAPPA/tmp2^2 -
    RHO_AS[indhet]*tauB[indhet]*tmp1/tmp2^2 +
    2*tmp.indhet.Kappa*RHO_AS[indhet]*tauB[indhet]*KAPPA*tmp1/tmp2^3 ## done check
  d2pi.dETA.dGAMMA[indhet]   = -2*KAPPA*(1-RHO_AS[indhet])*RHO_AS[indhet]*tauB[indhet]/tmp2^2 +
    2*KAPPA*(1-RHO_AS[indhet])*RHO_AS[indhet]*tauB[indhet]*tmp1/tmp2^3 ## done check

  d2lASE.dKAPPA.dETA = sum(d2lASE.dpi*dpi.dKAPPA*dpi.dETA +
                             dlASE.dpi*d2pi.dKAPPA.dETA)
  d2lASE.dKAPPA.dGAMMA = sum(d2lASE.dpi*dpi.dKAPPA*dpi.dGAMMA +
                               dlASE.dpi*d2pi.dKAPPA.dGAMMA)
  d2lASE.dETA.dGAMMA =sum(d2lASE.dpi*dpi.dETA*dpi.dGAMMA +
                            dlASE.dpi*d2pi.dETA.dGAMMA )

  ## Iee2
  Iee[1,1] = Iee[1,1] + d2lASE.dKAPPA
  Iee[2,2] = Iee[2,2] + d2lASE.dETA
  Iee[3,3] = Iee[3,3] + d2lASE.dGAMMA
  Iee[1,2] = Iee[2,1] = Iee[2,1] + d2lASE.dKAPPA.dETA
  Iee[1,3] = Iee[3,1] = Iee[3,1] + d2lASE.dKAPPA.dGAMMA
  Iee[2,3] = Iee[3,2] = Iee[2,3] + d2lASE.dETA.dGAMMA
  ## Iee similar to Obs

  ## Iaa alpha, alpha
  Iaa = matrix(0, 2,2)
  Iaa[1,1] = d2lASE.dETA
  Iaa[2,2] = d2lASE.dGAMMA
  Iaa[1,2] = Iaa[2,1] = d2lASE.dETA.dGAMMA
  ## Iaa similar to Obs

  ##Iae alpha_ETA, alpha_GAMMA * KAPPA, ETA, GAMMA
  Iae = matrix(0, 2,3)
  Iae[1,1] = d2lASE.dKAPPA.dETA
  Iae[1,2] = d2lASE.dETA
  Iae[1,3] = Iae[2,2] = d2lASE.dETA.dGAMMA
  Iae[2,1] = d2lASE.dKAPPA.dGAMMA
  Iae[2,3] = d2lASE.dGAMMA
  ## Iae[1,1] -1.483045, 0.7432159[obs]

  ##Itt theta
  W = -(digamma(D+vTHETA)+digamma(pis*vTHETA)*pis+digamma(pis_1*vTHETA)*pis_1 -
          Expmat[1,]*pis-Expmat[2,]*pis_1-
          digamma(vTHETA))
  dW.dTHETA = -vTHETA^2*(pis^2*(Expmat[3,] - trigamma(pis*vTHETA)) +
                           pis_1^2*(Expmat[4,] -
                                      trigamma(vTHETA*pis_1)) +
                           trigamma(vTHETA) - trigamma(vTHETA+D)) ## similar to obs
  dlASE.dTHETA = -sum(vTHETA^2*W) ##??
  Itt = sum(2*vTHETA^3*W - vTHETA^2*dW.dTHETA) ## -889.6219/ -870.7088(obs)

  ##Iat  alpha, theta (2,1) same as ETA,GAMMA x Theta
  d2lASE.dTHETA.dpi = -vTHETA^2*(Expmat[1,] - digamma(vTHETA*pis) -
                                   Expmat[2,] +
                                   digamma(vTHETA*pis_1)) -
    vTHETA^3*pis*(Expmat[3,] - trigamma(vTHETA*pis)) +
    vTHETA^3*pis_1*(Expmat[4,]- trigamma(vTHETA*pis_1)) ##mathmetica

  d2lASE.dTHETA.dKAPPA = sum(d2lASE.dTHETA.dpi*dpi.dKAPPA)
  d2lASE.dTHETA.dETA = sum(d2lASE.dTHETA.dpi*dpi.dETA)
  d2lASE.dTHETA.dGAMMA = sum(d2lASE.dTHETA.dpi*dpi.dGAMMA)## -5.373169/8.379897(obs)

  Iat = rbind(d2lASE.dTHETA.dETA, d2lASE.dTHETA.dGAMMA)

  ##Iet KEG, theta
  Iet = rbind(d2lASE.dTHETA.dKAPPA, Iat)

  # ---------------------------------
  # Final Information matrix
  # ---------------------------------
  np = length(betas)
  M1 = matrix(0, np+4, np+4)
  M1[1:np, 1:np] = Ibb
  M1[1:np, (np+1):(np+3)] = Ibe
  M1[(np+1):(np+3),1:np] = t(Ibe)
  M1[(np+1):(np+3),(np+1):(np+3)] = Iee
  M1[(np+4),(np+4)] = Itt
  M1[(np+4),(np+1):(np+3)] = Iet
  M1[(np+1):(np+3),(np+4)] = t(Iet)
  M2 = matrix(0, np+4,2)
  M2[1:3+np,]  = t(Iae)
  M2[4+np, ] = t(Iat)
  M1 = -M1
  M2 = -M2
  Iaa = -Iaa

  OImat = matrix(0, np+6, np+6)
  OImat[1:(np+4),1:(np+4)] = M1
  OImat[(np+5):(np+6),(np+5):(np+6)] = Iaa
  OImat[1:(np+4),(np+5):(np+6)] = M2
  OImat[(np+5):(np+6),1:(np+4)] = t(M2)

  M2M1 = t(M2)%*%solve(M1)%*%M2
  score_alpha = data.matrix(score_alpha)
  Score = t(score_alpha) %*% solve(Iaa - M2M1) %*% (score_alpha)
  p.score = pchisq(Score, 2, lower.tail = F)

  if(any(eigen(OImat)$values<0)){
    return(c(Score = -1, p.score = p.score)) # Info matrix is not positive definite
  }

  return(c(Score = Score, p.score = p.score))
}



# ----------------------------------------------------------------------
# ASE
# ----------------------------------------------------------------------


# updata kappa gamma eta in ASE
update_keg_ase <- function(para, H0, z_AS, RHO_AS, A, D, THETA, tauB, tau){
  if(H0 == 0){
    lbfgs = optim(par = para, fn = loglikBB, gr = Grad_BB_keg, H0 = 0,
                  z_AS=z_AS, RHO_AS=RHO_AS, A=A, D=D, THETA=THETA,
                  tauB=tauB, tau=tau,
                  lower = -50, upper = 10,method = "L-BFGS-B")
  }
  if(H0 == 1){
    lbfgs = optim(par = para, fn = loglikBB, gr = Grad_BB_keg, H0 = 1,
                  z_AS=z_AS, RHO_AS=RHO_AS, A=A, D=D, THETA=THETA,
                  tauB=tauB, tau=tau,
                  lower = -50, upper = 10,method = "L-BFGS-B")
  }
  if(H0 == 2){
    lbfgs = optim(par = para, fn = loglikBB, gr = Grad_BB_keg, H0 = 2,
                  z_AS=z_AS, RHO_AS=RHO_AS, A=A, D=D, THETA=THETA,
                  tauB=tauB, tau=tau,
                  lower = -50, upper = 10,method = "L-BFGS-B")
  }
  return(lbfgs)
}


ASE_sfit = function(para, H0, z_AS, RHO_AS, A, D, THETA, tauB, tau,
                    maxiter =500, tol = 1e-7){
  iter = 0
  #1 initialize all the parameters
  THETA = 0.05
  if(H0 == 0){
    KAPPA = exp(para[1])
    ETA = exp(para[2])
    GAMMA = exp(para[3])
  } else if(H0 == 1){
    KAPPA = exp(para[1])
    ETA = 1.0;
    GAMMA = exp(para[2])
  } else if(H0 == 2){
    KAPPA = exp(para[1])
    ETA = exp(para[2])
    GAMMA = 1.0
  }

  # update theta
  theta_lbfgs = update_theta(para=para, H0=H0, z_AS=z_AS, RHO_AS=RHO_AS,
                             A=A, D=D, THETA=THETA, tauB=tauB, tau=tau)

  THETA_cur = theta_lbfgs$par

  # update keg
  keg_lbfgs = update_keg_ase(para=para, H0=H0, z_AS=z_AS, RHO_AS=RHO_AS,
                             A=A, D=D, THETA=THETA, tauB=tauB, tau=tau)
  para_cur = keg_lbfgs$par

  if(any(c(keg_lbfgs$convergence,theta_lbfgs$convergence ) != 0))
    message(sprintf('keg and theta error code: %s, %s',
                    keg_lbfgs$convergence, theta_lbfgs$convergence))

  # update parameters
  while(iter < maxiter & (min(abs(para_cur- para)) > tol |
                          abs(THETA_cur - THETA) > tol)){
    para = para_cur
    THETA = THETA_cur

    if(H0 == 0){
      KAPPA = exp(para_cur[1])
      ETA = exp(para_cur[2])
      GAMMA = exp(para_cur[3])
    } else if(H0 == 1){
      KAPPA = exp(para_cur[1])
      ETA = 1.0;
      GAMMA = exp(para_cur[2])
    } else if(H0 == 2){
      KAPPA = exp(para_cur[1])
      ETA = exp(para_cur[2])
      GAMMA = 1.0
    }

    theta_lbfgs = update_theta(para=para, H0=H0, z_AS=z_AS,
                               RHO_AS=RHO_AS, A=A,D=D, THETA=THETA,
                               tauB=tauB, tau=tau)

    THETA_cur = theta_lbfgs$par

    # update keg
    keg_lbfgs = update_keg_ase(para=para, H0=H0, z_AS=z_AS, RHO_AS=RHO_AS,
                               A=A, D=D, THETA=THETA, tauB=tauB, tau=tau)
    para_cur = keg_lbfgs$par

    # if(any(c(keg_lbfgs$convergence,theta_lbfgs$convergence ) != 0))
    # message(sprintf('keg and theta error code: %s, %s',
    #                 keg_lbfgs$convergence, theta_lbfgs$convergence))

    iter = iter + 1
  }

  if(H0 == 0){
    KAPPA = exp(para_cur[1])
    ETA = exp(para_cur[2])
    GAMMA = exp(para_cur[3])
  } else if(H0 == 1){
    KAPPA = exp(para_cur[1])
    ETA = 1.0;
    GAMMA = exp(para_cur[2])
  } else if(H0 == 2){
    KAPPA = exp(para_cur[1])
    ETA = exp(para_cur[2])
    GAMMA = 1.0
  }
  return(list(KAPPA = KAPPA, ETA = ETA, GAMMA = GAMMA, THETA = THETA_cur,
              loglik =keg_lbfgs$value, iters = iter,
              conv = c(keg_lbfgs$convergence, theta_lbfgs$convergence)))
}

ASE_test = function(z_AS, RHO_AS, A, D, THETA, tauB, tau,
                    maxiter =500, tol = 1e-7){
  H0 =0
  para = c(0,0,0)
  res1 = try(ASE_sfit(para, H0, z_AS, RHO_AS, A, D, THETA, tauB, tau,
                      maxiter=maxiter, tol=tol))

  H0 =1
  para = c(0,0)
  res1.eta1 = ASE_sfit(para, H0, z_AS, RHO_AS, A, D, THETA, tauB, tau,
                       maxiter=maxiter, tol=tol)

  H0=2
  para = c(0,0)
  res1.gamma1 = ASE_sfit(para, H0, z_AS, RHO_AS, A, D, THETA, tauB, tau,
                         maxiter=maxiter, tol=tol)

  H0=0
  para = c(res1.gamma1$KAPPA, res1.gamma1$ETA, res1.eta1$GAMMA)
  res1b = try(ASE_sfit(para, H0, z_AS, RHO_AS, A, D, THETA, tauB, tau,
                       maxiter=maxiter, tol=tol))
  if(class(res1) != 'try-error' & class(res1b) != 'try-error'){
    if(res1b$loglik < res1$loglik)
      res1 = res1b
  }
  p.eta = pchisq(2*( res1.eta1$loglik- res1$loglik ), df =1, lower.tail = F)
  p.gamma = pchisq(2*( res1.gamma1$loglik- res1$loglik ), df =1, lower.tail = F)

  return(list(p.eta = p.eta, p.gamma = p.gamma, loglik.full = - res1$loglik,
              loglik.eta = - res1.eta1$loglik,
              loglik.gamma = - res1.gamma1$loglik,
              KAPPA = res1$KAPPA, ETA = res1$ETA, GAMMA = res1$GAMMA,
              phi = res1$THETA
              #,iters = c(res1$iters, res1.eta1$iters, res1.gamma1$iters)
  ))
}

# ----------------------------------------------------------------------
# Cis-Trans lrt
# ----------------------------------------------------------------------

# CisTrans_lrt <- function(trecase, trec, ase){
#   lrt.stat = -2*( trecase$loglik.full- trec$loglik.full - ase$loglik.full)
#   p.value = pchisq(lrt.stat, df =2, lower.tail = F)
#
#   return(c(lrt.stat = lrt.stat, p.value = p.value))
# }

CisTrans_lrt <- function(trecase, trec_ase_sep){
  lrt.stat = -2*( trecase$loglik.full- trec_ase_sep$loglik)
  p.value = pchisq(lrt.stat, df =2, lower.tail = F)

  return(c(lrt.stat = lrt.stat, p.value = p.value))
}


# ----------------------------------------------------------------------
# TRE & CASE
# ----------------------------------------------------------------------

#KEG_EaseGase: KAPPA, ETA, GAMMA, ETA_ase, GAMMA_ase at log scale

loglik.TReC_ASE_sep <- function(KEG_EaseGase, y,z, z_AS, X, betas, phi, RHO,
                                RHO_AS, tau1, tau2, A,D, THETA, tauB, tau){

  lNB <- loglikNB(para=KEG_EaseGase[1:3], H0=0, y=y, z=z, X=X, betas=betas, phi=phi,
                  RHO=RHO, tau1=tau1, tau2=tau2)
  lASE <- loglikBB(para=KEG_EaseGase[c(1,4:5)],H0=0,A=A,D=D,THETA=THETA, z_AS=z_AS,
                   RHO=RHO_AS, tauB=tauB, tau = tau)

  return(lNB + lASE)

}
grad.TRe_CASE.keg_sep <- function(KEG_EaseGase, y,z, z_AS, X, betas, phi, RHO,
                                  RHO_AS, tau1, tau2, A,D, THETA, tauB, tau){
  grad.nb = Grad_NB(para=KEG_EaseGase[1:3], H0=0, y=y, z=z, X=X, betas=betas, phi=phi,
                    RHO=RHO, tau1=tau1, tau2=tau2)
  grad.bb = Grad_BB_keg(para=KEG_EaseGase[c(1,4:5)],H0=0,A=A,D=D,THETA=THETA, z_AS=z_AS,
                        RHO=RHO_AS, tauB=tauB, tau = tau)
  return(c(grad.nb[1] + grad.bb[1], grad.nb[2:3], grad.bb[2:3]))
}



update_keg_trecase_sep <- function(KEG_EaseGase, H0, y,z, z_AS, X, betas, phi, RHO,
                                   RHO_AS, tau1, tau2, A,D, THETA, tauB, tau){

  lbfgs = optim(par = KEG_EaseGase, fn = loglik.TReC_ASE_sep,
                gr = grad.TRe_CASE.keg_sep,
                y =y, z = z, z_AS=z_AS, X=X, betas = betas, phi = phi,
                RHO = RHO, RHO_AS=RHO_AS,  tau1= tau1,
                tau2 = tau2,A=A,D=D,THETA=THETA,tauB=tauB, tau=tau,
                lower = -50, upper = 10,method = "L-BFGS-B")

  return(lbfgs)
}

#jointly estimate KAPPA
TReCASE_sep_sfit = function(KEG_EaseGase=rep(0,5), y,z,z_AS, X, RHO, RHO_AS,
                            tau1, tau2, A, D, tauB, tau, maxiter =500, tol = 1e-7){
  iter = 0
  #1 initialize all the parameters
  betas =  0
  THETA = phi = 0.05

  KAPPA = exp(KEG_EaseGase[1])
  ETA = exp(KEG_EaseGase[2])
  GAMMA = exp(KEG_EaseGase[3])
  ase_keg = KEG_EaseGase[c(1,4:5)]

  offsets <- compute_offset(z=z, RHO=RHO ,KAPPA= KAPPA, ETA=ETA,
                            GAMMA=GAMMA, tau1=tau1, tau2=tau2)
  if(is.null(dim(X)))
    X = data.matrix(X)

  # update beta phi
  bp = update_beta_phi(y,X, offsets = offsets )
  beta_cur = bp$betas
  phi_cur =  bp$phi

  # update theta
  theta_lbfgs = update_theta(para=ase_keg, H0=0, z_AS=z_AS, RHO_AS=RHO_AS,
                             A=A, D=D, THETA=THETA, tauB=tauB, tau=tau)

  THETA_cur = theta_lbfgs$par

  # update keg
  keg_lbfgs = update_keg_trecase_sep(KEG_EaseGase=KEG_EaseGase, y =y, z=z, z_AS=z_AS, X=X,
                                     betas=beta_cur, phi=phi_cur, RHO=RHO, RHO_AS=RHO_AS,
                                     tau1=tau1, tau2=tau2, A=A,D=D,
                                     THETA=THETA_cur, tauB=tauB, tau=tau)
  KEG_EaseGase_cur = keg_lbfgs$par

  if(any(c(keg_lbfgs$convergence,theta_lbfgs$convergence ) != 0))
    message(sprintf('keg and theta error code: %s, %s',
                    keg_lbfgs$convergence, theta_lbfgs$convergence))

  # update parameters
  while(iter < maxiter & (min(abs(KEG_EaseGase_cur- KEG_EaseGase)) > tol |
                          min(abs(beta_cur - betas)) > tol |
                          abs(phi_cur - phi) >tol |
                          abs(THETA_cur - THETA) > tol)){
    KEG_EaseGase = KEG_EaseGase_cur
    betas = beta_cur
    phi = phi_cur
    THETA = THETA_cur

    KAPPA = exp(KEG_EaseGase_cur[1])
    ETA = exp(KEG_EaseGase_cur[2])
    GAMMA = exp(KEG_EaseGase_cur[3])
    ase_keg = KEG_EaseGase_cur[c(1,4:5)]


    offsets <- compute_offset(z=z, RHO=RHO ,KAPPA= KAPPA, ETA=ETA,
                              GAMMA=GAMMA, tau1=tau1, tau2=tau2)
    bp = update_beta_phi(y,X, offsets = offsets )
    beta_cur = bp$betas
    phi_cur =  bp$phi

    theta_lbfgs = update_theta(para=ase_keg, H0=0, z_AS=z_AS,
                               RHO_AS=RHO_AS, A=A,D=D, THETA=THETA,
                               tauB=tauB, tau=tau)

    THETA_cur = theta_lbfgs$par

    # update keg
    # update keg
    keg_lbfgs = update_keg_trecase_sep(KEG_EaseGase=KEG_EaseGase, y =y, z=z, z_AS=z_AS, X=X,
                                       betas=beta_cur, phi=phi_cur, RHO=RHO, RHO_AS=RHO_AS,
                                       tau1=tau1, tau2=tau2, A=A,D=D,
                                       THETA=THETA_cur, tauB=tauB, tau=tau)
    KEG_EaseGase_cur = keg_lbfgs$par

    # if(any(c(keg_lbfgs$convergence,theta_lbfgs$convergence ) != 0))
    # message(sprintf('keg and theta error code: %s, %s',
    #                 keg_lbfgs$convergence, theta_lbfgs$convergence))

    iter = iter + 1
  }

  return(list(KEG = exp(KEG_EaseGase_cur), phi = phi_cur,
              betas = beta_cur, THETA = THETA_cur,
              loglik =-keg_lbfgs$value, iters = iter,
              conv = c(keg_lbfgs$convergence,theta_lbfgs$convergence)))
}


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
