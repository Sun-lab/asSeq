/*
 *  glm.c
 *  
 *
 *  Created by Wei Sun on 5/11/10.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "utility.h"
#include "glm.h"

/**********************************************************************
 *
 * Some very simple linear algebra functions
 *
 * Originally from mla.cpp of R/CNVTools 
 * 
 
 Package: CNVtools
 Type: Package
 Title: A package to test genetic association with CNV data
 Version: 1.42.0
 Date: 2009-10-06
 Author: Chris Barnes <christopher.barnes@imperial.ac.uk> and 
 Vincent Plagnol <vincent.plagnol@cimr.cam.ac.uk>
 Maintainer: Chris Barnes <christopher.barnes@imperial.ac.uk>
 Description: This package is meant to facilitate the testing of Copy Number 
 Variant data for genetic association, typically in case-control studies.
 License: GPL-3
 Depends: survival
 biocViews: DNACopyNumber,GeneticVariability
 Packaged: 2010-04-23 00:36:41 UTC; biocbuild
 
 **********************************************************************/

/**********************************************************************
 *
 * calculate mean value of y, return residuals (if resid!=0) or center
 *
 * the part about stratum is removed from the original code of CNVtools
 * 
 **********************************************************************/

int wcenter(double *y, int n, double *weight, int resid, double *ynew) 
{
  int i = 0, s=0, empty = 0;
  double swt=0.0, swy=0.0, wi, epsilon=1e-8;
  
  if (weight) {
    for (i=0; i<n; i++) {
      wi   = weight[i];
      swt += wi;
      swy += wi*y[i];
    }
  }else {
    for (i=0; i<n; i++) {
      swy += y[i];
    }
    swt = (double) n;
  }
  swy /= swt;
  
  if (swt>epsilon) {
    if(resid){
      for (i=0; i<n; i++){ ynew[i] = y[i] - swy; }
    }else {
      for (i=0; i<n; i++){ ynew[i] = swy; }
    }
  }else {
    empty = 1;
  }
  
  return(empty);
}

/**********************************************************************
 *
 * Replace y by residual from (weighted) regression through the origin
 * 
 **********************************************************************/

int wresid(double *y, int n, double *weight, double *x, 
           double *ynew, double *beta) 
{
  double swxy, swxx, wi, xi, wx;
  int i;
  
  swxy = swxx = 0.0;
  if (weight) {
    for (i=0; i<n; i++) {
      wi = weight[i];
      xi = x[i];
      wx = wi*xi;
      swxy += wx*y[i];
      swxx += wx*xi;
    }
  } else {
    for (i=0; i<n; i++) {
      xi = x[i];
      swxy += xi*y[i];
      swxx += xi*xi;
    }
  }
  
  if (swxx>0) {
    swxy /= swxx;
    *beta = swxy;
    for (i=0; i<n; i++) {
      if (weight[i] > 0.0) {
        ynew[i] = y[i] - swxy*x[i];
      } else {
        ynew[i] = y[i];
      }
    }
    return(n);
  }else{ 
    return(0);
  }
}

/**********************************************************************
 *
 * Weighted sum of squares
 * 
 **********************************************************************/

double wssq(double *y, int n, double *weights) 
{  
  double res = 0.0, yi;
  int i;
  
  if (weights) {
    for (i=0; i<n; i++) {
      yi = y[i];
      res += weights[i]*yi*yi;
    }
  }else {
    for (i=0; i<n; i++) {
      yi = y[i];
      res += yi*yi;
    }
  }
  
  return(res);
}

/**********************************************************************
 *
 * utilit functions for glm_fit
 *
 * Originally from R/family and glm_test.cpp of R/CNVTools 
 *
 * Codes invovled strata are removed
 * 
 **********************************************************************/

/* 
 
 Variance function
 
 family:
 1    Binomial
 2    Poisson
 3    Gaussian
 4    Gamma
 5    Negagtive Binomial
 */

double varfun(int family, double mu, double phi){
  switch (family) {
    case 1: return((mu*(1.0-mu)));  /* Binomial */
    case 2: return(mu);             /* Poisson */
    case 3: return(1.0);            /* Gaussian */
    case 4: return(mu*mu);          /* Gamma */
    case 5: return(mu + mu*mu*phi); /* Negative Binomial */
    default: return(0.0);
  }
}

/* Valid values for fitted value, mu. 
 
 If, during iteration, an invalid value is returned, the case is omitted 
 
 */

int muvalid(int family, double mu) {
  double minb = 0.0001, maxb = 0.9999, minp = 0.0001;
  double gammaMin = 0.001;
  switch (family) {
    case 1: return(mu>minb && mu<maxb);    /* Binomial */
    case 2: return(mu>minp);               /* Poisson  */
    case 4: return(mu>gammaMin);           /* Gamma    */
    case 5: return(mu>minp);               /* Negative Binomial */
    default: return(1);                    /* Gaussian */
  }
}

/* Link function
 
 Link
 1    Logit
 2    Log
 3    Identity
 4    Inverse
 
 Note that a canonical link shares the code of the corresponding family
 so that the test for canonical link is (link==family)
 
 */

double linkfun(int link, double mu) {
  switch (link) {
    case 1: {     /* Logit */
      if (mu == 1) return  HUGE_VAL;
      if (mu == 0) return -HUGE_VAL;
      return(log(mu/(1.0-mu)));
    }
    case 2: return(log(mu));             /* Log */
    case 3: return(mu);                  /* Identity */
    case 4: return(-1.0/mu);             /* Inverse */
    default: return 0.0;
  }
}

double invlink(int link, double eta) {
  switch (link) {
    case 1: {  /* Logit */
      /*
       if(eta > 30){ return( 1 - 2e-16) }
       if(eta < -30){ return( 2e-16) }
       */
      if (eta ==  HUGE_VAL) return (1); 
      if (eta == -HUGE_VAL) return (0); 
      return(exp(eta)/(1.0+exp(eta))); 
    }
    case 2: return(exp(eta));                /* Log */
    case 3: return(eta);                     /* Identity */
    case 4: return(-1.0/eta);                /* Inverse */
    default: return(0.0);
  }
}

/* dlink = d eta / dmu */
double dlink(int link, double mu) {
  switch (link) {
    case 1: return(-1.0/(mu*(1.0-mu)));    /* Logit */
    case 2: return(1.0/mu);                /* Log */
    case 3: return(1.0);                   /* Identity */
    case 4: return(-1.0/(mu*mu));          /* Inverse */
    default: return 0.0;
  }
}

/**********************************************************************
 *
 * utility functions for glm_fit
 *
 * following the implementation as R/stats/family.R 
 * 
 **********************************************************************/

void initialize(int family, double* y, double* mu, int N, double* nTotal_binom) 
{
  int i;
  
  if (family==BINOMIAL) {         /* Binomial */
    for (i=0; i<N; i++) {
      if (y[i] <0 || nTotal_binom[i] < 0) {
        error("negative values not allowed for the Binomial family");
      }
      if (y[i] > nTotal_binom[i]) {
        error("# of success is larger than # of total trials in Binomial family");
      }
      mu[i] = (y[i] + 0.5)/(nTotal_binom[i]+1.0);
    }
  }else if (family==POISSON) {    /* Poisson */
    for (i=0; i<N; i++) {
      if (y[i] < 0) {
        error("negative values not allowed for the Poisson family");
      }      
      mu[i] = y[i] + 0.1;
    }
  }else if (family==GAUSSIAN) {   /* Gaussian*/
    for (i=0; i<N; i++) {
      mu[i] = y[i];
    }
  }else if (family==GAMMA){       /* Gamma */
    for (i=0; i<N; i++) {
      if (y[i] <= 0) {
        error("non-poistive values not allowed for the Gamma family");
      }      
      mu[i] = y[i] + 0.1;
    }
  }else if (family==NB){          /* Negagtive Binomial */
    for (i=0; i<N; i++) {
      if (y[i] < 0) {
        error("negative values not allowed for the Negative Binomial family");
      }else if (y[i] < 0.01) {
        mu[i] = y[i] + 0.1667;
      }else{
        mu[i] = y[i];
      }
    }
  }else {
    error("invaid family");
  }
}

/**********************************************************************
 *
 * glmFit
 *
 * Originally from glm_test.cpp of R/CNVTools 
 *
 * There are several channges:
 * 
 * (1)  For binomial distribution, allow two sets of inputs, 
 *      y and nTotal_binom, i.e. the number of succes and the total 
 *      number of trials
 * 
 * (2)  Codes invovled strata are removed
 * 
 * (3)  The way to handle the situation of no covaraite
 * 
 * (4)  The way to handel least square problem with Guassian family
 * 
 * (5)  Initial values of mu
 *
 * (6)  Add support for negative binomial distribution
 *
 
 Package: CNVtools
 Type: Package
 Title: A package to test genetic association with CNV data
 Version: 1.42.0
 Date: 2009-10-06
 Author: Chris Barnes <christopher.barnes@imperial.ac.uk> and 
 Vincent Plagnol <vincent.plagnol@cimr.cam.ac.uk>
 Maintainer: Chris Barnes <christopher.barnes@imperial.ac.uk>
 Description: This package is meant to facilitate the testing of Copy Number 
 Variant data for genetic association, typically in case-control studies.
 License: GPL-3
 Depends: survival
 biocViews: DNACopyNumber,GeneticVariability
 Packaged: 2010-04-23 00:36:41 UTC; biocbuild
 
 
 Input:
 
 family       GLM family (see below)
 link         Link function (see below)
 N            # units
 M            # X variables
 y            y-variable (N-vector)
 X            If M>0, N*M matrix of X variables
 maxit        Maximum number of iterations of IRLS algorithm
 conv         Proportional change in weighted sum of squares residuals to
              declare convergence
 init         If true (non-zero), the iteration starts from initial estimates 
              of fitted values (see below). This option has no effect if
              no iteration is required
 
 Output:
 
 rank         rank of X after regression on strata
 Xb           orthogonal basis for X space (N*rank matrix)
 fitted       fitted values 
 resid        working residuals (on linear predictor scale) (N-vector)
 weights      weights (N-vector)
 scale        scale factor (scalar)
 df_resid     residual degrees of freedom
 beta         regression coeficien of the last covariate

 Return
 
 0            convergence
 1            no convergence after maxit iterations
 
 **********************************************************************/

int glmFit(int* familyR, int* linkR, int* dims, int* nIter,
           double *y, double *offset, double *z,
           double *X, double *nTotal_binom, double *convR, 
           int *rank, double *Xb, double *fitted, double *resid, 
           double *weights, double *phi, int* trace, 
           double *scale, int *df_resid, double* beta) 
{
  double epsilon = 1e-8;       /* Singularity threshold */
  int N, M, maxit, init, useOffset;
  int i = 0, j=0, Nu, dfr = 0, irls, empty = 0;
  int x_rank = 0, convg = 0, iter = 0;
  
  N = dims[0];
  M = dims[1];
  maxit = dims[2];
  init  = dims[3];
  useOffset = dims[4];
  
  int family  = *familyR;
  int link    = *linkR;
  double conv = *convR;
  
  double mu, ri, wi, D, Vmu, wss, wss_last=0.0, ssx=0.0, ssr=0.0;
  double *xi, *xbi, *xbj, zmu, wsum;
    
  if(family > 6 || family < 1){
    Rprintf("family=%d, ", family);
    error("Invalid family!\n");
  }
  
  /* Is iteration necessary? */
  irls = !((family==GAUSSIAN) && (link==IDENTITY));
  
  /* ----------------------------------------------------------------*
   * by default, initialize mu (fitted) by y itself, with neccesary 
   * modification, e.g., y + 0.1 to avoid log(0) for Poisson family
   * ----------------------------------------------------------------*/
  if (!init) {
    initialize(family, y, fitted, N, nTotal_binom);
  }
  
  /* ----------------------------------------------------------------*
   * Initialize wi (weights) and (standardized) residual 
   * In IRLS, z_i = eta_i + (y_i - mu_i)*(d_eta/d_mu)_i
   *
   * In the following code:
   * ri = (y_i - mu_i)*(d_eta/d_mu)_i
   * wi = Vmu is the weight,  
   *
   * for those invlaid values, set their weight to be 0.
   * ----------------------------------------------------------------*/
  Nu      = 0;
  wsum    = 0.0;
  
  for (i=0; i<N; i++) {

    mu = fitted[i];
    
    if (!muvalid(family, mu)) {
      wi = ri = 0.0;
    }else {
      Nu ++;
      Vmu = varfun(family, mu, *phi);
      
      if (link == family) {
        ri = (y[i] - mu)/Vmu;
        wi = Vmu;
      }else {
        D  = dlink(link, mu);
        ri = D*(y[i] - mu);
        wi = 1.0/(D*D*Vmu);
      }
    }
    
    weights[i] = wi;
    resid[i]   = ri;
    if (weights[i] < epsilon) weights[i] = 0.;
    
    wsum += weights[i];
  }
  
  /* ----------------------------------------------------------------*
   * If summation of all weights is too small, stop 
   * ----------------------------------------------------------------*/
  
  if (wsum < epsilon) {
    if(*trace)
      Rprintf("  glmFit: empty model, summation of all weights are small!\n");
    
    /* set M = -1, so that no extra computation is done */
    M = -1;
  }
  
  if(*trace > 3){
    Rprintf("\n  glmFit: finish initialization, N=%d, M=%d, family=%d, and irls=%d\n", 
            N, M, family, irls);
  }
  
  /* ----------------------------------------------------------------*
   * If M=0, there is only an intercept 
   * ----------------------------------------------------------------*/
  
  if (M == 0){
    zmu = 0.0;
    
    /* all the observations have the same fitted value  */
    for (i=0; i<N; i++) {
      /**
       * z_i = current estimate of eta + (y-mu)/gradient
       * where
       * linkfun(link, fitted[i]) = eta
       * resid[i] = (y-mu)/gradien
       */
      z[i] = linkfun(link, fitted[i]) + resid[i];
      if (useOffset) { z[i] -= offset[i]; }
      
      zmu += z[i]*weights[i];
    }
    
    zmu = zmu/wsum;
    dfr = Nu;
    
    if (useOffset) {
      for (i=0; i<N; i++) {
        fitted[i] = invlink(link, zmu + offset[i]);
        resid[i]  = y[i] - fitted[i];
      }        
    }else {
      mu = invlink(link, zmu);
      for (i=0; i<N; i++) {
        fitted[i] = mu;
        resid[i]  = y[i] - mu;
      }        
    }
    
    if (family>2){
      *scale = wssq(resid, N, weights)/dfr;
    }else{
      *scale = 1.0;
    }
    
    x_rank = 0;
    
  }else if (M > 0) {
    /* ----------------------------------------------------------------*
     * If M>0, include covariates 
     * ----------------------------------------------------------------*/
    
    convg    = 0;
    iter     = 0;
    wss_last = 0.0;
    
    if (!irls) {
      /* Simple linear Gaussian case */
      xi  = X;
      xbi = Xb;
      x_rank = 0;
      
      for (i=0; i<M; i++, xi+=N) {
        wcenter(xi, N, weights, 1, xbi);
        ssx = wssq(xbi, N, weights);
        
        xbj = Xb;
        
        for (j=0; j<x_rank; j++, xbj+=N)  
          wresid(xbi, N, weights, xbj, xbi, beta);
        
        ssr = wssq(xbi, N, weights);
        
        if (ssr/ssx > epsilon) {
          wresid(resid, N, weights, xbi, resid, beta);
          x_rank++;
          xbi+=N;
        }
      }
      
      /* obtain the fitted values */
      for (i=0; i<N; i++) fitted[i] = y[i] - resid[i];
      
      wss_last = wssq(resid, N, weights);
      
    }else{
      
      /* IRLS algorithm */
      while(iter<maxit && !convg) {
        
        if (*trace > 9) {
          Rprintf("    glmFit: iteration %d: \n", iter);
        }
        
        for (i=0; i<N; i++) {
          /**
           * current estimate of eta + (y-mu)/gradient
           *
           * linkfun(link, fitted[i]) = eta
           * resid[i] = (y-mu)/gradien
           */
          z[i] = linkfun(link, fitted[i]) + resid[i];
        }
        
        if (useOffset) {
          for (i=0; i<N; i++) z[i] -= offset[i];
        }
                
        empty = wcenter(z, N, weights, 1, resid);  //removes the mean from z
        
        if (empty == 1) {
          if(*trace > 0)
            Rprintf("  glmFit: empty model, summation of all weights are small!\n");
          
          break;
        }
        
        if (*trace > 9) {
          Rprintf("    glmFit: iteration %d:, initialized z\n", iter);
        }
        
        
        /**
         * tries to fit the regression line (no intercept) to the residuals
         */
        
        xi  = X;
        xbi = Xb;
        x_rank = 0;
        
        for (i=0; i<M; i++, xi+=N) {
          wcenter(xi, N, weights, 1, xbi);
          ssx = wssq(xbi, N, weights);
          xbj = Xb;
          
          for (j=0; j<x_rank; j++, xbj+=N) wresid(xbi, N, weights, xbj, xbi, beta);
          
          ssr = wssq(xbi, N, weights);
                    
          if (ssr/ssx > epsilon) {
            /**
             * takes the residuals after fitting the regression line (no intercept) 
             * to the mean value 
             */
            wresid(resid, N, weights, xbi, resid, beta); 
            x_rank++;
            xbi+=N;
          }
        }
        
        if (*trace > 9) {
          Rprintf("    glmFit: iteration %d:, got Xb\n", iter);
        }
        
        /* well, it is question whether we should give error or just warning here */
        if(x_rank < M){ 
          if(*trace > 1){
            Rprintf("  glmFit: x_rank=%d, M=%d, X is not full rank\n", x_rank, M);
          }
          
          break;
        }
        
        wss = 0.0;
        Nu  = 0;
        
        for (i=0; i<N; i++) {
          
          if (useOffset) {
            mu = invlink(link, z[i] + offset[i] - resid[i]);
          }else {
            mu = invlink(link, z[i] - resid[i]);
          }
          
          fitted[i] = mu;
          
          if (weights[i] <= 0.0) {
            wi = ri = 0.0;
          } else {
            
            Vmu = varfun(family, mu, *phi);
            Nu ++;
            
            if (link == family) {
              ri = (y[i] - mu)/Vmu;
              wi = Vmu;
            }else {
              D = dlink(link, mu);
              ri = D*(y[i] - mu);
              wi = 1.0/(D*D*Vmu);
            }
            wss += wi*ri*ri;
            
            weights[i] = wi;
            resid[i]   = ri;
            if (weights[i] < epsilon) weights[i] = 0.;
          }
        }
        
        if(wss > 1.0/epsilon){
          if(*trace > 1)
            Rprintf("  glmFit: huge wss, indicting failt to fit the model!\n");
          
          break;
        }
        
        if (*trace > 5) {
          Rprintf("    glmFit: iteration %d, wss=%.3f\n", iter, wss);
        }
        
        convg = (Nu<=0) || (iter && (fabs(wss-wss_last)/(wss_last + 0.1) < conv));
        wss_last = wss;
        iter ++;
      }
    }
    
    if (convg) {
      /* assume there is an intercept */
      dfr = Nu - 1 - x_rank;
      
      if (family > 2) {
        *scale = wss_last/(dfr);
      }else{
        *scale = 1.0;
      }
    }else {
      dfr    = 0.0;
      *scale = 1.0;
    }
  }
  
  *df_resid = dfr>0? dfr : 0;
  *rank     = x_rank;
  *nIter    = iter;
    
  return(irls && convg);
}

/**********************************************************************
 
 GLM score test 
 
 Input:
 
 P         Number of new explanatory variables to be added 
 Z         N*P matrix containing covariates
 
 For all other input arguments, see glm_fit, but note that M now coincides 
 with rank -- the number of columns in Xb
 
 Output:
 
 chi2  Score test 
 df    Degrees of freedom for asymptotic chi-squared distribution
 
 **********************************************************************/
            
void glm_score_test(int* dims, double *Z, double *resid, 
                    double *weights, double *Xb, double* scaleR,
                    double *chi2, int *df) 
{
  
  int i = 0, j = 0, rank, N, M, P;
  
  double epsilon1 = 1.e-8;   /* First stage singularity test */
  double *Zi = Z;
  double *Zr, *Zri, ssz, ssr, *Zrj, Zrij, ws, wss, wz, test = 0.0;
  double *Xbj, beta;
  double scale = *scaleR;
  
  N = dims[0];
  M = dims[1];
  P = dims[2];
  
  /* Work array */
  Zr  = (double *) Calloc(N*P, double);
  Zri = Zr;
  
  /* Main algorithm */
  rank = 0;
  
  for (i=0; i<P; i++, Zi+=N) {
    /* Regress each column of Z on X basis */
    wcenter(Zi, N, weights, 1, Zri);
    ssz = wssq(Zri, N, weights);
    Xbj = Xb;
    
    for (j=0; j<M; j++, Xbj+=N){
      wresid(Zri, N, weights, Xbj, Zri, &beta);
    }
    
    ssr = wssq(Zri, N, weights);
    
    if (ssr/ssz > epsilon1) {     /* First singularity test */
      Zrj = Zr;
      for (j=0; j<rank; j++, Zrj+=N){
        wresid(Zri, N, weights, Zrj, Zri, &beta);
      }
      
      /* Sum and sum of squares */
      ws = 0.0, wss = 0.0;
      for (j=0; j<N; j++) {
        Zrij = Zri[j];
        wz   = weights[j]*Zrij;
        ws  += wz*resid[j];
        wss += Zrij*wz;
      }
      
      /* Second singularity test */
      if (wss/ssr > epsilon1) {
        test += ws*ws/(scale*wss);
        rank++;
        Zri  += N;
      }else {
        error("colinearity in added covaraites Z\n");
      }
      
    }
  }
  
  *chi2 = test;
  *df   = rank;
  
  Free(Zr);
}

/**********************************************************************
 *
 * glmNB
 *
 * Originally from MASS/negbin.R 
 * 
 
 Package: MASS
 Priority: recommended
 Version: 7.3-5
 Date: 2010-01-03
 Depends: R (>= 2.10.1), grDevices, graphics, stats, utils
 Suggests: lattice, nlme, survival
 Author: S original by Venables & Ripley. R port by Brian Ripley
 <ripley@stats.ox.ac.uk>, following earlier work by Kurt Hornik
 and Albrecht Gebhardt.
 Maintainer: Brian Ripley <ripley@stats.ox.ac.uk>
 Description: Functions and datasets to support Venables and Ripley,
 'Modern Applied Statistics with S' (4th edition).
 Title: Main Package of Venables and Ripley's MASS
 License: GPL-2 | GPL-3
 URL: http://www.stats.ox.ac.uk/pub/MASS4/
 LazyLoad: yes
 LazyData: yes
 Packaged: 2010-01-03 10:50:27 UTC; ripley
 Repository: CRAN
 Date/Publication: 2010-01-03 14:05:40
 
 **********************************************************************/

/**********************************************************************
 *
 * log likelihood of Poisson
 *
 **********************************************************************/

double loglik_Poisson(int N, double* mu, double* y){
  int i;
  double yi, mui, logL = 0.0;
  
  for (i=0; i<N; i++) {
    yi  = y[i];
    mui = mu[i];
    
    logL += (yi*log(mui) - mui - lgammafn(yi + 1.0));
  }
  
  return(logL);
}

/**********************************************************************
 *
 * log likelihood of negative binomial
 *
 **********************************************************************/

double loglik_NB(int N, double phi, double* mu, double* y){
  int i;
  double logL1, logL, yi, mui;
  double th = 1.0/phi;
  double logL0 = th*log(th) - lgammafn(th);
  
  logL = 0.0;
  
  for (i=0; i<N; i++) {
    yi  = y[i];
    mui = mu[i];
    
    if (yi==0) {
      logL1  = th*log(th) - th*log(th + mui);
    }else {
      logL1  = lgammafn(th + yi) - lgammafn(yi + 1.0) + yi*log(mui) - (th + yi)*log(th + mui);
      logL1 += logL0;
    }

    logL += logL1;
  }
  
  return(logL);
}

/**********************************************************************
 *
 * score_info
 *
 * score and Fisher information, i.e., the first and negative second 
 * derivative of likelihood of phi 
 *
 **********************************************************************/

void score_info(int N, double theta, double* mu, double* y, 
                double* score, double* info)
{
  int i;
  double score1=0.0, info1=0.0;
  double mui, yi, scorei, infoi, thMui;
  
  for (i=0; i<N; i++) {
    yi  = y[i];
    mui = mu[i];
    
    thMui   = theta + mui;
    scorei  = digamma(yi + theta) - digamma(theta) - (theta + yi)/thMui;
    score1 += (scorei - log(thMui) + 1 + log(theta));
    
    infoi   = trigamma(theta) - trigamma(yi + theta) + (mui - yi)/(thMui*thMui);
    info1  += (infoi + 1/thMui - 1/theta);
  }
  
  *score = score1;
  *info  = info1;
}

/**********************************************************************
 *
 * phi_ml
 *
 * MLE of phi (over-dispersion parameter), given mu 
 *
 * Actually we find MLE of 1/phi here and then take inverse
 *
 **********************************************************************/

int phi_ml(double* y, double* mu, int N, int limit, double eps, 
           double* phi, int initPhi, int trace)
{
  double theta0, del, tmp;
  double score=0.0;
  double info=0.0;
  int i, it=0;
  double minTheta = 1e-5;
  double maxTheta = 1.0/minTheta;
  int tryPoisson = 0;
  int tryZINB = 0;
  int fail = 0;
  
  if(initPhi){
    theta0 = 1.0/(*phi);
  }else{
    theta0 = 0.0;
    for (i=0; i<N; i++) {
      tmp = y[i]/mu[i] - 1.0;
      theta0 += tmp*tmp;
    }
    theta0 = (double)N/theta0;
  }
  
  it  = 0;
  del = 1.0;
  
  if(trace > 5) Rprintf("  phi.ml: initial phi = %.2e\n", 1/theta0);
  
  while(it < limit && fabs(del) > eps) {
    score_info(N, theta0, mu, y, &score, &info);
    del     = score/info;
    theta0 += del;
    it     += 1;
    
    if(trace > 5) Rprintf("  phi.ml: iter %d, phi=%.2e, score=%.2e, info=%.2e\n", 
                          it,  1/theta0, score, info);
    
    if (theta0 > maxTheta) {
      theta0 = maxTheta;
      if(trace > 3)
        Rprintf("    phi is truncated at %.2e, no overDispersion?\n", 1/maxTheta);
      
      tryPoisson = 1;
      break;
    }
    
    if(theta0 < minTheta) {
      theta0 = minTheta;
      if(trace > 3)
        Rprintf("    phi is truncated at %.2e, too much overDispersion?\n", 1/minTheta);
      
      tryZINB = 1;
      break;
    }
    
  }
  
  if(it == limit) {
    fail = 1;
    if(trace > 3)
      Rprintf("  phi.ml: iteration limit reached in phi_ml\n");
  }
  
  *phi = 1/theta0;
  
  return(tryPoisson + 2*tryZINB + 4*fail);
}

/**********************************************************************
 *
 * main function of glmNB
 *
 **********************************************************************/

int glmNB(int *dims, int *nIter, double *y, double *z, 
          int *linkR, double *offset, double *X, double *convR, 
          int *rank, double *Xb, double *fitted, double *resid, 
          double *weights, double *phi, double *scale, 
          int *df_resid, int* family, double *twologlik, 
          double *scoreTestP, int *trace, double *beta)
{
  int N, M, maxit, init, useOffset;
  int i, succeed = 0, iter = 0, convged=0;
  double conv = *convR, initPhi;
  int fam0, cv = 0; /* cv is convergence indicator */
  double nTotal_binom=0.0;  /* only used for binomial link, NOT useful here */
  double del, Lm, Lm0, phi0;
  double score, scoreNum, scoreDen, scorePval, yi, mui;
  
  /* convergence indicator for phi 
   * if cvPhi = 0, NB model is OK. 
   * if cvPhi = 1, suggest we need to use Poisson
   * if cvPhi = 2, suggest we need to use ZINB
   */
  int cvPhi;
  
  N = dims[0];
  M = dims[1];
  maxit = dims[2];
  init  = dims[3];
  useOffset = dims[4];
  
  if(*trace > 3) 
    Rprintf("\n  glmNB: N=%d, M=%d, maxit=%d, init=%d, useOffset=%d\n", 
            N, M, maxit, init, useOffset);

  /**
   * if there is no initial values, ignore the parameter *family 
   * and start with a Poission model
   */
  if(!init){
    /* Initial fit */
    
    fam0 = POISSON;
    
    cv = glmFit(&fam0, linkR, dims, nIter, y, offset, z, 
                X, &nTotal_binom, convR, rank, Xb, fitted, resid, 
                weights, phi, trace, scale, df_resid, beta);
    
    if (cv==0) {
      if(*trace){
        Rprintf("\n  glmNB: fail to converge in initial glmFit by Poission regression\n");
      }
      
      maxit   = -1;
      succeed = 0;
      return(succeed);
    }else {
      
      /* test for overdispersion by Dean's Score test */
      scoreNum = 0.0;
      scoreDen = 0.0;
      
      for (i=0; i<N; i++) {
        yi  = y[i];
        mui = fitted[i];
        scoreNum += (yi - mui)*(yi - mui) - yi;
        scoreDen += mui*mui;
      }
      
      score = scoreNum/sqrt(2.0*scoreDen);
      
      /**
       * double pnorm(double x, double mu, double sigma, int lower_tail, int give_log);
       */
      scorePval = pnorm(score, 0.0, 1.0, 0, 0);
      
      if(*trace > 3) 
        Rprintf("\n  overdispersion score = %.2e, p-value = %.2e\n\n", score, scorePval);
      
      if(scorePval > *scoreTestP){
        *family    = POISSON;
        Lm         = loglik_Poisson(N, fitted, y);
        *twologlik = 2.0*Lm;
        *phi       = 0.0;
        
        maxit   = -1;
        succeed = 1;
        return(succeed);

      }else {
        fam0    = NB;
        *family = NB;
      
        /**
         * calculate phi by MLE, without initial values of phi
         */
        cvPhi = phi_ml(y, fitted, N, maxit, conv, phi, 0, *trace);
        
        if(cvPhi==0){
          if(*trace > 3) 
            Rprintf("\n  initial value for phi: %e\n", *phi);
        }else if (cvPhi==1){
          if(*trace > 3) 
            Rprintf("\n  Choose Poisson model due to small phi: %e\n", *phi);
          
          *family    = POISSON;
          Lm         = loglik_Poisson(N, fitted, y);
          *twologlik = 2.0*Lm;
          *phi       = 0.0;
          
          maxit   = -1;
          succeed = 1;
          return(succeed);

        }else if(cvPhi==2){
          if(*trace > 1) 
            Rprintf("\n  The overdispersion parameter is too large: phi=%e\n", *phi);
          
          *family = ZINB;
          
          maxit   = -1;
          succeed = 0;
          return(succeed);

        }else { /* estimation of phi fail to converge */
          if(*trace > 1) 
            Rprintf("\n  Estimation of phi fail to converge: phi=%e\n", *phi);
          
          maxit   = -1;
          succeed = 0;
          return(succeed);
        }
      }
    }
  }else {
    /**
     * if there is initial values, 
     * there must be both initial fitted values and inital phi
     */
    fam0 = *family;
    
    cv = glmFit(&fam0, linkR, dims, nIter, y, offset, z, 
                X, &nTotal_binom, convR, rank, Xb, fitted, resid, 
                weights, phi, trace, scale, df_resid, beta);
    
    if (cv==0) {
      if(*trace > 1){
        Rprintf("\n  glmNB: fail to converge using initial values, fam0=%d\n", fam0);
      }
      
      maxit   = -1;
      succeed = 0;
      return(succeed);
      
    }else{
      /**
       * no need to go further if this is a Poission regression
       */
      
      if(fam0 == POISSON){
        Lm         = loglik_Poisson(N, fitted, y);
        //Rprintf("\n===--->>> fam0=%d, cv=%d, Lm=%e\n", fam0, cv, Lm);

        *twologlik = 2.0*Lm;
        *phi       = 0.0;

        maxit   = -1;
        succeed = 1;
        return(succeed);

      }else {
        cvPhi = phi_ml(y, fitted, N, maxit, conv, phi, 0, *trace);
        
        if(cvPhi==0){
          if(*trace > 3) 
            Rprintf("\n  glmNB: initial value for phi: %e\n", *phi);
        }else {
          if(*trace > 1) 
            Rprintf("\n  glmNB: fail to estimate phi\n");
          
          maxit   = -1;
          succeed = 0;
          return(succeed);
        }
      }
    
    }
  }
  
  if(maxit > 0){
    Lm   = loglik_NB(N, *phi, fitted, y);
  }else{
    Lm   = 0.0;
  }

  del  = 1.0;
  Lm0  = Lm + 1.0;
  iter = 0;
  convged = 0;
  succeed = 0;
  
  while (iter < maxit && (!convged)) {
    
    dims[3] = 1; /* use initial values */
    
    cv = glmFit(&fam0, linkR, dims, nIter, y, offset, z, 
                X, &nTotal_binom, convR, rank, Xb, fitted, resid, 
                weights, phi, trace, scale, df_resid, beta);
    
    if (cv==0) { 
      if(*trace > 1) 
        Rprintf("\n  glmNB: fail to converge in glmFit\n");
      
      break;
    }
    
    phi0  = *phi;
    cvPhi = phi_ml(y, fitted, N, maxit, conv, phi, 1, *trace);
    
    if(cvPhi==0){
      if(*trace > 3) 
        Rprintf("\n  finish phi_ml, cvPhi=%d, phi=%e\n", cvPhi, *phi);
    }else {
      if(*trace > 1) 
        Rprintf("  glmNB: fail in phi_ml\n");
      
      break;
    }
    
    del = phi0 - *phi;
    Lm0 = Lm;
    Lm  = loglik_NB(N, *phi, fitted, y);
        
    if (*trace > 3) {
      Rprintf("\n  Phi(%d) = (%e, %e), logLik = (%e, %e)\n\n", 
              iter, phi0, *phi, Lm0, Lm);
    }
    
    convged = fabs(Lm0 - Lm) + fabs(del) < conv;
    
    iter++;
  }

  if(iter == maxit) {
    if (*trace) {
      Rprintf("\n  glmNB: Alternation limit reached: iter=%d\n", iter);
    }
  }else if (convged) {
    succeed = 1;
  }

  *twologlik = 2.0*Lm;
  
  return(succeed);
}

/**********************************************************************
 *
 * logL_b_TReC: log likelihood of TReC model regarding to b
 *
 **********************************************************************/

double logL_b_TReC(double b1, int N, int fam, double b0, double phi,  
                   double* y, double* x, double* mu, double* mu1)
{
  int i;
  double logL, xi, cst0, cst1;
  
  cst0 = exp(b1 - b0);
  cst1 = (1.0 + exp(b1))/(1.0 + exp(b0));
    
  for (i=0; i<N; i++) {
    xi = x[i];
    
    if (fabs(xi - 0.0) < 0.01) {
      mu1[i] = mu[i];
    }else if (fabs(xi - 2.0) < 0.01) {
      mu1[i] = mu[i]*cst0;
    }else if (fabs(xi - 1.0) < 0.01) {
      mu1[i] = mu[i]*cst1;
    }else {
      error("invalid genotype\n");
    }
  }
  
  logL = 0.0;
  if (fam == POISSON) {
    logL = loglik_Poisson(N, mu1, y);
  }else if (fam == NB) {
    logL = loglik_NB(N, phi, mu1, y);
  }else {
    error("wrong family\n");
  }
  
  return(logL);
}

/**********************************************************************
 *
 * grad_b_TReC: 1st and 2nd gradient of logL_b_TReC
 *
 **********************************************************************/

void grad_b_TReC(double b1, int N, int fam, double b0, double phi, 
                 double* y, double* x, double* mu, double* gr)
{
  int i;
  double ti, mui, xi, grad1, grad2, df_dmu1, df_dmu2, dmu_db;
  double cst0, cst1, cst2;
  
  cst0 = exp(b1 - b0);
  cst1 = (1.0 + exp(b1))/(1.0 + exp(b0));
  cst2 = exp(b1)/(1.0 + exp(b1));
  
  grad1 = grad2 = 0.0;
  
  for (i=0; i<N; i++) {
    xi = x[i];
    ti = y[i];
    
    if (fabs(xi - 0.0) < 0.01) {
      continue;
    }else if (fabs(xi - 2.0) < 0.01) {
      mui    = mu[i]*cst0;
      dmu_db = mui;
    }else if (fabs(xi - 1.0) < 0.01) {
      mui    = mu[i]*cst1;
      dmu_db = mui*cst2;
    }else {
      error("invalid genotype\n");
    }
    
    df_dmu1 = ti/mui - (1.0 + phi*ti)/(1.0 + phi*mui);        
    grad1  += df_dmu1*dmu_db;
    
    df_dmu2 = -ti/(mui*mui) + phi*(1.0 + phi*ti)/(1.0 + phi*mui)/(1.0 + phi*mui);        
    grad2  += (df_dmu2*dmu_db*dmu_db + df_dmu1*dmu_db);
  }
  
  gr[0] = grad1;
  gr[1] = grad2;
}

/**********************************************************************
 * b_TReC_ml: find the MLE of b in TReC model by 
 * a Newton-Raphson algrorithm
 **********************************************************************/

void b_TReC_ml(double* b_xj, int N, int fam, double b0, double phi, 
               double* y,double* x, double* mu, int limit, 
               double eps, int trace, int* failR)
{
  int it, fail;
  double logL, b1, del, gr[2];
  double min_b = -1e5;
  double max_b =  1e5;
  
  it   = 0;
  del  = 1.0;
  fail = 0;
  b1   = b0;
  gr[0] = 1.0;
  gr[1] = 1.0;
    
  if(trace > 5){
    Rprintf("b_TReC_ml: limit=%d, eps=%e\n", limit, eps);
  }
  
  while (it < limit && fabs(del) > eps) {
    
    grad_b_TReC(b1, N, fam, b0, phi, y, x, mu, gr);
    
    del = gr[0]/gr[1];
    b1 -= del;
    it += 1;
    
    if (it > 10 && fabs(del) > 0.1) {
      fail = 1;
      break;
    }
    
    if(trace > 5){
      Rprintf("      b_TReC_ml: it=%d, b1=%e, del=%e, gr[0]=%e, gr[1]=%e\n", 
              it, b1, del, gr[0], gr[1]);
    }
    
    if (b1 > max_b) {
      b1 = max_b;
      if(trace > 3)
        Rprintf("    Estimate of b_xj is truncated at %.2e\n", max_b);
      
      break;
    }
    
    if(b1 < min_b) {
      b1 = min_b;
      if(trace > 3)
        Rprintf("    Estimate of b_xj is truncated at %.2e\n", min_b);
      
      break;
    }  
  }
  
  if (it == limit) {
    // fail = 1;
    
    if(trace > 1)
      Rprintf("  b_TReC_ml: iteration limit reached in b_TReC_ml\n");
  }
  
  if (fabs(gr[0]) > 0.01) {
    fail = 1;
    
    if(trace > 1)
      Rprintf("  b_TReC_ml: 1st derivative = %.2e\n", gr[0]);
  }
  
  if (gr[1] > 1e-16) {
    fail = 1;
    
    if(trace > 1)
      Rprintf("  b_TReC_ml: 2st derivative = %.2e\n", gr[1]);
  }
    
  *b_xj  = b1;
  *failR = fail;
}

/**********************************************************************
 *
 * main function of glmNBlog
 * 
 * the last covariate does not follow a log linear model 
 * but a linear model. So iterately estimate all other co
 *
 **********************************************************************/

int glmNBlog(int *dimsNew, int *nIter, double *pY, double *z, 
             int *linkR, double *offset, double *pX, double *conv, 
             double *convGLM, int *rank, double *Xb, double *fitted, 
             double *resid, double *weights, double *phi, double *scale, 
             int *df_resid, int* family, double *twoLL_trec, 
             double *scoreTestP, int *trace, double *beta,
             double *fitted2, double *offsetN)
{
  int i, g, k, k1, k2, k3, k4, df1, convBase, convSNPj, succeed, fail;
  int N     = dimsNew[0];
  int nX    = dimsNew[1];
  int maxit = dimsNew[2];
  int init  = dimsNew[3]; /* whetehr to use initial values */
  int useOffset = dimsNew[4];
  double gr[2];
  
  double *pZ, pZk, twoLL_trec0, twoLL_trec1, nless5;
  double bxj_old, bxj_new, btmp, phi_old, paraDiff;
  
  succeed = 1;
  for (k=0; k<N; k++) { fitted2[k] = fitted[k]; }
  
  pZ = pX + N*nX;
  
  phi_old = *phi;
  bxj_old = bxj_new = 0.0;

  /*************************************************
   * must use initial values
   *************************************************/
  
  if(!init){
    error("we expect to use initial values\n");  
  }
  
  /*************************************************
   * phi = 0 for Poisson distribution
   *************************************************/
  
  if(*family==POISSON && phi_old > 0){
    error("phi=%e while family is Poisson\n", phi_old);  
  }
  
  /*************************************************
   * calculate old likelihood       
   *************************************************/
  if (*family==NB) {
    twoLL_trec0 = 2.0*loglik_NB(N, phi_old, fitted2, pY);
  }else if (*family==POISSON) {
    twoLL_trec0 = 2.0*loglik_Poisson(N, fitted2, pY);
  }else {
    error("invalid family\n");
  }
  
  for(g=0; g<maxit; g++){
    
    paraDiff = 0.0;

    /* --------------------------------------------------------
     * 1. Given theta, phi, b_0, b_k, and b_u, estimate b_{x_j}
     * -------------------------------------------------------*/
        
    fail = 0;

    b_TReC_ml(&bxj_new, N, *family, bxj_old, phi_old, pY, pZ, 
              fitted2, maxit, *convGLM, *trace, &fail);
    
    if(fail){
      if(*trace){
        Rprintf("\n  g=%d, fail b_TReC_ml in glmNBlog, fail=%d\n", g, fail);
      }
      
      succeed=0;
      break;
    }
    
    if (paraDiff < fabs(bxj_new - bxj_old)) paraDiff = fabs(bxj_new - bxj_old);
    
    /**********************************
     * calculate new likelihood       
     **********************************/
    if (*family==NB) {
      twoLL_trec1 = 2.0*loglik_NB(N, phi_old, fitted2, pY);
    }else if (*family==POISSON) {
      twoLL_trec1 = 2.0*loglik_Poisson(N, fitted2, pY);
    }else {
      error("invalid family\n");
    }
    
    // Rprintf("\n===--->>> family=%d, twoLL_trec1=%e\n", *family, twoLL_trec1);

    if (*trace > 1) {
      Rprintf("\n  g=%d, phi=%e, b_old=%e, b_new=%e\n", g, phi_old, bxj_old, bxj_new);
      Rprintf("  twoLL_trec(old, new)=(%e, %e)\n", twoLL_trec0, twoLL_trec1);
    }
    
    if((twoLL_trec1 - twoLL_trec0)/(fabs(twoLL_trec0) + 1.0) < -1e-4){
      
      nless5 = 0.0;
      for (k=0; k<N; k++) { 
        if (pY[k] < 5) nless5 += 1.0;
      }
      
      Rprintf("  g=%d, nless5=%.1f, twoLL(old, new)=(%e, %e)\n", 
              g, nless5, twoLL_trec0, twoLL_trec1);

      if (nless5 > 0.75*N) {
        succeed=0;
        break;
      }
      
      error("\n  likelihood decreases during update of b_xj in glmNBlog\n");
    }
    
    bxj_old = bxj_new;
    
    /* --------------------------------------------------------
     * 2. Given b_{x_j} and phi, estimate b_0, b_k, b_u, phi
     * in fact, we do not estimate b_0, b_k, and b_u, 
     * instead we only need to estimate mu, or fitted2.
     * -------------------------------------------------------*/
    
    twoLL_trec0 = twoLL_trec1;
    
    dimsNew[1] = nX;
    dimsNew[3] = 1; /* use initial values */
    dimsNew[4] = 1; /* use offset         */
    
    for (k=0; k<N; k++){
      pZk = pZ[k];
      
      if (fabs(pZk) < 0.01) {
        offsetN[k] = offset[k];
      }else if (fabs(pZk - 1.0) < 0.01) {
        offsetN[k] = offset[k] + log(0.5*(1 + exp(bxj_old)));
      }else if (fabs(pZk - 2.0) < 0.01) {
        offsetN[k] = offset[k] + bxj_old;
      }else {
        error("invalid value of genotype.");
      }
    }
    
    convSNPj = glmNB(dimsNew, nIter, pY, z, linkR, offsetN, pX, convGLM, 
                     rank, Xb, fitted2, resid, weights, phi, scale, &df1, 
                     family, &twoLL_trec1, scoreTestP, trace, &btmp);

    // Rprintf("\n===--->>> family=%d, twoLL_trec1=%e\n", *family, twoLL_trec1);

    if (convSNPj == 0) {
      if(*trace > 1)
        Rprintf("\n  g=%d, fail to estimate phi in glmNBlog\n", g);
      
      succeed = 0;
      break;
    }
    
    /* need to substract the extra degree for freedom in the offset */
    df1 -= 1;

    if (paraDiff < fabs(*phi - phi_old)) paraDiff = fabs(*phi - phi_old);
    
    /**********************************
     * compare likelihood       
     **********************************/
    
    if((twoLL_trec1 - twoLL_trec0)/(fabs(twoLL_trec0) + 1.0) < -1e-4){
      
      nless5 = 0.0;
      for (k=0; k<N; k++) { 
        if (pY[k] < 5) nless5 += 1.0;
      }
      
      Rprintf("  g=%d, nless5=%.1f, twoLL(old, new)=(%e, %e)\n", 
              g, nless5, twoLL_trec0, twoLL_trec1);

      if (nless5 > 0.75*N) {
        succeed=0;
        break;
      }
      
      warning("\n  likelihood decreases during update of phi in glmNBlog\n");
      
      succeed = 0;
      break;
    }
    
    phi_old = *phi;
    twoLL_trec0 = twoLL_trec1;
    
    if (*trace > 3){      
      Rprintf("\n  g=%d, parDiff=%e", g, paraDiff);
      Rprintf("\n  bxj=%e, phi=%e, twoLL_trec1=%e\n", bxj_old, phi_old, twoLL_trec1);
    }
    
    /* --------------------------------------------------------
     * 3. check convergence
     * -------------------------------------------------------*/
    
    if(paraDiff < *conv){
      
      if (*trace > 1){
        Rprintf("\n  converged in glmNBlog\n");
        Rprintf("  bxj=%e, phi=%e, df1=%d, twoLL1=%e\n", bxj_old, phi_old, df1, twoLL_trec1);
      }
      
      break;
    }
  }
  
  if (g >= maxit) {
    if (*trace){
      Rprintf("\n  g=%d, paraDiff=%.3e, reach max iteration in glmNBlog", g, paraDiff);
    }
    succeed = 0;
  }
  
  if (succeed) {
    *beta = bxj_old;
    *twoLL_trec = twoLL_trec1;
    *df_resid   = df1;
    for (i=0; i<N; i++) fitted[i] = fitted2[i];    
  }else {
    *beta       = 0.0;
    *twoLL_trec = -1e6;
    *df_resid   = 0;
  }

  return(succeed);
}
