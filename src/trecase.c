/*
 *  trecase.c
 *
 *  Created by Wei Sun on 6/02/2010.
 *  Modified by Vasyl Zhabotynsky on 04/05/2011
 *
 */
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <R.h>
#include "utility.h"
#include "glm.h"
#include "ase.h"

#define LMM   5
#define NPARA 2

/**********************************************************************
 *
 * fprintRow: print out one row
 *
 **********************************************************************/

void fprintRow(FILE* fo, double beta, double chisq,  double dfr, double pval){
  
  if (chisq < 0) {
    fprintf(fo, "NA\tNA\t");
  }else{
    fprintf(fo, "%.4e\t%.3f\t", beta, chisq);
  }
  
  if (dfr < 0) {
    fprintf(fo, "NA\t");
  }else{
    fprintf(fo, "%.2f\t", dfr);
  }

  if (pval > 1.0001) {
    fprintf(fo, "NA\t");
  }else{
    fprintf(fo, "%.2e\t", pval);
  }
}

/**********************************************************************
 *
 * logL_b: joint likelihood of b_{x_j}
 *
 **********************************************************************/

double logL_b(double b1, int N, int h, int fam, double b0, double phi,
              double theta, double* y, double*x, double* mu, 
              double* mu1, double* nA, double* nTotal, double* zeta)
{
  int i, k;
  double ti, ni, ni0, mui, xi, pi; 
  double cst0, cst1, logL, sumL, piI;
  
  pi   = exp(b1)/(1.0 + exp(b1));
  cst0 = exp(b1 - b0);
  cst1 = (1.0 + exp(b1))/(1.0 + exp(b0));
  
  /* -------------------------------------------------------
   * first part, the log likelihood of TReC model
   * ------------------------------------------------------*/

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
  }
  
  /* -------------------------------------------------------
   * second part, the log likelihood of ASE model
   * ------------------------------------------------------*/
  
  sumL = 0.0;
  
  for(i=0; i<h; i++){
    ni   = nTotal[i];
    ni0  = nA[i];
    
    sumL += lchoose(ni, ni0);
    
    if(zeta[i] > 0){ piI = pi; }else { piI = 0.5; }
    
    if(ni0 > 0){
      for(k=0; k<ni0; k++) sumL += log(piI + k*theta);
    }
    
    if(ni0 < ni){
      for(k=0; k<ni-ni0; k++) sumL += log(1 - piI + k*theta);
    }
    
    for(k=0; k<ni; k++){
      sumL -= log(1 + k*theta);
    }
  }
  
  logL += sumL;
  
  return(logL);
}

/**********************************************************************
 *
 * grad_b: 1st and 2nd gradients of joint likelihood with 
 * respect to b_{x_j}
 *
 * here phi and theta are the dispersion paramters for 
 * TReC and ASE models, respectively.
 *
 **********************************************************************/

void grad_b(double b1, int N, int h, double b0, double phi,
            double theta, double* y, double*x, double* mu, 
            double* nA, double* nTotal, double* zeta, 
            double* gr, int trace) 
{
  int i, k;
  double grad1, grad2, df_dmu1, df_dmu2, dmu_db, dh_dpi1, dh_dpi2, dpi_db;
  double ti, ni, ni0, mui, xi, pi, tmp; 
  double cst0, cst1;
  

  if (trace > 9) {
    Rprintf("\n  grad_b\n");
    Rprintf("  b1=%e, b0=%e, phi=%e, theta=%e\n  mu=\n", b1, b0, phi, theta);
    
    for (i=0; i<5; i++) {
      Rprintf("%.2e ", mu[i]);
    }
    Rprintf("\n");    
  }

  /* ----------------------------------------------------------
   * first, calculate the 1st/2nd derivative of the TReC model
   * ---------------------------------------------------------*/
    
  pi    = exp(b1)/(1.0 + exp(b1));
  grad1 = 0.0;
  grad2 = 0.0;
  cst0  = exp(b1 - b0);
  cst1  = (1.0 + exp(b1))/(1.0 + exp(b0));
  
  for (i=0; i<N; i++) {
    xi = x[i];
    ti = y[i];
    mui = mu[i];
    
    if (ti < -1e-10) {
      error("invalid ti: %e @ grad_b\n", ti);
    }
    
    if (fabs(xi) < 0.01) {
      continue;
    }else if (fabs(xi - 2.0) < 0.01) {
      mui    = mui*cst0;
      dmu_db = mui;
    }else if (fabs(xi - 1.0) < 0.01) {
      mui    = mui*cst1;
      dmu_db = mui*pi;
    }else {
      error("invalid genotype\n");
    }
    
    if (mui < 1e-8) {
      Rprintf("invalid mui: %e @ grad_b, set to 1e-8\n", mui);
      mui = 1e-8;
    }
    
    df_dmu1 = ti/mui - (1.0 + phi*ti)/(1.0 + phi*mui);        
    grad1  += df_dmu1*dmu_db;
    
    df_dmu2 = -ti/(mui*mui) + phi*(1.0 + phi*ti)/(1.0 + phi*mui)/(1.0 + phi*mui);        
    grad2  += (df_dmu2*dmu_db*dmu_db + df_dmu1*dmu_db);
  }
  
  if (trace > 9) {
    Rprintf("grad1=%e, grad2=%e\n", grad1, grad2);
  }
  
  /* ----------------------------------------------------------
   * next, calculate the 1st/2nd derivative of the ASE model
   * ---------------------------------------------------------*/
  
  dh_dpi1 = 0.0;
  dh_dpi2 = 0.0;

  for (i=0; i<h; i++) {
    ni  = nTotal[i];
    ni0 = nA[i];
    
    if (zeta[i] > 0) {
      if(ni0 > 0){
        for(k=0; k<ni0; k++){
          tmp = 1.0/(pi + k*theta);
          dh_dpi1 += tmp;
          dh_dpi2 -= tmp*tmp;
        }
      }
      
      if(ni0 < ni){
        for(k=0; k<ni-ni0; k++){
          tmp = 1.0/(1.0 - pi + k*theta);
          dh_dpi1 -= tmp;
          dh_dpi2 -= tmp*tmp;
        }
      }
    }
  }
  
  if (trace > 9) {
    Rprintf("dh_dpi1=%e, dh_dpi2=%e\n", dh_dpi1, dh_dpi2);
  }
  
  /* ----------------------------------------------------------
   * add up the 1st/2nd derivatives from the TReC and ASE models
   * ---------------------------------------------------------*/
  
  dpi_db = pi/(1.0 + exp(b1));
  
  tmp    = dpi_db/(1.0 + exp(b1)) - pi*dpi_db;
  grad1 += dh_dpi1*dpi_db;
  grad2 += dh_dpi1*tmp + dh_dpi2*dpi_db;
  
  if (trace > 9) {
    Rprintf("grad1=%e, grad2=%e\n", grad1, grad2);
  }
  
  gr[0] = grad1;
  gr[1] = grad2;
}

/**********************************************************************
 * b_TReC_ml: find the MLE of b in TReC model by 
 * a Newton-Raphson algrorithm
 **********************************************************************/

void b_ml(double* b_xj, int N, int h, double b0, double phi, 
          double theta, double* y, double* x, double* mu, 
          double* nA, double* nTotal, double* zeta, 
          int limit, double eps, int trace, int* failR)
{
  int it, fail;
  double b1, del, gr[2];
  double min_b = -1e5;
  double max_b =  1e5;
  
  it   = 0;
  del  = 1.0;
  fail = 0;
  b1   = b0;
  gr[0] = 1.0;
  gr[1] = 1.0;
  
  if(trace > 3){
    Rprintf("b_ml: limit=%d, eps=%e\n", limit, eps);
  }
  
  while (it < limit && fabs(del) > eps) {
   
    grad_b(b1, N, h, b0, phi, theta, y, x, mu, nA, nTotal, zeta, gr, trace);
    
    del = gr[0]/gr[1];
    b1 -= del;
    it += 1;
    
    if (it > 10 && fabs(del) > 0.1) {
      fail = 1;
      break;
    }
    
    if(trace > 5){
      Rprintf("      b_ml: it=%d, b1=%e, del=%e, gr[0]=%e, gr[1]=%e\n", 
              it, b1, del, gr[0], gr[1]);
    }
    
    if (b1 > max_b) {
      fail = 1;
      b1 = max_b;
      if(trace > 3)
        Rprintf("    Estimate of b_xj is truncated at %.2e\n", max_b);
      break;
    }
    
    if(b1 < min_b) {
      fail = 1;
      b1 = min_b;
      if(trace > 3)
        Rprintf("    Estimate of b_xj is truncated at %.2e\n", min_b);
      break;
    }  
  }
  
  
  if (it == limit) {
    // fail = 1;
    
    if(trace > 3)
      Rprintf("  b_ml: iteration limit reached in b_ml\n");
  }
  
  if (fabs(gr[0]) > 0.01) {
    fail = 1;
    
    if(trace > 1)
      Rprintf("  b_ml: 1st derivative = %.2e\n", gr[0]);
  }
  
  if (gr[1] > 1e-16) {
    fail = 1;
    
    if(trace > 1)
      Rprintf("  b_ml: 2st derivative = %.2e\n", gr[1]);
  }
  
  *b_xj  = b1;
  *failR = fail;
}


/**********************************************************************
 *
 * trecase
 *
 * map eQTL by joinly modeling TReC and ASE. 
 *
 **********************************************************************/

void trecase (int* dims, double* Y, double* X, double* Z, double* z1, 
              double* Y1, double* Y2, double* Zh, double* offset, 
              char** output, double* RP_cut, int* cis_only,  
              int* cis_distance, int* eChr, int* ePos, int* mChr,  
              int* mPos, double* conv, double* convGLM, 
              int* yFailBaselineModel, double* scoreTestP, 
              double* transTestP, int* trace, int* succeed)
{
  int i, j, k, g, nIter, convBase, convSNPj, h0, h1;
  int df0, df1, df1_join;
  double dfr_TReC, dfr_ASE, dfr_Joint, pZk, nless5;
  double dfr_joint_ASE, dfr_joint_TReC;
  double phi0=0.0, phi1=0.0, scale=1.0, nT1;
  double twoLL0, twoLL1, twoLL_bxj0, twoLL_bxj1, btmp, paraDiff;
  double twoLL_ase0, twoLL_ase1, twoLL_trec0=0.0, twoLL_trec1=0.0;
  double twoLL_trec_joint0, twoLL_ase_joint0;
  double twoLL_trec_joint1, twoLL_ase_joint1;
  double th0=0.0; /* theta0 estimated by H0 */
  double *exPara, *nA, *nTotal, *zeta, *offsetN;
  int useASE, useTReC, useASE_j, useTReC_j, useJointModel, adjZ;
  double gr[2];
  
  int family = NB;   /* family will be decided by baseline model later */
  int linkR  = LOG;  /* link function */

  /* parameters used to estimate b_{x_j} */
  double bxj_trec, bxj_ase, pi_old, pi_ase, theta_ase;

  /* parameters used to in iterative updating */
  double bxj_old, bxj_new, theta_old, theta_new, phi_old, phi_new;
  
  /* pointers to Y, the last column of X, Z, and Zh */
  double *pY, *pY1, *pY2, *pXlast, *pZ, *pZh;
  
  /* Work array */
  double *Xb, *fitted0, *fitted1, *fitted2, *resid, *weights;
  
  /* rank of model matrix */
  int rank0, rank1;
    
  /* position difference between gene and marker */
  int pos_diff;

  /* grid used to output frequency */
  double grid;

  /* para for optimizaton of ASE model */
  double para_ase[2];

  /* summary statistics */
  double chisqTrans, chisqTReC, chisqASE, chisqJoint;
  double pvalTrans, pvalTReC, pvalASE, pvalJoint;
  
  /* p-value frequencies */
  unsigned long freqTReC[101];
  for(i=0; i<=100; i++){ freqTReC[i] = 0;  }
  
  unsigned long freqASE[101];
  for(i=0; i<=100; i++){ freqASE[i] = 0;   }

  unsigned long freqJoint[101];
  for(i=0; i<=100; i++){ freqJoint[i] = 0; }

  /* **********************************************************
   * parameters for function lbfgsb, which will be used to
   * obtain MLE for H1: with allelic imbalance
   * **********************************************************/
  
  int npara, lmm, fail, fncount, grcount, maxit1, nREPORT;
  int nbd[2];


  npara   = NPARA;
  lmm     = LMM;
  fail    = 0;
  fncount = 0;
  grcount = 0;
  maxit1  = 100;
  nREPORT = 5;
  nbd[0]  = 2;
  nbd[1]  = 2;
  //technical parameters below:
  double *wa, *g1;
  int *iwa;
  SEXP x1;
  PROTECT(x1 = allocVector(REALSXP, npara));
  //consider replacing with simple Calloc
  wa  = (double *) S_alloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm,sizeof(double));
  iwa = (int*) R_alloc(3*npara,sizeof(int));
  g1 = (double *)R_alloc(npara, sizeof(double));
  //g1 = vect(npara);
  //end technical parameters
  
  double initPara[2];
  double lower[2];
  double upper[2];
  double Fmin, factr, pgtol;
  
  /* initPara = c(theta, pi) */
  
  initPara[0] = 0.1; 
  initPara[1] = 0.5;
  
  lower[0] = 0.0 + 1e-8;
  upper[0] = 1e8;
  
  lower[1] = 0.0 + 1e-8;
  upper[1] = 1.0 - 1e-8;
  
  factr = 1e7;
  pgtol = 0.0;
  
  char msg[1023];
  
  /**********************************************************/

  int nY = dims[0];
  int nX = dims[1];
  int nZ = dims[2];
  int N  = dims[3];
  int maxit = dims[4];
  int useOffset = dims[5];
  
  /* minimum of total number of reads */
  int min_nT    = dims[6];
  
  /* minimum # of samples with >= min_nT reads */  
  int min_N     = dims[7];
  
  /* minimum # of sample with heterzygous genotypes */  
  int min_Nhet  = dims[8];
  
  /* p-value cutoff */  
  double P_cut  = *RP_cut;

  /* ======================================================== */

  exPara = (double *) Calloc(3*N+7, double);
  
  exPara[0] = (double) N; 
  exPara[1] = 0.0; /* h: sample size for ASE model */
  exPara[2] = 0.0; /* family, which will be updated in baseline GLM */
  exPara[3] = 0.0; /* bxj_old */
  exPara[4] = 0.0; /* phi_old */
  exPara[5] = 0.0; /* theta   */
  exPara[6] = 0.0; /* pi      */
  
  nA     = exPara + 7;
  nTotal = nA + N;
  zeta   = nTotal + N;
  
  /* initial zeta */
  for (k=0; k<N; k++) {
    zeta[k] = -1.0;
  }
  
  /* ======================================================== */
  
  /* dimsNew is passed to function glmNB, and then function glmFit */
  int dimsNew[5];
  dimsNew[0] = N;
  dimsNew[1] = nX;
  dimsNew[2] = maxit;
  dimsNew[3] = 0; /* whetehr to use initial values */
  dimsNew[4] = useOffset;
  
  /* allocate memory. The extra column in Xb is used to store 
   * one SNP, the same as for X
   */

  Xb      = (double *) Calloc(N*(nX+1), double);
  fitted0 = (double *) Calloc(N, double); /* for Null TReC model */
  fitted1 = (double *) Calloc(N, double); /* for Alternative TReC model */
  fitted2 = (double *) Calloc(N, double); /* for Alternative TReC model */
  resid   = (double *) Calloc(N, double);
  weights = (double *) Calloc(N, double);
  offsetN = (double *) Calloc(N, double);
  
  /* point to the last column of X */
  pXlast  = X + N*nX;
  
  if(*trace){
    Rprintf("\n--------------------------------------------------\n");
    Rprintf("(nY, nZ, nX, N, maxit) = (%d, %d, %d, %d, %d)\n", nY, nZ, nX, N, maxit);
    Rprintf("(useOffset, min_nT, min_N, min_Nhet) = (%d, %d, %d, %d)\n", 
            useOffset, min_nT, min_N, min_Nhet);
    Rprintf("P_cut=%e", P_cut);
    Rprintf("\n--------------------------------------------------\n");
  }
  
  /* output file handles */
  FILE *fo, *ff;
  
  /* time records */
  time_t sec_s;
  time_t sec_e;
  
  /* starting time */
  sec_s = time(NULL);
  
  /* output file for the eQTL mapping results */
  fo = fopen (output[0], "w");
  
  /**
   * write out the header in the output file
   */
//  fprintf(fo, "GeneRowID\tMarkerRowID\tTReC_b\tTReC_Chisq\tTReC_df\tTReC_Pvalue\t");
//  fprintf(fo, "ASE_b\tASE_Chisq\tASE_df\tASE_Pvalue\t");
//  fprintf(fo, "Joint_b\tJoint_Chisq\tJoint_df\tJoint_Pvalue\t");
//  fprintf(fo, "n_TReC\tn_ASE\tn_ASE_Het\ttrans_Chisq\ttrans_Pvalue\tfinal_Pvalue\n");
  fprintf(fo, "GeneRowID\tMarkerRowID\t");
  fprintf(fo, "NBod\tBBod\t");
  fprintf(fo, "TReC_b\tTReC_Chisq\tTReC_df\tTReC_Pvalue\t");
  fprintf(fo, "ASE_b\tASE_Chisq\tASE_df\tASE_Pvalue\t");
  fprintf(fo, "Joint_b\tJoint_Chisq\tJoint_df\tJoint_Pvalue\t");
  fprintf(fo, "n_TReC\tn_ASE\tn_ASE_Het\ttrans_Chisq\ttrans_Pvalue\tfinal_Pvalue\n");

  /* **********************************************************
   * identifify eQTL gene by gene
   * **********************************************************/
  
  /* pY/pY1/pY2 are pointers to gene expression data */
  pY  = Y;
  pY1 = Y1;
  pY2 = Y2;

  for(i=0; i<nY; i++,pY+=N,pY1+=N,pY2+=N){
    
    useASE  = 1;
    useTReC = 1;
    adjZ    = 1;
    
    if(*trace){
      Rprintf("\ni=%d\n", i);
    }

    /* **********************************************************
     * fit a baseline model using only the confouding covariates 
     * family is assigned to a value of either Poisson or NB
     * **********************************************************/
    
    dimsNew[1] = nX;
    dimsNew[3] = 0; /* no initial values. NOTE, dimsNew[3] will be 
                       changed within glmNB, dumb. I know */

    convBase = glmNB(dimsNew, &nIter, pY, z1, &linkR, offset, X, convGLM, 
                     &rank0, Xb, fitted0, resid, weights, &phi0, &scale, 
                     &df0, &family, &twoLL_trec0, scoreTestP, trace, &btmp);
    
    if (!convBase){
      useTReC = 0;
      if (*trace)
        Rprintf("  i=%d, fail to fit baseline TReC model\n", i);
    }
    
    /* ======================================================== */
    exPara[2] = (double)family;
    /* ======================================================== */

    /* *****************************************************
     * calculate nA and nTotal
     * *****************************************************/
    
    h0 = 0;
    
    for (k=0; k<N; k++) {
      nT1 = pY1[k] + pY2[k];
      
      if(nT1 < min_nT){ continue; }
      
      nTotal[h0] = nT1;
      nA[h0]     = pY1[k];
      
      h0++;
    }
    
    /* if sample size is not enough */
    if(h0 < min_N) useASE = 0;
    
    if(useASE){
      /* *****************************************************
       * obtain MLE for H0: no allelic imbalance
       * *****************************************************/
      
      exPara[1]   = (double) h0;
      exPara[6]   = 0.5; /* fixed pi */
      
      initPara[0] = 0.1; 
      initPara[1] = 0.5; /* value for pi, not used for H0 */
      npara       = 1;
      
      // Rprintf("trecase: i=%d, h0=%d, nbd[0]=%d, range=(%.2e,%.2e), starting lbfgsb\n", 
      //        i, h0, nbd[0], lower[0], upper[0]);
	  /* *****************************************************
       * use slight modification of lbfgsb - lbfgsb1: passes two arrays
	   * instead of recreation of them at each iteration
       * *****************************************************/
      /*lbfgsb(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
      *       negLogH0, negGradLogH0, &fail, (void*)exPara, factr, pgtol,  
      *       &fncount, &grcount, maxit, msg, 0, nREPORT);*/

      lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
             negLogH0, negGradLogH0, &fail, (void*)exPara, factr, pgtol,  
             &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,x1);
	  
      
      if (fail) {
        if (*trace)
          Rprintf("  i=%d, fail to fit baseline ASE model\n", i);

        useASE = 0;
        th0    = -1000.00;
        twoLL_ase0 = 1e16;
      }else{
        th0 = initPara[0];
        twoLL_ase0 = -2.0*Fmin;

        if(*trace > 1){
          negGradLogH0(npara, initPara, gr, (void*)exPara,x1);
          Rprintf("\n  Obtained MLE for H0: twologlik0=%.4e, ", twoLL_ase0);
          Rprintf("theta0=%.4e, gradience=%.4e\n", th0, *gr);
        }
      }
    }
    
    yFailBaselineModel[i] = 1 - useTReC + 2*(1 - useASE);

    if (*trace > 1)
      Rprintf("  i=%d, h0=%d, useTReC=%d, useASE=%d, phi0=%e, theta0=%e\n", 
              i, h0, useTReC, useASE, phi0, th0);
    
    /* **********************************************************
     * Now start to fit models with both confounding covariates
     * and each of the SNPs 
     * **********************************************************/

    pZ  = Z;
    pZh = Zh;
    
    for(j=0; j<nZ; j++,pZ+=N,pZh+=N){
      
      adjZ = 1;

      if(*cis_only){
        if(eChr[i] != mChr[j]) continue;
        
        pos_diff = abs(ePos[i] - mPos[j]);
        
        if(pos_diff > *cis_distance) continue;
      }
      
      if(*trace > 2){
        Rprintf("\ni=%d, j=%d\n", i, j);
      }
      
      useTReC_j = useTReC;
      useASE_j  = useASE;

      h1 = 0;
      
      /* **********************************************************
       * Initial fitting with TReC and ASE
       * **********************************************************/
      
      if (useTReC){
        /* *
         * fill the last column of X by Z[j], genotype of one SNP
         */
        
        for (k=0; k<N; k++) {
          pXlast[k]  = pZ[k];
          fitted1[k] = fitted0[k];
        }
        
        phi1  = phi0;
        scale = 1.0;
        
        dimsNew[1] = nX;
        dimsNew[3] = 1; /* use initial values */
                
        convSNPj = glmNBlog(dimsNew, &nIter, pY, z1, &linkR, offset, X, conv, convGLM, 
                            &rank1, Xb, fitted1, resid, weights, &phi1, &scale, 
                            &df1, &family, &twoLL_trec1, scoreTestP, trace, 
                            &bxj_trec, fitted2, offsetN);
        
        if(*trace > 1){
          Rprintf("\n  convSNPj@glmNBlog = %d\n", convSNPj);
        }
        
        if(convSNPj == 0){
          
          for (k=0; k<N; k++) {
            fitted1[k] = fitted0[k];
          }
          
          adjZ = 0;
          dimsNew[1] = nX + 1;
          dimsNew[3] = 1; /* use initial values */

          convSNPj = glmNB(dimsNew, &nIter, pY, z1, &linkR, offset, X, convGLM, 
                           &rank1, Xb, fitted1, resid, weights, &phi1, &scale, 
                           &df1, &family, &twoLL_trec1, scoreTestP, trace, &bxj_trec);
          if(*trace){
            Rprintf("\n  convSNPj@glmNB = %d\n", convSNPj);
          }
          
        }
        
        if(convSNPj == 0){
          if(*trace){
            Rprintf("\n  Fail TReC: i=%d, j=%d, family=%d\n", i, j, family);
          }
          useTReC_j = 0;
        }else if(df0 - df1 != 1){
          if(*trace){
            Rprintf("\n  i=%d, j=%d, dfs=(%d, %d), ranks=(%d, %d)\n", 
                  i, j, df0, df1, rank0, rank1);
          }
          if(df0 - df1 < 0.5) useTReC_j = 0;
        }
        
        if (useTReC_j && *trace > 1) {
          Rprintf("\n  finish TReC: i=%d, j=%d, phi=%e, b=%e, ", i, j, phi1, bxj_trec);
          Rprintf("twologlike=%e\n", twoLL_trec1);
        }
        
      }
      
      if (useASE) {
        /* *****************************************************
         * calculate nA and nTotal
         * *****************************************************/
        h0 = 0;
        h1 = 0;
        
        for (k=0; k<N; k++) {
          nT1 = pY1[k] + pY2[k];
          
          if(nT1 < min_nT) continue;
          
          if(nTotal[h0] != nT1) error("mismatch ;( \n");
          
          /* pZ[k] = 0 or 4 if homozygous */
          if (fabs(pZh[k] - 2.0) > 1.99) {
            zeta[h0] = -1.0;
            nA[h0]   = pY2[k];
          }else {
            zeta[h0] = 1.0;
            
            if(fabs(pZh[k] - 1.0) < 0.01){
              nA[h0] = pY2[k];
            }else if(fabs(pZh[k] - 3.0) < 0.01){
              nA[h0] = pY1[k];
            }else {
              error("invalid values for Z\n");
            }
            
            h1++;
          }
          h0 ++;
        }
        
        /* if sample size of heterzygous genotype is not enough */
        if(h1 < min_Nhet){
          useASE_j=0;
        }else{
          /* *****************************************************
           * obtain MLE for H1: with allelic imbalance
           * *****************************************************/
          
          exPara[1]   = (double) h0;
          initPara[0] = th0; 
          initPara[1] = 0.5;
          npara       = 2;

          // Rprintf("trecase: i=%d, j=%d, h0=%d, th0=%e, nbd=(%d, %d) starting lbfgsb\n", 
          //        i, j, h0, th0, nbd[0], nbd[1]);

          //lbfgsb(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
          //       negLogH1, negGradLogH1, &fail, (void*)exPara, factr, pgtol,  
          //       &fncount, &grcount, maxit1, msg, 0, nREPORT);
          lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
             negLogH1, negGradLogH1, &fail, (void*)exPara, factr, pgtol,  
             &fncount, &grcount, maxit1, msg, 0, nREPORT, wa, iwa, g1,x1);
          
          if (fail){
            useASE_j=0;
            
            if (*trace)
              Rprintf("  i=%d, j=%d, fail ASE model\n", i, j);

          }else{
            twoLL_ase1 = -2.0*Fmin;
            
            theta_ase = initPara[0];
            pi_ase    = initPara[1];
            bxj_ase   = log(pi_ase) - log(1.0-pi_ase);
          }
        }
        
        if (*trace > 2){
          Rprintf("  i=%d, j=%d, h0=%d, h1=%d, useASE_j=%d\n", 
                  i, j, h0, h1, useASE_j);
          
          if(useASE_j)
            Rprintf("  pi_ase=%e, theta_ase=%e, bxj_ase=%e\n", 
                    pi_ase, theta_ase, bxj_ase);
        }
      }
      
      /* **********************************************************
       * Iterative updating
       * **********************************************************/
      
      useJointModel = 0;

      if(useTReC_j && useASE_j && adjZ){
        
        useJointModel = 1;
        
        bxj_old   = 0.0;
        theta_old = th0;
        phi_old   = phi0;
        pi_old    = 0.5;
        
        if (*trace > 2)
          Rprintf("\n  i=%d, j=%d, use joint models\n\n", i, j);
        
        for (k=0; k<N; k++) {
          fitted2[k] = fitted0[k];
        }
                
        /**********************************
         * calculate old likelihood       
         **********************************/
        if (family==NB) {
          twoLL_trec_joint0 = 2.0*loglik_NB(N, phi_old, fitted2, pY);
        }else if (family==POISSON) {
          twoLL_trec_joint0 = 2.0*loglik_Poisson(N, fitted2, pY);
        }else {
          error("invalid family\n");
        }
        
        para_ase[0] = theta_old;
        para_ase[1] = pi_old;
        
        twoLL_ase_joint0 = - 2.0*negLogH1(2, para_ase, exPara,x1);
        
        twoLL0 = twoLL_trec_joint0 + twoLL_ase_joint0;
        twoLL1 = twoLL0;
        
        /* --------------------------------------------------------
         * iterations to estimate bxj
         * -------------------------------------------------------*/
        
        for(g=0; g<maxit; g++){
          
          paraDiff = 0.0;

          /* --------------------------------------------------------
           * 1. Given theta, phi, b_0, b_k, and b_u, estimate b_{x_j}
           * -------------------------------------------------------*/

          /*******************************************************
           * estimate b_xj      
           *******************************************************/
          twoLL_bxj0  = twoLL1;

          fail = 0;
                    
          b_ml(&bxj_new, N, h0, bxj_old, phi_old, theta_old, pY, 
               pZ, fitted2, nA, nTotal, zeta, maxit, *convGLM, 
               *trace, &fail);
          
          if(fail){
            if(*trace){
              Rprintf("\n  i=%d, j=%d, fail to estimate bxj in b_ml (trecase)\n", i, j);
            }
            
            useJointModel = 0;
            break;
          }
                    
          if (paraDiff < fabs(bxj_new - bxj_old)) 
            paraDiff = fabs(bxj_new - bxj_old);
          
          /********************************************************
           * calculate new likelihood according to new bxj estimate     
           *******************************************************/
          pi_old = exp(bxj_new)/(1.0 + exp(bxj_new));
          
          for (k=0; k<N; k++){
            pZk = pZ[k];
            
            if (fabs(pZk - 1.0) < 0.01) {
              fitted2[k] = exp(log(fitted2[k]) - log(1.0 + exp(bxj_old)) + log(1.0 + exp(bxj_new)));
            }else if (fabs(pZk - 2.0) < 0.01) {
              fitted2[k] = exp(log(fitted2[k]) + (bxj_new - bxj_old));
            }
          }
          
          if (family==NB) {
            twoLL_trec_joint1 = 2.0*loglik_NB(N, phi_old, fitted2, pY);
          }else if (family==POISSON) {
            twoLL_trec_joint1 = 2.0*loglik_Poisson(N, fitted2, pY);
          }else {
            error("invalid family\n");
          }

          para_ase[0] = theta_old;
          para_ase[1] = pi_old;
          
          twoLL_ase_joint1 = - 2.0*negLogH1(2, para_ase, exPara,x1);

          twoLL_bxj1 = twoLL_trec_joint1 + twoLL_ase_joint1;
          
          if (*trace > 1) {
            Rprintf("\n  i=%d, j=%d, g=%d\n", i, j, g);
            Rprintf("  theta=%e, phi=%e, b_old=%e, b_new=%e\n", 
                    theta_old, phi_old, bxj_old, bxj_new);
            
            Rprintf("  twoLL(old, new)=(%e, %e)\n", twoLL_bxj0, twoLL_bxj1);
          }
          
          if((twoLL_bxj1 - twoLL_bxj0)/fabs(twoLL_bxj0) < -0.01){
            Rprintf("\n  i=%d, j=%d, g=%d\n", i, j, g);
            error("  likelihood decreases during update of b_xj\n");
          }
          
          bxj_old = bxj_new;
          
          /* --------------------------------------------------------
           * 2. Given b_{x_j}, estimate theta
           * -------------------------------------------------------*/
          twoLL_ase_joint0 = twoLL_ase_joint1;

          exPara[1]   = (double) h0;
          exPara[6]   = pi_old; /* fixed pi */
          initPara[0] = theta_old; 
          initPara[1] = pi_old; /* value for pi, not used for H0 */
          npara       = 1;
          
          // Rprintf("trecase: i=%d, j=%d, g=%d, h0=%d, theta=%e, pi=%e, nbd=(%d, %d) starting lbfgsb\n", 
          //        i, j, g, h0, theta_old, pi_old, nbd[0], nbd[1]);

          //lbfgsb(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
          //       negLogH0, negGradLogH0, &fail, (void*)exPara, factr, pgtol,  
          //       &fncount, &grcount, maxit, msg, 0, nREPORT);
          lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
             negLogH0, negGradLogH0, &fail, (void*)exPara, factr, pgtol,  
             &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,x1);

          twoLL_ase_joint1 = -2.0*Fmin;
                    
          if (fail) {
            if(*trace){
              Rprintf("\n  i=%d, j=%d, h0=%d, fail to estimate theta in joint model\n", 
                      i, j, h0);
              negGradLogH0(npara, initPara, gr, (void*)exPara,x1);
              Rprintf("  theta=%.4e, gradience=%.4e, fail=%d\n", initPara[0], gr[0], fail);
            }
            useJointModel = 0;
            break;
          }
                    
          theta_new = initPara[0];
          
          if (paraDiff < fabs(theta_new - theta_old)) 
          	paraDiff = fabs(theta_new - theta_old);
          
          /**********************************
           * compare old and new likelihood       
           **********************************/
          
          if((twoLL_ase_joint1 - twoLL_ase_joint0)/fabs(twoLL_ase_joint0) < -0.01 ){
            Rprintf("\n  i=%d, j=%d, g=%d ", i, j, g);
            Rprintf("  theta_old=%e, theta_new=%e\n", theta_old, theta_new);
            Rprintf("  twoLL(old, new)=(%e, %e)", twoLL_ase_joint0, twoLL_ase_joint1);
            error("\n  likelihood decreases during update of theta in join modeling\n");
          }
          
          theta_old = theta_new;
          
          /* --------------------------------------------------------
           * 3. Given b_{x_j} and phi, estimate b_0, b_k, b_u, phi
           * in fact, we do not estimate b_0, b_k, and b_u, 
           * instead we only need to estimate mu, or fitted1.
           * -------------------------------------------------------*/
          twoLL_trec_joint0 = twoLL_trec_joint1;

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
          
          phi_new  = phi_old;
          convSNPj = glmNB(dimsNew, &nIter, pY, z1, &linkR, offsetN, X, convGLM, 
                           &rank1, Xb, fitted2, resid, weights, &phi_new, &scale, 
                           &df1_join, &family, &twoLL_trec_joint1, scoreTestP, trace, 
                           &btmp);
          
          if (convSNPj == 0) {
            if(*trace)
              Rprintf("\n  i=%d, j=%d, fail to estimate phi in joint model\n", i, j);
            
            useJointModel = 0;
            break;
          }
          
          /* need to substract the extra degree for freedom in the offset */
          df1_join -= 1;
          
          if (paraDiff < fabs(phi_new - phi_old)) paraDiff = fabs(phi_new - phi_old);
          
          /**********************************
           * compare likelihood       
           **********************************/
          
          if((twoLL_trec_joint1 - twoLL_trec_joint0)/fabs(twoLL_trec_joint0) < -0.01){
            Rprintf("\n  i=%d, j=%d, g=%d ", i, j, g);
            Rprintf("  twoLL(old, new)=(%e, %e)", twoLL_trec_joint0, twoLL_trec_joint1);
            error("\n  likelihood decreases during update of phi in join modeling\n");
          }
                    
          phi_old = phi_new;
          twoLL_trec_joint0 = twoLL_trec_joint1;
          
          if (*trace > 1){
            Rprintf("\n  i=%d, j=%d, g=%d, parDiff=%e", i, j, g, paraDiff);
            Rprintf("\n  bxj=%e, theta=%e, phi=%e\n", bxj_old, theta_old, phi_old);
          }
          
          /* --------------------------------------------------------
           * check convergence
           * -------------------------------------------------------*/
          
          if(paraDiff < *conv){
            twoLL1  = twoLL_trec_joint1 + twoLL_ase_joint1;
            
            if (*trace > 1){
              Rprintf("\n  i=%d, j=%d, converged using joint model\n", i, j);
              Rprintf("  bxj=%e, theta=%e, phi=%e, twoLL0=%e, twoLL1=%e\n", 
                      bxj_old, theta_old, phi_old, twoLL0, twoLL1);
            }
            
            break;
          }
        }
        
        if (g >= maxit) {
          if (*trace > 1){
            Rprintf("\n  i=%d, j=%d, reach max iteration using joint model. ", i, j);
            Rprintf("g=%d, maxit=%d, paraDiff=%.3e\n", g, maxit, paraDiff);
          }
          useJointModel = 0;
        }
              
      }
      
      /* *********************************************************
       * calculate p-value for TReC model
       * *********************************************************/
      
      if (useTReC_j) {
        chisqTReC = twoLL_trec1 - twoLL_trec0;
        
        if (chisqTReC < -0.01) {
          nless5 = 0.0;
          for (k=0; k<N; k++) { 
            if (pY[k] < 5) nless5 += 1.0;
          }
          
          Rprintf("  g=%d, nless5=%f, N=%d, 2logL(old, new)=(%e, %e)\n", 
                  g, nless5, N, twoLL_trec0, twoLL_trec1);
          
          useTReC_j=0;
          
          if (nless5 <= 0.75*N) {
            error("likelihood decreases for TReC model: i=%d, j=%d, 2logL=(%.2e, %.2e), chisq=%.3e\n", 
                  i, j, twoLL_trec0, twoLL_trec1, chisqTReC);
          }
          
        }else{
        
          dfr_TReC = (double)(df0 - df1);

          if(dfr_TReC < 0.5){
            pvalTReC = 995.0;
            if (*trace > 1) {
              Rprintf("dfr_TReC=%e, df0=%d, df1=%d\n", dfr_TReC, df0, df1);
            }
          }else{
          
            if (chisqTReC < 1e-5) { 
              pvalTReC = 1.0; 
            }else{
              pvalTReC  = pchisq(chisqTReC, dfr_TReC, 0, 0);
            }

            k = (int) (pvalTReC / 0.01);
            freqTReC[k] += 1;
          }
          
        }
      }
      
      if(! useTReC_j){
        useJointModel = 0;
        chisqTReC = -995.0;
        dfr_TReC  = -1.0;
        pvalTReC  =  995.0;
      }
      
      /* *********************************************************
       * calculate p-value for ASE model
       * *********************************************************/
      
      if(useASE_j) {
        chisqASE = twoLL_ase1 - twoLL_ase0;
        
        if (chisqASE < -0.01) {
          error("likelihood decreases for ASE model: i=%d, j=%d, twoLL=(%.2e, %.2e), chisq=%.3e\n", 
                i, j, twoLL_ase0, twoLL_ase1, chisqASE);
        }
        
        dfr_ASE = 1.0;
        if (fabs(th0) < 1e-7)  dfr_ASE += 1.0;
        if (fabs(theta_ase) < 1e-7)  dfr_ASE -= 1.0;
        
        if(dfr_ASE < 0.5){
          pvalASE = 995.0;
          if (*trace > 1) {
            Rprintf("dfr_ASE=%e, df0=%d, df1=%d\n", dfr_ASE, df0, df1);
          }          
        }else{
          if (chisqASE < 1e-5) { 
            pvalASE = 1.0; 
          }else{
            pvalASE = pchisq(chisqASE, dfr_ASE, 0, 0);
          }
          
          k  = (int) (pvalASE / 0.01);
          freqASE[k] += 1;
        }
        
      }else {
        chisqASE = -995.0;
        dfr_ASE  = -1.0;
        pvalASE  =  995.0;
      }
      
      /* *********************************************************
       * calculate p-value for joint model
       * *********************************************************/
      
      if(useJointModel){

        chisqTrans = twoLL_ase1 + twoLL_trec1 - twoLL1;
        
        if (chisqTrans < -0.01) {
          Rprintf("twoLL(ASE, TReC, TReCASE)=(%.3e, %.3e)\n",  twoLL_ase1, twoLL_ase1, twoLL1);
          warning("testing for cis v.s. trans, i=%d, j=%d, chisq=%.3e\n", i, j, chisqTrans);
        }
        
        if (chisqTrans < 1e-5) { 
          pvalTrans = 1.0; 
        }else{
          pvalTrans = pchisq(chisqTrans, 1.0, 0, 0);
        }
        
        chisqJoint = twoLL1 - twoLL0;

        if (chisqJoint < -0.01) {
          Rprintf("twoLL_trec=(%.3e, %.3e) ",  twoLL_trec0, twoLL_trec_joint1);
          Rprintf("twoLL_ase=(%.3e, %.3e) \n", twoLL_ase0,  twoLL_ase_joint1);
          
          error("wrong twoLL for joint model i=%d, j=%d, chisq=%.3e\n", 
                i, j, chisqJoint);
        }

        dfr_joint_ASE = 1.0;
        if (fabs(th0) < 1e-7)  dfr_joint_ASE += 1.0;
        if (fabs(theta_new) < 1e-7)  dfr_joint_ASE -= 1.0;
        
        dfr_joint_TReC = (double)(df0 - df1_join);
        
        /* we need to subtrace 1 here because in joint model, 
         * b (of TReC model) and pi (of ASE model) are the same */
        dfr_Joint = dfr_joint_ASE + dfr_joint_TReC - 1.0;
        
        if(dfr_Joint < 0.5){
          pvalJoint = 995.0;
          if (*trace > 1) {
            Rprintf("dfr_Joint=%e, dfr_joint_ASE=%e, df0=%d, df1_join=%d\n", 
                    dfr_Joint, dfr_joint_ASE, df0, df1);
          }          
        }else{
          if (chisqJoint < 1e-5) { 
            pvalJoint = 1.0; 
          }else{
            pvalJoint = pchisq(chisqJoint, dfr_Joint, 0, 0);
          }
          
          k  = (int) (pvalJoint / 0.01);
          freqJoint[k] += 1;
        }
        
      }else {
        chisqJoint = -995.0;
        dfr_Joint  = -1.0;
        pvalJoint  =  995.0;
        
        chisqTrans = -995.0;
        pvalTrans  =  995.0;        
      }
      
      if(pvalTReC < P_cut || pvalASE < P_cut || pvalJoint < P_cut){
        
        if (pvalTReC  < P_cut) freqTReC[100] += 1;
        
        if (pvalASE   < P_cut) freqASE[100]  += 1;

        if (pvalJoint < P_cut) freqJoint[100] += 1;

        /* gene ID and SNP ID */
//        fprintf(fo, "%d\t%d\t", i+1, j+1);
//        fprintRow(fo, bxj_trec, chisqTReC,  dfr_TReC,  pvalTReC);
//        fprintRow(fo, bxj_ase,  chisqASE,   dfr_ASE,   pvalASE);
//        fprintRow(fo, bxj_old,  chisqJoint, dfr_Joint, pvalJoint);
//        fprintf(fo, "%d\t%d\t%d\t", df0+nX+1, h0, h1);
        fprintf(fo, "%d\t%d\t", i+1, j+1);
        fprintf(fo, "%.2e\t%.2e\t", phi_new, theta_new);
        fprintRow(fo, bxj_trec, chisqTReC,  dfr_TReC,  pvalTReC);
        fprintRow(fo, bxj_ase,  chisqASE,   dfr_ASE,   pvalASE);
        fprintRow(fo, bxj_old,  chisqJoint, dfr_Joint, pvalJoint);
        fprintf(fo, "%d\t%d\t%d\t", df0+nX+1, h0, h1);
        
        if (chisqTrans < 0.0) {
          fprintf(fo, "NA\t");
        }else{
          fprintf(fo, "%.3f\t", chisqTrans);
        }

        if (pvalTrans > 1.0001) {
          fprintf(fo, "NA\t");
        }else{
          fprintf(fo, "%.2e\t", pvalTrans);
        }

        if(pvalTrans < *transTestP || pvalTrans > 1.0001){
          if (pvalTReC > 1.0001) {
            fprintf(fo, "NA\n");
          }else{
            fprintf(fo, "%.2e\n", pvalTReC);
          }
        }else {
          if (pvalJoint > 1.0001) {
            fprintf(fo, "NA\n");
          }else{
            fprintf(fo, "%.2e\n", pvalJoint);
          }
        }

      }
      
    }
  }
  
  fclose(fo);
  
  /* *********************************************************
   * print out p-value frequencies of TReC model
   * *********************************************************/
  
  ff   = fopen(output[1], "w");
  grid = 0.0;
  
  fprintf(ff, "Model\t<%.2e", P_cut);
  for(i=0; i<100; i++){
    fprintf(ff, "\t%.2f-%.2f", grid, grid+0.01);
    grid += 0.01;
  }
  fprintf(ff, "\n");
  
  fprintf(ff, "TReC\t%lu", freqTReC[100]);
  for(i=0; i<100; i++){
    fprintf(ff, "\t%lu", freqTReC[i]);
  }
  fprintf(ff, "\n");
  
  fprintf(ff, "ASE\t%lu", freqASE[100]);
  for(i=0; i<100; i++){
    fprintf(ff, "\t%lu", freqASE[i]);
  }
  fprintf(ff, "\n");

  fprintf(ff, "Joint\t%lu", freqJoint[100]);
  for(i=0; i<100; i++){
    fprintf(ff, "\t%lu", freqJoint[i]);
  }
  fprintf(ff, "\n");
  
  fclose(ff);
  
  Free(Xb);
  Free(fitted0);
  Free(fitted1);
  Free(fitted2);
  Free(resid);
  Free(weights);
  Free(offsetN);
  Free(exPara);
  UNPROTECT(1);
  /* end time */
  sec_e  = time(NULL);
  if(*trace){
    Rprintf("\n------------------------------------------------------\n");
    Rprintf("total time spent is %ld secs\n", sec_e-sec_s);
    Rprintf("------------------------------------------------------\n");
  }
  
  *succeed = 1;
}

/**********************************************************************
 *
 * trecase_max1
 *
 * map eQTL by joinly modeling TReC and ASE. 
 * 
 * This function is a simplified version of trecase to be used by 
 * function trecase_permute
 *
 **********************************************************************/

void trecase_max1 (int* dims, double* Y, double* X, double* Z, 
                   double* z1, double* Y1, double* Y2, double* Zh, 
                   double* offset, int* cis_only, int* cis_distance, 
                   int* eChr, int* ePos, int* mChr, int* mPos, 
                   double* conv, double* convGLM, int* yFailBaselineModel, 
                   double* scoreTestP, double* transTestP, 
                   int* best_m, double* pval, int* trace, int* succeed, 
                   double *Xb, double *fitted0, double *fitted1, 
                   double *fitted2, double *resid, double *weights, 
                   double *offsetN, double *exPara,
				   double *wa, int *iwa, double *g1, SEXP x1)
{
  int i, j, k, g, nIter, convBase, convSNPj, h0, h1;
  int df0, df1, df1_join;
  double dfr_TReC, dfr_ASE, dfr_Joint, pZk, nless5;
  double dfr_joint_ASE, dfr_joint_TReC;
  double phi0=0.0, phi1=0.0, scale=1.0, nT1;
  double twoLL0, twoLL1, twoLL_bxj0, twoLL_bxj1, btmp, paraDiff;
  double twoLL_ase0, twoLL_ase1, twoLL_trec0=0.0, twoLL_trec1=0.0;
  double twoLL_trec_joint0, twoLL_ase_joint0;
  double twoLL_trec_joint1, twoLL_ase_joint1;
  double th0=0.0; /* theta0 estimated by H0 */
  double *nA, *nTotal, *zeta;
  int useASE, useTReC, useASE_j, useTReC_j, useJointModel, adjZ;
  double gr[2];
  
  int family = -1;   /* family will be decided by baseline model later */
  int linkR  = LOG;  /* link function */
  
  /* parameters used to estimate b_{x_j} */
  double bxj_trec, bxj_ase, pi_old, pi_ase, theta_ase;
  
  /* parameters used to in iterative updating */
  double bxj_old, bxj_new, theta_old, theta_new, phi_old, phi_new;
  
  /* pointers to Y, the last column of X, Z, and Zh */
  double *pY, *pY1, *pY2, *pXlast, *pZ, *pZh;
  
  /* rank of model matrix */
  int rank0=0, rank1=0;
  
  /* position difference between gene and marker */
  int pos_diff;
  
  /* grid used to output frequency */
  double grid;
  
  /* para for optimizaton of ASE model */
  double para_ase[2];
  
  /* summary statistics */
  double chisqTrans, chisqTReC, chisqASE, chisqJoint;
  double pvalTrans, pvalTReC, pvalASE, pvalJoint;

  /* p-value for each gene */
  double pv1=0.0, pv_min[4];
  
  /* index of the marker with best association */
  int bestM[4];
  
  /* **********************************************************
   * parameters for function lbfgsb, which will be used to
   * obtain MLE for H1: with allelic imbalance
   * **********************************************************/
  
  int npara, lmm, fail, fncount, grcount, maxit1, nREPORT;
  int nbd[2];
  
  npara   = NPARA;
  lmm     = LMM;
  fail    = 0;
  fncount = 0;
  grcount = 0;
  maxit1  = 100;
  nREPORT = 5;
  nbd[0]  = 2;
  nbd[1]  = 2;

  double initPara[2];
  double lower[2];
  double upper[2];
  double Fmin, factr, pgtol;
  
  /* initPara = c(theta, pi) */
  
  initPara[0] = 0.1; 
  initPara[1] = 0.5;
  
  lower[0] = 0.0 + 1e-8;
  upper[0] = 1e8;
  
  lower[1] = 0.0 + 1e-8;
  upper[1] = 1.0 - 1e-8;
  
  factr = 1e7;
  pgtol = 0.0;
  
  char msg[1023];
  
  /**********************************************************/
  
  int nY = dims[0];
  int nX = dims[1];
  int nZ = dims[2];
  int N  = dims[3];
  int maxit = dims[4];
  int useOffset = dims[5];
  
  /* minimum of total number of reads */
  int min_nT    = dims[6];
  
  /* minimum # of samples with >= min_nT reads */  
  int min_N     = dims[7];
  
  /* minimum # of sample with heterzygous genotypes */  
  int min_Nhet  = dims[8];
  
  /* ======================================================== */
  
  // exPara = (double *) Calloc(3*N+7, double);
  
  exPara[0] = (double) N; 
  exPara[1] = 0.0; /* h: sample size for ASE model */
  exPara[2] = 0.0; /* family, which will be updated in baseline GLM */
  exPara[3] = 0.0; /* bxj_old */
  exPara[4] = 0.0; /* phi_old */
  exPara[5] = 0.0; /* theta   */
  exPara[6] = 0.0; /* pi      */
  
  /* nA and nTotal will be initialized later */
  nA     = exPara + 7;
  nTotal = nA + N;
  zeta   = nTotal + N;
  
  /* initial zeta */
  for (k=0; k<N; k++) {
    zeta[k] = -1.0;
  }
  
  /* ======================================================== */
  
  /* dimsNew is passed to function glmNB, and then function glmFit */
  int dimsNew[5];
  dimsNew[0] = N;
  dimsNew[1] = nX;
  dimsNew[2] = maxit;
  dimsNew[3] = 0; /* whetehr to use initial values */
  dimsNew[4] = useOffset;
      
  /* pointer to the last column of X */
  pXlast  = X + N*nX;
  
  if(*trace){
    Rprintf("\n--------------------------------------------------\n");
    Rprintf("(nY, nZ, nX, N, maxit) = (%d, %d, %d, %d, %d)\n", nY, nZ, nX, N, maxit);
    Rprintf("(useOffset, min_nT, min_N, min_Nhet) = (%d, %d, %d, %d)", 
            useOffset, min_nT, min_N, min_Nhet);
    Rprintf("\n--------------------------------------------------\n");
  }
  
  /* time records */
  time_t sec_s;
  time_t sec_e;
  
  /* starting time */
  sec_s = time(NULL);
  
  /* **********************************************************
   * identifify eQTL gene by gene
   * **********************************************************/
  
  /* pY/pY1/pY2 are pointers to gene expression data */
  pY  = Y;
  pY1 = Y1;
  pY2 = Y2;
  
  for(i=0; i<nY; i++,pY+=N,pY1+=N,pY2+=N){
    
    for (k=0; k<4; k++) {
      pv_min[k] = 1.0;
      bestM[k]  = -9;
    }
    
    useASE  = 1;
    useTReC = 1;
    adjZ    = 1;
    
    if(*trace){
      Rprintf("\ni=%d\n", i);
    }
    
    /* **********************************************************
     * fit a baseline model using only the confouding covariates 
     * family is assigned to a value of either Poisson or NB
     * **********************************************************/
    
    dimsNew[1] = nX;
    dimsNew[3] = 0; /* no initial values. NOTE, dimsNew[3] will be 
                       changed within glmNB, dumb. I know */
    
    convBase = glmNB(dimsNew, &nIter, pY, z1, &linkR, offset, X, convGLM, 
                     &rank0, Xb, fitted0, resid, weights, &phi0, &scale, 
                     &df0, &family, &twoLL_trec0, scoreTestP, trace, &btmp);
    
    if (!convBase){
      useTReC = 0;
      if (*trace)
        Rprintf("  i=%d, fail to fit baseline TReC model\n", i);
    }
    
    if(*trace > 1){
      Rprintf("\n  finish fitting baseline TreC model twologlik0=%.4e, ", twoLL_trec0);
      Rprintf("family=%d, phi0=%.4e\n", family, phi0);
    }
    
    /* ======================================================== */
    exPara[2] = (double)family;
    /* ======================================================== */
    
    /* *****************************************************
     * calculate nA and nTotal
     * *****************************************************/
    
    h0 = 0;
    
    for (k=0; k<N; k++) {
      nT1 = pY1[k] + pY2[k];

      if(nT1 < min_nT){ continue; }
      
      nTotal[h0] = nT1;
      nA[h0]     = pY1[k];
      
      h0++;
    }
    
    /* if sample size is not enough */
    if(h0 < min_N) useASE = 0;
    
    if(useASE){
      /* *****************************************************
       * obtain MLE for H0: no allelic imbalance
       * *****************************************************/
      
      exPara[1]   = (double) h0;
      exPara[6]   = 0.5; /* fixed pi */
      
      initPara[0] = 0.1; 
      initPara[1] = 0.5; /* value for pi, not used for H0 */
      npara       = 1;
      
      // Rprintf("trecase_max1: i=%d, h0=%d, nbd[0]=%d, range=(%.2e,%.2e), starting lbfgsb\n", 
      //        i, h0, nbd[0], lower[0], upper[0]);

      //lbfgsb(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
      //       negLogH0, negGradLogH0, &fail, (void*)exPara, factr, pgtol,  
      //       &fncount, &grcount, maxit, msg, 0, nREPORT);
      lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
             negLogH0, negGradLogH0, &fail, (void*)exPara, factr, pgtol,  
             &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,x1);      
      if (fail) {
        if (*trace)
          Rprintf("  i=%d, fail to fit baseline ASE model\n", i);
        
        useASE = 0;
        th0    = -1000.00;
        twoLL_ase0 = 1e16;
      }else{
        th0 = initPara[0];
        twoLL_ase0 = -2.0*Fmin;
        
        if(*trace > 1){
          negGradLogH0(npara, initPara, gr, (void*)exPara,x1);
          Rprintf("\n  Obtained MLE for H0: twologlik0=%.4e, ", twoLL_ase0);
          Rprintf("theta0=%.4e, gradience=%.4e\n", th0, *gr);
        }
      }
    }
    
    yFailBaselineModel[i] = 1 - useTReC + 2*(1 - useASE);
    
    if (*trace > 1)
      Rprintf("  i=%d, h0=%d, useTReC=%d, useASE=%d, phi0=%e, theta0=%e\n", 
              i, h0, useTReC, useASE, phi0, th0);
    
    
    /* **********************************************************
     * Now start to fit models with both confounding covariates
     * and each of the SNPs 
     * **********************************************************/
    
    pZ  = Z;
    pZh = Zh;
    
    for(j=0; j<nZ; j++,pZ+=N,pZh+=N){
      
      adjZ = 1;

      if(*cis_only){
        if(eChr[i] != mChr[j]) continue;
        
        pos_diff = abs(ePos[i] - mPos[j]);
        
        if(pos_diff > *cis_distance) continue;
      }
      
      if(*trace > 2){
        Rprintf("\ni=%d, j=%d\n", i, j);
      }
      
      useTReC_j = useTReC;
      useASE_j  = useASE;
      
      h1 = 0;
      
      /* **********************************************************
       * Initial fitting with TReC and ASE
       * **********************************************************/
      
      if (useTReC){
        /* -------------------------------------------------------
         * now fit a glmNBlog model, with adjustment of Z
         * ----------------------------------------------------- */
        
        for (k=0; k<N; k++) {
          pXlast[k]  = pZ[k];
          fitted1[k] = fitted0[k];
        }
        
        phi1  = phi0;
        scale = 1.0;
        
        dimsNew[1] = nX;
        dimsNew[3] = 1; /* use initial values */
        
        convSNPj = glmNBlog(dimsNew, &nIter, pY, z1, &linkR, offset, X, conv, convGLM, 
                            &rank1, Xb, fitted1, resid, weights, &phi1, &scale, 
                            &df1, &family, &twoLL_trec1, scoreTestP, trace, 
                            &bxj_trec, fitted2, offsetN);
        
        if(convSNPj == 0){
          for (k=0; k<N; k++) {
            fitted1[k] = fitted0[k];
          }
          
          adjZ = 0;
          dimsNew[1] = nX + 1;
          dimsNew[3] = 1; /* use initial values */
                    
          convSNPj = glmNB(dimsNew, &nIter, pY, z1, &linkR, offset, X, convGLM, 
                           &rank1, Xb, fitted1, resid, weights, &phi1, &scale, 
                           &df1, &family, &twoLL_trec1, scoreTestP, trace, 
                           &bxj_trec);
        }          
        
        if(convSNPj == 0){
          if(*trace)
            Rprintf("\n  Fail TReC: i=%d, j=%d, family=%d\n", i, j, family);
          
          useTReC_j = 0;
        }else if(df0 - df1 != 1){
          if(*trace){
            Rprintf("\n  i=%d, j=%d, dfs=(%d, %d), ranks=(%d, %d)\n", 
                  i, j, df0, df1, rank0, rank1);
          }
          if(df0 - df1 < 0.5) useTReC_j = 0;
        }
        
        if (useTReC_j && *trace > 1) {
          Rprintf("\n  finish TReC: i=%d, j=%d, phi=%e, b=%e, ", i, j, phi1, bxj_trec);
          Rprintf("twologlike=%e\n\n", twoLL_trec1);
        }
        
      }
      
      if (useASE) {
        /* *****************************************************
         * calculate nA and nTotal
         * *****************************************************/
        h0 = 0;
        h1 = 0;
        
        for (k=0; k<N; k++) {
          nT1 = pY1[k] + pY2[k];
          
          if(nT1 < min_nT) continue;
          
          if(nTotal[h0] != nT1) error("mismatch ;( \n");
          
          /* pZ[k] = 0 or 4 if homozygous */
          if (fabs(pZh[k] - 2.0) > 1.99) {
            zeta[h0] = -1.0;
            nA[h0]   = pY2[k];
          }else {
            zeta[h0] = 1.0;
            
            if(fabs(pZh[k] - 1.0) < 0.01){
              nA[h0] = pY2[k];
            }else if(fabs(pZh[k] - 3.0) < 0.01){
              nA[h0] = pY1[k];
            }else {
              error("invalid values for Z\n");
            }
            
            h1++;
          }
          h0 ++;
        }
        
        /* if sample size of heterzygous genotype is not enough */
        if(h1 < min_Nhet){
          useASE_j=0;
        }else{
          /* *****************************************************
           * obtain MLE for H1: with allelic imbalance
           * *****************************************************/
          
          exPara[1]   = (double) h0;
          initPara[0] = th0; 
          initPara[1] = 0.5;
          npara       = 2;
          
          // Rprintf("trecase_max1: i=%d, j=%d, h0=%d, th0=%e, nbd[0]=(%d, %d), starting lbfgsb\n", 
          //        i, j, h0, th0, nbd[0], nbd[1]);
                    
          //lbfgsb(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
          //       negLogH1, negGradLogH1, &fail, (void*)exPara, factr, pgtol,  
          //       &fncount, &grcount, maxit1, msg, 0, nREPORT);
          lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
                 negLogH1, negGradLogH1, &fail, (void*)exPara, factr, pgtol,  
                 &fncount, &grcount, maxit1, msg, 0, nREPORT, wa, iwa, g1,x1);
          
          if (fail){
            useASE_j=0;
            
            if (*trace)
              Rprintf("  i=%d, j=%d, fail ASE model\n", i, j);
            
          }else{
            twoLL_ase1 = -2.0*Fmin;
            
            theta_ase = initPara[0];
            pi_ase    = initPara[1];
            bxj_ase   = log(pi_ase) - log(1.0-pi_ase);
          }
        }
        
        if (*trace > 2){
          Rprintf("  i=%d, j=%d, h0=%d, h1=%d, useASE_j=%d\n", 
                  i, j, h0, h1, useASE_j);
          
          if(useASE_j)
            Rprintf("  pi_ase=%e, theta_ase=%e, bxj_ase=%e\n", 
                    pi_ase, theta_ase, bxj_ase);
          
        }
      }
      
      /* **********************************************************
       * Iterative updating
       * **********************************************************/
      
      useJointModel = 0;
      
      if(useTReC_j && useASE_j && adjZ){
        
        useJointModel = 1;
        
        bxj_old   = 0.0;
        theta_old = th0;
        phi_old   = phi0;
        pi_old    = 0.5;
        
        if (*trace > 2)
          Rprintf("\n  i=%d, j=%d, use joint models\n", i, j);
        
        for (k=0; k<N; k++) {
          fitted2[k] = fitted0[k];
        }
        
        /**********************************
         * calculate old likelihood       
         **********************************/
        
        if (family==NB) {
          twoLL_trec_joint0 = 2.0*loglik_NB(N, phi_old, fitted2, pY);
        }else if (family==POISSON) {
          twoLL_trec_joint0 = 2.0*loglik_Poisson(N, fitted2, pY);
        }else {
          error("invalid family\n");
        }
        
        para_ase[0] = theta_old;
        para_ase[1] = pi_old;
        
        twoLL_ase_joint0 = - 2.0*negLogH1(2, para_ase, exPara,x1);
        
        twoLL0 = twoLL_trec_joint0 + twoLL_ase_joint0;
        twoLL1 = twoLL0;
        
        /* --------------------------------------------------------
         * iterations to estimate bxj
         * -------------------------------------------------------*/
        
        for(g=0; g<maxit; g++){
          
          paraDiff = 0.0;
          
          /* --------------------------------------------------------
           * 1. Given theta, phi, b_0, b_k, and b_u, estimate b_{x_j}
           * -------------------------------------------------------*/
          
          /*******************************************************
           * estimate b_xj      
           *******************************************************/
          twoLL_bxj0  = twoLL1;
          
          fail = 0;
          
          b_ml(&bxj_new, N, h0, bxj_old, phi_old, theta_old, pY, 
               pZ, fitted2, nA, nTotal, zeta, maxit, *convGLM, 
               *trace, &fail);
          
          if(fail){
            if(*trace){
              Rprintf("\n  i=%d, j=%d, fail to estimate bxj in b_ml (trecase_max1)\n", i, j);
            }
            
            useJointModel = 0;
            break;
          }
                    
          if (paraDiff < fabs(bxj_new - bxj_old)) 
            paraDiff = fabs(bxj_new - bxj_old);
          
          /**********************************
           * calculate new likelihood       
           **********************************/
          pi_old = exp(bxj_new)/(1.0 + exp(bxj_new));
          
          for (k=0; k<N; k++){
            pZk = pZ[k];
            
            if (fabs(pZk - 1.0) < 0.01) {
              fitted2[k] = exp(log(fitted2[k]) - log(1.0 + exp(bxj_old)) + log(1.0 + exp(bxj_new)));
            }else if (fabs(pZk - 2.0) < 0.01) {
              fitted2[k] = exp(log(fitted2[k]) + (bxj_new - bxj_old));
            }
          }
          
          if (family==NB) {
            twoLL_trec_joint1 = 2.0*loglik_NB(N, phi_old, fitted2, pY);
          }else if (family==POISSON) {
            twoLL_trec_joint1 = 2.0*loglik_Poisson(N, fitted2, pY);
          }else {
            error("invalid family\n");
          }
          
          para_ase[0] = theta_old;
          para_ase[1] = pi_old;
          
          twoLL_ase_joint1 = - 2.0*negLogH1(2, para_ase, exPara,x1);
          
          twoLL_bxj1 = twoLL_trec_joint1 + twoLL_ase_joint1;
          
          if (*trace > 1) {
            Rprintf("\n  i=%d, j=%d, g=%d\n", i, j, g);
            Rprintf("  theta=%e, phi=%e, b_old=%e, b_new=%e\n", 
                    theta_old, phi_old, bxj_old, bxj_new);
            
            Rprintf("  twoLL(old, new)=(%e, %e)\n", twoLL_bxj0, twoLL_bxj1);
          }
          
          if((twoLL_bxj1 - twoLL_bxj0)/fabs(twoLL_bxj0) < -0.01){
            Rprintf("\n  i=%d, j=%d, g=%d\n", i, j, g);
            error("  likelihood decreases during update of b_xj\n");
          }
          
          bxj_old = bxj_new;
          
          /* --------------------------------------------------------
           * 2. Given b_{x_j}, estimate theta
           * -------------------------------------------------------*/
          twoLL_ase_joint0 = twoLL_ase_joint1;
          
          exPara[1]   = (double) h0;
          exPara[6]   = pi_old; /* fixed pi */
          initPara[0] = theta_old; 
          initPara[1] = pi_old; /* value for pi, not used for H0 */
          npara       = 1;
          
          // Rprintf("trecase_max1: i=%d, j=%d, g=%d, h0=%d, pi=%e, th0=%e, nbd[0]=(%d, %d), starting lbfgsb\n", 
          //        i, j, g, h0, pi_old, theta_old, nbd[0], nbd[1]);

          //lbfgsb(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
          //       negLogH0, negGradLogH0, &fail, (void*)exPara, factr, pgtol,  
          //       &fncount, &grcount, maxit, msg, 0, nREPORT);
          lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
                 negLogH0, negGradLogH0, &fail, (void*)exPara, factr, pgtol,  
                 &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1,x1);          

          twoLL_ase_joint1 = -2.0*Fmin;
          
          if (fail) {
            if(*trace){
              Rprintf("\n  i=%d, j=%d, h0=%d, fail estimating theta in joint model\n", 
                      i, j, h0);
              negGradLogH0(npara, initPara, gr, (void*)exPara,x1);
              Rprintf("  theta=%.4e, gradience=%.4e, fail=%d\n", initPara[0], gr[0], fail);
            }
            useJointModel = 0;
            break;
          }
          
          theta_new = initPara[0];
          
          if (paraDiff < fabs(theta_new - theta_old)) 
            paraDiff = fabs(theta_new - theta_old);
          
          /**********************************
           * compare old and new likelihood       
           **********************************/
          
          if((twoLL_ase_joint1 - twoLL_ase_joint0)/fabs(twoLL_ase_joint0) < -0.01 ){
            Rprintf("\n  i=%d, j=%d, g=%d ", i, j, g);
            Rprintf("  theta_old=%e, theta_new=%e\n", theta_old, theta_new);
            Rprintf("  twoLL(old, new)=(%e, %e)", twoLL_ase_joint0, twoLL_ase_joint1);
            error("\n  likelihood decreases during update of theta in join modeling\n");
          }
          
          theta_old = theta_new;
          
          /* --------------------------------------------------------
           * 3. Given b_{x_j} and phi, estimate b_0, b_k, b_u, phi
           * in fact, we do not estimate b_0, b_k, and b_u, 
           * instead we only need to estimate mu, or fitted1.
           * -------------------------------------------------------*/
          twoLL_trec_joint0 = twoLL_trec_joint1;
          
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
          
          phi_new = phi_old;
          convSNPj = glmNB(dimsNew, &nIter, pY, z1, &linkR, offsetN, X, convGLM, 
                           &rank1, Xb, fitted2, resid, weights, &phi_new, &scale, 
                           &df1_join, &family, &twoLL_trec_joint1, scoreTestP, trace, 
                           &btmp);
          
          if (convSNPj == 0) {
            if(*trace)
              Rprintf("\n  i=%d, j=%d, fail to estimate phi in joint model\n", i, j);
            
            useJointModel = 0;
            break;
          }
          
          /* need to substract the extra degree for freedom in the offset */
          df1_join -= 1;
          
          if (paraDiff < fabs(phi_new - phi_old)) paraDiff = fabs(phi_new - phi_old);
          
          /**********************************
           * compare likelihood       
           **********************************/
          
          if((twoLL_trec_joint1 - twoLL_trec_joint0)/fabs(twoLL_trec_joint0) < -0.01){
            Rprintf("\n  i=%d, j=%d, g=%d ", i, j, g);
            Rprintf("  twoLL(old, new)=(%e, %e)", twoLL_trec_joint0, twoLL_trec_joint1);
            error("\n  likelihood decreases during update of phi in join modeling\n");
          }
          
          phi_old = phi_new;
          twoLL_trec_joint0 = twoLL_trec_joint1;
          
          if (*trace > 1){
            Rprintf("\n  i=%d, j=%d, g=%d, parDiff=%e", i, j, g, paraDiff);
            Rprintf("\n  bxj=%e, theta=%e, phi=%e\n", bxj_old, theta_old, phi_old);
          }
          
          /* --------------------------------------------------------
           * check convergence
           * -------------------------------------------------------*/
          
          if(paraDiff < *conv){
            twoLL1  = twoLL_trec_joint1 + twoLL_ase_joint1;
            
            if (*trace > 1){
              Rprintf("\n  i=%d, j=%d, converged using joint model\n", i, j);
              Rprintf("  bxj=%e, theta=%e, phi=%e, twoLL0=%e, twoLL1=%e\n", 
                      bxj_old, theta_old, phi_old, twoLL0, twoLL1);
            }
            
            break;
          }
        }
        
        if (g >= maxit) {
          if (*trace > 1){
            Rprintf("\n  i=%d, j=%d, reach max iteration using joint model. ", i, j);
            Rprintf("g=%d, maxit=%d, paraDiff=%.3e\n", g, maxit, paraDiff);
          }
          useJointModel = 0;
        }
        
      }
      
      /* *********************************************************
       * calculate p-value for TReC model
       * *********************************************************/
      
      if (useTReC_j) {
        chisqTReC = twoLL_trec1 - twoLL_trec0;
        
        
        if (chisqTReC < -0.01) {
          nless5 = 0.0;
          for (k=0; k<N; k++) { 
            if (pY[k] < 5) nless5 += 1.0;
          }
          
          Rprintf("  g=%d, nless5=%f, N=%d, 2logL(old, new)=(%e, %e)\n", 
                  g, nless5, N, twoLL_trec0, twoLL_trec1);
          
          useTReC_j=0;

          if (nless5 <= 0.75*N) {
            error("likelihood decreases for TReC model: i=%d, j=%d, 2logL=(%.2e, %.2e), chisq=%.3e\n", 
                  i, j, twoLL_trec0, twoLL_trec1, chisqTReC);
          }
          
        }else{
        
          dfr_TReC = (double)(df0 - df1);
          
          if(dfr_TReC < 0.5){
            pvalTReC = 995.0;
            if (*trace > 1) {
              Rprintf("dfr_TReC=%e, df0=%d, df1=%d\n", dfr_TReC, df0, df1);
            }
          }else{
            
            if (chisqTReC < 1e-5) { 
              pvalTReC = 1.0; 
            }else{
              pvalTReC  = pchisq(chisqTReC, dfr_TReC, 0, 0);
            }
            
          }
        }
      }
      
      if(! useTReC_j){
        chisqTReC = -995.0;
        dfr_TReC  = -1.0;
        pvalTReC  =  995.0;
      }
      
      /* *********************************************************
       * calculate p-value for ASE model
       * *********************************************************/
      
      if(useASE_j) {
        chisqASE = twoLL_ase1 - twoLL_ase0;
        
        if (chisqASE < -0.01) {
          error("likelihood decreases for ASE model: i=%d, j=%d, 2logL=(%.2e, %.2e), chisq=%.3e\n", 
                i, j, twoLL_ase0, twoLL_ase1, chisqASE);
        }
        
        dfr_ASE = 1.0;
        if (fabs(th0) < 1e-7)  dfr_ASE += 1.0;
        if (fabs(theta_ase) < 1e-7)  dfr_ASE -= 1.0;
        
        if(dfr_ASE < 0.5){
          pvalASE = 995.0;
          if (*trace > 1) {
            Rprintf("dfr_ASE=%e, df0=%d, df1=%d\n", dfr_ASE, df0, df1);
          }          
        }else{
          if (chisqASE < 1e-5) { 
            pvalASE = 1.0; 
          }else{
            pvalASE = pchisq(chisqASE, dfr_ASE, 0, 0);
          }
        }
        
      }else {
        chisqASE = -995.0;
        dfr_ASE  = -1.0;
        pvalASE  =  995.0;
      }
      
      /* *********************************************************
       * calculate p-value for joint model
       * *********************************************************/
      
      if(useJointModel){
        
        chisqTrans = twoLL_ase1 + twoLL_trec1 - twoLL1;
        
        if (chisqTrans < -0.01) {
          Rprintf("twoLL(ASE, TReC, TReCASE)=(%.3e, %.3e)\n",  twoLL_ase1, twoLL_ase1, twoLL1);
          warning("testing for cis v.s. trans, i=%d, j=%d, chisq=%.3e\n", i, j, chisqTrans);
        }
        
        if (chisqTrans < 1e-5) { 
          pvalTrans = 1.0; 
        }else{
          pvalTrans = pchisq(chisqTrans, 1.0, 0, 0);
        }
        
        chisqJoint = twoLL1 - twoLL0;
        
        if (chisqJoint < -0.01) {
          Rprintf("twoLL_trec=(%.3e, %.3e) ",  twoLL_trec0, twoLL_trec_joint1);
          Rprintf("twoLL_ase=(%.3e, %.3e) \n", twoLL_ase0,  twoLL_ase_joint1);
          
          error("wrong twoLL for joint model i=%d, j=%d, chisq=%.3e\n", 
                i, j, chisqJoint);
        }
        
        dfr_joint_ASE = 1.0;
        if (fabs(th0) < 1e-7)  dfr_joint_ASE += 1.0;
        if (fabs(theta_new) < 1e-7)  dfr_joint_ASE -= 1.0;
        
        dfr_joint_TReC = (double)(df0 - df1_join);
        
        /* we need to subtract 1 here because in the joint model, 
         * b (of TReC model) and pi (of ASE model) are the same */
        dfr_Joint = dfr_joint_ASE + dfr_joint_TReC - 1.0;
        
        if(dfr_Joint < 0.5){
          pvalJoint = 995.0;
          if (*trace > 1) {
            Rprintf("dfr_Joint=%e, dfr_joint_ASE=%e, df0=%d, df1_join=%d\n", 
                    dfr_Joint, dfr_joint_ASE, df0, df1);
          }          
        }else{
          if (chisqJoint < 1e-5) { 
            pvalJoint = 1.0; 
          }else{
            pvalJoint = pchisq(chisqJoint, dfr_Joint, 0, 0);
          }
        }
        
      }else {
        chisqJoint = -995.0;
        dfr_Joint  = -1.0;
        pvalJoint  =  995.0;
        
        chisqTrans = -995.0;
        pvalTrans  =  995.0;
      }
      
      /* minimum p-value for TReC model */
      if(pvalTReC < pv_min[0]){
        pv_min[0] = pvalTReC;
        bestM[0]  = j; /* here j is markerID - 1*/
      }

      /* minimum p-value for ASE model */
      if(pvalASE < pv_min[1]){
        pv_min[1] = pvalASE;
        bestM[1]  = j; /* here j is markerID - 1*/
      }
      
      /* minimum p-value for Joint model */
      if(pvalJoint < pv_min[2]){
        pv_min[2] = pvalJoint;
        bestM[2]  = j; /* here j is markerID - 1*/
      }
            
      /* p-value of TReC or TReCASE model, depends on  
       * it is cis-eQTL or not */
      pv1 = 1.0;      
      if(pvalTrans < *transTestP || pvalTrans > 1.0001){
        pv1 = pvalTReC;
      }else {
        pv1 = pvalJoint;
      }
      
      if(pv1 < pv_min[3]){
        pv_min[3] = pv1;
        bestM[3]  = j; /* here j is markerID - 1 */
      }
      
    }
    
    for (k=0; k<4; k++) {
      pval[i + k*nY]   = pv_min[k];
      best_m[i + k*nY] = bestM[k];
    }
  }
  /* end time */
  sec_e  = time(NULL);
  if(*trace){
    Rprintf("\n------------------------------------------------------\n");
    Rprintf("total time spent is %ld secs\n", sec_e-sec_s);
    Rprintf("------------------------------------------------------\n");
  }
  
  *succeed = 2;
}

/**********************************************************************
 *
 * trecase_permute
 *
 * map eQTL by TReC + ASE model for cis- only eQTL mapping. 
 *
 * Carry out fixed number of permutations, and calculate permutation
 * p-value of the best association of each gene expression trait
 *
 **********************************************************************/
//may improve a bit, by moving here and passing to trecase_max1 , double *wa, int *iwa, double *g, SEXP x1

void trecase_permute (int* dims, double* Y, double* X, double* Z,  
                      double* z1, double* Y1, double* Y2, double* Zh,  
                      int* npermute, double* offset, int* cis_only, 
                      int* cis_distance, int* eChr, int* ePos, 
                      int* mChr, int* mPos, double* conv, 
                      double* convGLM, int* yFailBaselineModel, 
                      double* scoreTestP, double* transTestP, 
                      int* best_m, double* pval, double* perPval, 
                      int* trace, int* succeed)
{
  int i, j, p, q, gap;
  time_t timer;
  
  int nY = dims[0];
  int nX = dims[1];
  int nZ = dims[2];
  int N  = dims[3];
  
  int nPer = *npermute;
  
  if(*trace){
    Rprintf("nY=%d, nX=%d, nZ=%d, N=%d, nPer=%d\n", nY, nX, nZ, N, nPer);
  }
  
  int *best_m0, *perm1;
  double *pval0, *perY, *perY1, *perY2, *perX, *perm2;
  double tmp=0.0, *ppY1, *ppY2;
  
  /* Work array */
  double *Xb, *fitted0, *fitted1, *fitted2, *resid;
  double *weights, *offsetN, *exPara;
  
  int npara,lmm;
  npara   = NPARA;
  lmm     = LMM;
  //technical parameters below:
  double *wa, *g1;
  int *iwa;
  SEXP x1;
  PROTECT(x1 = allocVector(REALSXP, npara));
  //consider replacing with simple Calloc
  wa  = (double *) S_alloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm,sizeof(double));
  iwa = (int*) R_alloc(3*npara,sizeof(int));
  g1 = (double *)R_alloc(npara, sizeof(double));
  //end technical parameters

  /***
   * best (associated) marker and p-value for each gene 
   * in each permutation
   */  
  best_m0 = (int *) Calloc(4*nY, int);
  pval0   = (double *) Calloc(4*nY, double);
  
  /* pernutation index */
  perm1   = (int *)Calloc(N, int);
  perm2   = (double *) Calloc(N, double);

  /* permuted expression values */
  perY    = (double *) Calloc(nY*N, double);
  perY1   = (double *) Calloc(nY*N, double);
  perY2   = (double *) Calloc(nY*N, double);
  perX    = (double *) Calloc((nX+1)*N, double);

  /* initial permutation p-value */
  for(q=0; q<4*nY; q++){ 
    perPval[q] = 0.0; 
    pval[q]    = 1.0;
    best_m[q]  = -9;
  }
  
  /* allocate memory. The extra column in Xb is used to store 
   * one SNP, the same as for X
   */
  
  Xb      = (double *) Calloc(N*(nX+1), double);
  fitted0 = (double *) Calloc(N, double); /* for Null TReC model */
  fitted1 = (double *) Calloc(N, double); /* for Alternative TReC model */
  fitted2 = (double *) Calloc(N, double); /* for Alternative TReC model */
  resid   = (double *) Calloc(N, double);
  weights = (double *) Calloc(N, double);
  offsetN = (double *) Calloc(N, double);
  
  exPara = (double *) Calloc(3*N+7, double);

  /***
   * find the smallest p-value for each gene with unpermuted data 
   */
  timer=time(NULL);
  
  Rprintf("\nun-permuted data %s\n", asctime(localtime(&timer)));
  
  *succeed = 0;
  
  trecase_max1 (dims, Y, X, Z, z1, Y1, Y2, Zh, offset, 
                cis_only, cis_distance, eChr, ePos, mChr, mPos, 
                conv, convGLM, yFailBaselineModel, scoreTestP, 
                transTestP, best_m, pval, trace, succeed, 
                Xb, fitted0, fitted1, fitted2, resid, 
                weights, offsetN, exPara,
				wa,iwa,g1,x1);

  if(*succeed != 2){
    error("fail to finish trecase_max1 in trecase permute\n");
  }
  
  /***
   * carry out permutations 
   */
  
  GetRNGstate();
  
  /**
   * gap to control the print out tracing information.
   */
  gap = ceil(500.0/nY);
  gap = 10*gap;
  
  for(p=1; p<=nPer; p++){
    
    if (*trace) {
      Rprintf("p=%d\n", p);
    }
    
    if(p % gap == 0){
      timer=time(NULL);
      Rprintf("\n%d-th permutation %s\n", p, asctime(localtime(&timer)));
    }
        
    getPermute(perm1, N);
    if(*trace){
      Rprint_vi(perm1, 0, N-1);
    }
    
    permute(Y, perY, perm1, nY, N);
    permute(X, perX, perm1, nX, N);
    permute(Y1, perY1, perm1, nY, N);
    permute(Y2, perY2, perm1, nY, N);

    /**
     * permute Y1 and Y2
     */
    for(j=0; j<N; j++){ 
      perm2[j] = unif_rand(); 
    }
    
    if(*trace){
      Rprint_v(perm2, 0, N-1);
    }
    
    ppY1 = perY1;
    ppY2 = perY2;

    for(i=0; i<nY; i++,ppY1+=N,ppY2+=N){
      for(j=0; j<N; j++){

        if(perm2[j] > 0.5){
          tmp     = ppY1[j];
          ppY1[j] = ppY2[j];
          ppY2[j] = tmp;
        }
        
      }
    }
    
    for(q=0; q<4*nY; q++){ 
      pval0[q]    = 1.0;
      best_m0[q]  = -9;
    }    
    
    *succeed = 0;
    
    trecase_max1 (dims, perY, perX, Z, z1, perY1, perY2, Zh, offset, 
                  cis_only, cis_distance, eChr, ePos, mChr, mPos, 
                  conv, convGLM, yFailBaselineModel, scoreTestP, 
                  transTestP, best_m0, pval0, trace, succeed,
                  Xb, fitted0, fitted1, fitted2, resid, 
                  weights, offsetN, exPara,
				  wa,iwa,g1,x1);
    
    if(*succeed != 2){
      error("fail to finish trecase_max1 in trecase permute\n");
    }
    
    for(q=0; q<4*nY; q++){
      if(pval0[q] <= pval[q]) perPval[q] += 1.0; 
    }
    
  }
  
  PutRNGstate();
  
  for(q=0; q<4*nY; q++){ perPval[q] /= nPer; }
  
  *succeed = 0;

  Free(best_m0);
  Free(pval0);
  Free(perm1);
  Free(perm2);
  Free(perY);
  Free(perY1);
  Free(perY2);
  Free(perX);

  Free(Xb);
  Free(fitted0);
  Free(fitted1);
  Free(fitted2);
  Free(resid);
  Free(weights);
  Free(offsetN);
  Free(exPara);
  UNPROTECT(1);
    
  *succeed = 1;
}
