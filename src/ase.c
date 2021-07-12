/*
 *  ase.c
 *
 *  Created by Wei Sun on 5/25/2010.
 *  Modified by Vasyl Zhabotynsky on 04/05/2011 
 *
 */

#include "ase.h"

/**********************************************************************
 *
 * negative log likelihood and gradient function under H0
 *
 * H0: pi = 0.5
 *
 * The first parameter (n) is the number of paramters, n=1
 * Sample size N is included in the parameter "ex".
 *
 **********************************************************************/

double negLogH0 (int n, double* para, void* ex, SEXP x1){
  int i, k, N, h;
  double sumL, ni, ni0, pi0, piI, theta;
  double *exPara, *nA, *nTotal, *zeta;
  
  exPara = (double *) ex;
  N      = ceil(exPara[0]-0.5);
  h      = ceil(exPara[1]-0.5);
  pi0    = exPara[6];
  nA     = exPara + 7;
  nTotal = nA + N;
  zeta   = nTotal + N;
  
  theta  = para[0];
  
  sumL = 0.0;
  
  for(i=0; i<h; i++){
    ni   = nTotal[i];
    ni0  = nA[i];
    
    sumL += lchoose(ni, ni0);
    
    if(zeta[i] > 0){ piI = pi0; }else { piI = 0.5; }
    
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
  
  return(-sumL);
}


/**********************************************************************
 *
 * negGradLogH0A: gradient negLogH0
 *
 **********************************************************************/

void negGradLogH0(int n, double* para, double* gr, void* ex, SEXP x1)
{
  double grad, ni, ni0, pi0, piI, theta;
  int i, k, N, h;
  double *exPara, *nA, *nTotal, *zeta;
  
  exPara = (double *) ex;
  N = ceil(exPara[0]-0.5);
  h = ceil(exPara[1]-0.5);
  pi0    = exPara[6];
  nA     = exPara + 7;
  nTotal = nA + N;
  zeta   = nTotal + N;
  
  theta = para[0];

  grad = 0.0;
  
  for (i=0; i<h; i++) {
    ni  = nTotal[i];
    ni0 = nA[i];
    
    if(zeta[i] > 0){ piI = pi0; }else { piI = 0.5; }

    if(ni0 > 0){
      for(k=1; k<ni0; k++) grad += k/(piI + k*theta);
    }
    
    if(ni0 < ni){
      for(k=1; k<ni-ni0; k++) grad += k/(1.0 - piI + k*theta);
    }

    for(k=1; k<ni; k++) grad -= k/(1.0 + k*theta);
  }
  
  gr[0] = -grad;
}

/**********************************************************************
 *
 * negative log likelihood and gradient function under H1
 *
 * H1: with alleleic imbalance, and pi_0 = 0.5
 *
 * Note the first parameter (n) is the number of paramters, n=2
 * Sample size N is included in the parameter "ex".
 *
 **********************************************************************/

double negLogH1 (int n, double* para, void* ex, SEXP x1){
  int i, k, N, h;
  double sumL, ni, ni0;
  double pi1, piI, theta;
  double *exPara, *nA, *nTotal, *zeta;
  
  exPara = (double *) ex;
  N      = ceil(exPara[0]-0.5);
  h      = ceil(exPara[1]-0.5);

  nA     = exPara + 7;
  nTotal = nA + N;
  zeta   = nTotal + N;
  
  theta = para[0];
  pi1   = para[1];

  sumL = 0.0;
  
  for(i=0; i<h; i++){
    ni   = nTotal[i];
    ni0  = nA[i];
    
    sumL += lchoose(ni, ni0);
    
    if(zeta[i] > 0){ piI = pi1; }else { piI = 0.5; }
  
    if(ni0 > 0){
      for(k=0; k<ni0; k++) sumL += log(piI + k*theta);
    }
    
    if(ni0 < ni){
      for(k=0; k<ni-ni0; k++) sumL += log(1.0 - piI + k*theta);
    }
    
    for(k=0; k<ni; k++){
      sumL -= log(1.0 + k*theta);
    }
    
  }
  
  return(-sumL);
}

/**********************************************************************
 *
 * negGradLogH1: gradient negLogH1
 *
 **********************************************************************/

void negGradLogH1 (int n, double* para, double* gr, void* ex, SEXP x1)
{
  int i, k, N, h;
  double gradPi1, gradTh, pi1, piI, theta, tmp;
  double *exPara, *nA, *nTotal, *zeta, ni, ni0;
  
  exPara = (double *) ex;
  N      = ceil(exPara[0]-0.5);
  h      = ceil(exPara[1]-0.5);

  nA     = exPara + 7;
  nTotal = nA + N;
  zeta   = nTotal + N;

  theta = para[0];
  pi1   = para[1];
  
  gradPi1 = 0.0;
  gradTh  = 0.0;
  
  for(i=0; i<h; i++){
    ni   = nTotal[i];
    ni0  = nA[i];
    
    if(zeta[i] > 0){
      piI = pi1;
      
      if(ni0 > 0){
        for(k=0; k<ni0; k++){
          tmp = 1.0/(piI + k*theta);
          
          gradPi1 += tmp;
          gradTh  += k*tmp;
        }
      }
      
      if(ni0 < ni){
        for(k=0; k<ni-ni0; k++){
          tmp = 1.0/(1.0 - piI + k*theta);
          
          gradPi1 -= tmp;
          gradTh  += k*tmp;
        }
      }
      
    }else {
      piI = 0.5;
      
      if(ni0 > 0){
        for(k=0; k<ni0; k++){
          tmp = 1.0/(piI + k*theta);
          
          gradTh += k*tmp;
        }
      }
      
      if(ni0 < ni){
        for(k=0; k<ni-ni0; k++){
          tmp = 1.0/(1.0 - piI + k*theta);
          
          gradTh += k*tmp;
        }
      }
      
    }
    
    for(k=0; k<ni; k++){
      gradTh -=  k/(1 + k*theta);
    }

  }
  
  gr[0] = -gradTh;
  gr[1] = -gradPi1;
}

/**********************************************************************
 *
 * ase
 *
 * allele specific expression. 
 *
  Input:
  
  Y1, Y2  Expression on two haplotypes
  Z       Covariates of interest
 **********************************************************************/


void ase (int* dims, double* Y1, double* Y2, double* Z, char** output, 
             double* RP_cut, int* cis_only, int* cis_distance, 
             int* eChr, int* ePos, int* mChr, int* mPos,  
             int* trace, int* succeed)
{
  int i, j, k, h0, h1;
  double chisq, pval, loglik0, loglik1, nT1, dfr;
  double pi0, th0, pi1, th1;
  double *exPara, *nA, *nTotal, *zeta;
  int nY = dims[0];
  int nZ = dims[1];
  int N  = dims[2];

  /* 
   * minimum of total number of reads to be considered 
   * If one sample has less than min_nT reads, it is disgarded
   */
  int min_nT = dims[3];
  
  /* 
   * minimum of sample size for testing
   */  
  int min_N  = dims[4];
  
  /* 
   * minimum of sample size for heterzygous genotypes
   */  
  int min_Nhet  = dims[5];
  
  /* 
   * p-value cutoff
   */  
  double P_cut = *RP_cut;
  
  /** 
   * we have to combine nA, nTotal, and zeta into one vector for it usage in 
   * sovling for MLE under H1
   */
  
  exPara = (double *) R_alloc(3*N+7, sizeof(double));
  nA     = exPara + 7;
  nTotal = nA + N;
  zeta   = nTotal + N;
  
  /* initial zeta */
  for (k=0; k<N; k++) {
    zeta[k] = -1.0;
  }
  
  /** 
   * exPara[0]  is sample size, fixed throught the computation 
   * exPara[1]  is sample size used in actual computation 
   *            which will be updated for each possible pair of gene and marker
   * exPara[6]  is pi, which will also be updated in each run
   */
  exPara[0] = (double) N; 
  exPara[1] = 0.0; 
  exPara[6] = 0.0; 

  /* pointers to Y and Z */
  double *pY1, *pY2, *pZ;
    
  /* position difference between gene and marker */
  int pos_diff;
  
  /* grid used to output frequency */
  double grid;
  
  /* 
   * p-value frequencies
   * freqs[100] = #{ pval < P_cut }
   * i = 0:99
   * freqs[i]   = #{ pval < [i/100, (i+1)/100) }
   */
  
  unsigned long freqs[101];
  for(i=0; i<=100; i++){ freqs[i] = 0; }
  
  /* 
   * variables for testing H0
   */
  double lower0, upper0, grH0;
  
  /* **********************************************************
   * parameters for function lbfgsb, which will be used to
     obtain MLE for H1: with allelic imbalance
   
   void lbfgsb(int n, int lmm, double *x, double *lower,
          double *upper, int *nbd, double *Fmin, optimfn fn,
          optimgr gr, int *fail, void *ex, double factr,
          double pgtol, int *fncount, int *grcount,
          int maxit, char *msg, int trace, int nREPORT);
   
   n:       the number of parameters
   
   lmm:     is an integer giving the number of BFGS updates 
            retained in the "L-BFGS-B" method, It defaults to 5.
   
   x:       starting parameters on entry and the final parameters on exit
   
   lower:   lower bounds
   
   upper:   upper bounds
   
   nbd:     specifies which bounds are to be used. 
            nbd(i)=0 if x(i) is unbounded,
            1 if x(i) has only a lower bound,
            2 if x(i) has both lower and upper bounds, and
            3 if x(i) has only an upper bound.
            On exit nbd is unchanged.
   
   Fmin:    final value of the function
   
   fn:      the function to be minimized
   
   gr:      the gradient function
   
   fail:    integer code, 0 for success, 51 for warning and 52 for error
   
   ex:      extra parameters for the function to be minimized
   
   factr:   controls the convergence of the "L-BFGS-B" method. 
            Convergence occurs when the reduction in the objective is 
            within this factor of the machine tolerance. Default is 1e7, 
            that is a tolerance of about 1e-8.
   
   pgtol:   helps control the convergence of the "L-BFGS-B" method. 
            It is a tolerance on the projected gradient in the current 
            search direction. This defaults to zero, when the check 
            is suppressed.
   
   fncount: the number of calls to fn 
   
   grcount: the number of calls to gr
   
   maxit:   maximum of iterations
   
   msg:     A character string giving any additional information 
            returned by the optimizer, or NULL
   
   trace:   Non-negative integer. If positive, tracing information 
            on the progress of the optimization is produced. 
            Higher values may produce more tracing information: 
            for method "L-BFGS-B" there are six levels of tracing. 
   
   nREPORT: The frequency of reports for the "BFGS", "L-BFGS-B" 
            and "SANN" methods if control$trace is positive. 
            Defaults to every 10 iterations for "BFGS" and "L-BFGS-B"
   
   * **********************************************************/
  
  int npara, lmm, fail, failA, failB, fncount, grcount, maxit, nREPORT;
  int nbd[2];

  npara   = 2;
  lmm     = 5;
  fail    = 0;
  failA   = 0;
  failB   = 0;
  fncount = 0;
  grcount = 0;
  maxit   = 100;
  nREPORT = 5;
  nbd[0]  = 1;
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

  double gr[2];
  double initPara[2];
  double lower[2];
  double upper[2];
  double Fmin, factr, pgtol;
  
  /* initPara = c(theta, pi0, pi1) */
  
  initPara[0] = 0.1; 
  initPara[1] = 0.5;

  lower[0] = 0.0 + 1e-16;
  upper[0] = 1e16;
  
  lower[1] = 0.0 + 1e-16;
  upper[1] = 1.0 - 1e-16;

  factr = 1e7;
  pgtol = 0.0;
  
  char msg[1023];
  
  /**********************************************************/

  if(*trace){
    Rprintf("\n--------------------------------------------------\n");
    Rprintf("(nY, nZ, N) = (%d, %d, %d), trace=%d\n", nY, nZ, N, *trace);
    Rprintf("p.cut=%.4e, min.nTotal=%d, min.N=%d, min.Nhet=%d", 
            P_cut, min_nT, min_N, min_Nhet);
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
  fprintf(fo, "GeneRowID\tMarkerRowID\tTheta0\tlogLik0\t");
  fprintf(fo, "Theta1\tPi1\tlogLik1\tChisq\tPvalue\tdf\tN\tNhet\n");
  
  /***
   * identifify eQTL gene by gene
   */
  
  /* pY1/pY2 is the pointer to gene expression data */
  pY1 = Y1;
  pY2 = Y2;
  
  for(i=0; i<nY; i++,pY1+=N,pY2+=N){

    if(*trace > 1){
      Rprintf("\ni=%d\n", i);
    }
    
    /* *****************************************************
     * calculate nA and nTotal
     * *****************************************************/
    h0 = 0;
    
    for (k=0; k<N; k++) {
      nT1 = pY1[k] + pY2[k];
      
      if(nT1 < min_nT){ continue; }
      
      nTotal[h0] = nT1;
      nA[h0]     = pY2[k];
      
      h0++;
    }
    
    /* if sample size is not enough */
    if(h0 < min_N){ continue; }
        
    if(*trace > 1){
      Rprintf("\ni=%d, h0=%d\n", i, h0);
    }
    
    /* *****************************************************
     * obtain MLE for H0
     * situation A: assume pi_0 = 0.5
     * *****************************************************/
    
    exPara[1]   = (double) h0;
    exPara[6]   = 0.5; /* fixed pi0 */
    initPara[0] = 0.1; 
    initPara[1] = 0.5; /* value for pi0, not used for H0A */
    npara       = 1;
    
    //lbfgsb(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
    //       negLogH0, negGradLogH0, &failA, (void*)exPara, factr, pgtol,  
    //       &fncount, &grcount, maxit, msg, 0, nREPORT);
    lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
           negLogH0, negGradLogH0, &failA, (void*)exPara, factr, pgtol,  
           &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1, x1);

    if (failA) {
      if (*trace)
        Rprintf("  i=%d, fail to fit baseline ASE model @ situation A\n", i);
    
      continue;
    }else{
      th0 = initPara[0];
      loglik0 = -Fmin;
      
      if(*trace > 1){
        Rprintf("\n  Obtained MLE for H0: twologlik0=%.4e", loglik0);
        Rprintf(" theta0A=%.4e", th0);
      }
    }
        
    /* *****************************************************
     * start eQTL mapping
     * *****************************************************/
    
    pZ = Z;
    
    for(j=0; j<nZ; j++,pZ+=N){
      
      if(*cis_only){
        if(eChr[i] != mChr[j]) continue;
        
        pos_diff = abs(ePos[i] - mPos[j]);
        
        if(pos_diff > *cis_distance) continue;
      }
            
      /* *****************************************************
       * calculate nA and nTotal
       * *****************************************************/
      h0 = 0;
      h1 = 0;
      
      for (k=0; k<N; k++) {
        nT1 = pY1[k] + pY2[k];
        
        if(nT1 < min_nT){ continue; }
        
        if(nTotal[h0] != nT1){
          error("mismatch ;( \n");
        }
        
        /* pZ[k] = 0 or 4 if homozygous */
        if (fabs(pZ[k] - 2.0) > 1.99) {
          zeta[h0] = -1.0;
          nA[h0]   = pY2[k];
        }else {
          zeta[h0] = 1.0;
          
          if(fabs(pZ[k] - 1.0) < 0.01){
            nA[h0] = pY2[k];
          }else if(fabs(pZ[k] - 3.0) < 0.01){
            nA[h0] = pY1[k];
          }else {
            error("invalid values for Z\n");
          }
          
          h1++;
        }
        h0 ++;
      }
            
      /* if sample size of heterzygous genotype is not enough */
      if(h1 < min_Nhet){ continue; }

      if(*trace > 1){
        Rprintf("\ni=%d, j=%d, h0=%d, h1=%d\n", i, j, h0, h1);
      }
      
      /* *****************************************************
       * obtain MLE for H1A: with allelic imbalance
       * *****************************************************/
      
      exPara[1]   = (double) h0;
      initPara[0] = th0;  /* theta */
      initPara[1] = 0.5;  /* pi1   */
      npara       = 2;
      
      //lbfgsb(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
      //       negLogH1, negGradLogH1, &fail, (void*)exPara, factr, pgtol,  
      //       &fncount, &grcount, maxit, msg, 0, nREPORT);
      lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin, 
             negLogH1, negGradLogH1, &fail, (void*)exPara, factr, pgtol,  
             &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1, x1);
      
      if (fail) {
        if(*trace){
          Rprintf("\n  i=%d, j=%d, h0=%d, h1=%d, fail for MLE of H1, fail=%d\n", 
                  i, j, h0, h1, fail);
        }
        continue;
      }        
      
      th1 = initPara[0];
      pi1 = initPara[1]; 
      loglik1 = -Fmin;
      
      if(*trace > 1){
        Rprintf("\n  Obtained MLE for H1: loglik1=%.4e fail=%d", loglik1, fail);
        Rprintf("\n  theta1=%.4e, pi1=%.4e", initPara[0], initPara[1]);
        negGradLogH1(2, initPara, gr, (void*)exPara, x1);
        Rprintf("\n  grad()=c(%.4e, %.4e)", gr[0], gr[1]);          
      }
            
      chisq = 2.0*(loglik1 - loglik0);

      if (chisq < -1e-5) {
        error("wrong loglik! i=%d, j=%d, loglik=(%.4e, %.4e)\n", 
              i, j, loglik0, loglik1);
      }
            
      dfr = 1.0;
      if (fabs(th0) < 1e-7)  dfr += 1.0;
      if (fabs(th1) < 1e-7)  dfr -= 1.0;

      if(fabs(dfr) < 0.01){
        if(*trace > 1){
          Rprintf("\n  i=%d, j=%d, h0=%d, h1=%d ", i, j, h0, h1);
          Rprintf("dfr=0, theta1=%.4e, pi1=%.4e\n", th1, pi1);
        }
        continue;
      }
      
      if (chisq < 1e-5) { 
        pval = 1.0; 
      }else{
        pval = pchisq(chisq, dfr, 0, 0);
      }
      
      k = (int) (pval / 0.01);
      freqs[k] += 1;
      
      if(pval < P_cut){
        freqs[100] += 1;
        
        /* gene ID and SNP ID */
        fprintf(fo, "%d\t%d\t%e\t%e\t", i+1, j+1, th0, loglik0);
        fprintf(fo, "%e\t%e\t%e\t", th1, pi1, loglik1);
        fprintf(fo, "%.3f\t%.2e\t%.0f\t%d\t%d\n", chisq, pval, dfr, h0, h1);
      }
      
    }
  }
  
  fclose(fo);
  
  // print out the frequencies
  ff   = fopen(output[1], "w");
  grid = 0.0;
  
  fprintf(ff, "<%.2e", P_cut);
  for(i=0; i<100; i++){
    fprintf(ff, "\t%.2f-%.2f", grid, grid+0.01);
    grid += 0.01;
  }
  fprintf(ff, "\n");
  
  fprintf(ff, "%lu", freqs[100]);
  for(i=0; i<100; i++){
    fprintf(ff, "\t%lu", freqs[i]);
  }
  fprintf(ff, "\n");  
  fclose(ff);
  UNPROTECT(1);
  /* end time */
  sec_e  = time(NULL);
  if(*trace){
    Rprintf("\ntotal time spent in ase is %ld secs\n", sec_e-sec_s);
    Rprintf("------------------------------------------------------\n");
  }
    
  *succeed = 1;
}
