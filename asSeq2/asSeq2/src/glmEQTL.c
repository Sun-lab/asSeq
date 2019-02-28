/*
 *  glmEQTL.c
 *
 *  Created by Wei Sun on 5/08/2010.
 *
 */

#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include "glm.h"
#include "utility.h"

/**********************************************************************
 *
 * glmEQTL
 *
 * map eQTL by generalized linear model, Could carry out 
 * cis- only eQTL mapping, resutls are wrttien into two files. 
 *
  Input:
  
  Y   Response variable
  X   Confounding factors
  Z   Covariates of interest
 **********************************************************************/


void glmEQTL (int* dims, double* Y, double* X, double* Z, double* z1, 
              int* link, double* offset, int* adjZ, char** output,  
              double* RP_cut, int* cis_only, int* cis_distance, 
              int* eChr, int* ePos, int* mChr, int* mPos, double* conv, 
              double* convGLM, int* yFailBaselineModel, double* scoreTestP, 
              int* trace, int* succeed)
{
  int i, j, k, c, nIter, df0, df1, convBase, convSNPj;
  double chisq, pval, phi0, phi1, scale, twologlik0, twologlik1, beta;
  int family = NB;
  int linkR  = *link;
  int adjZR  = *adjZ;
  
  /* pointers to Y, the last column of X, and Z */
  double *pY, *pXlast, *pZ;
  
  /* Work array */
  double *Xb, *fitted0, *fitted1, *fitted2, *resid, *weights, *offsetN;
  
  /* rank of model matrix */
  int rank0, rank1;
    
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
  
  int nY = dims[0];
  int nX = dims[1];
  int nZ = dims[2];
  int N  = dims[3];
  int maxit = dims[4];
  int useOffset = dims[5];

  double P_cut = *RP_cut;

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
  fitted0 = (double *) Calloc(N, double);
  fitted1 = (double *) Calloc(N, double);
  fitted2 = (double *) Calloc(N, double);
  resid   = (double *) Calloc(N, double);
  weights = (double *) Calloc(N, double);
  offsetN = (double *) Calloc(N, double);

  /* point to the last column of X */
  pXlast  = X + N*nX;
  
  if(*trace){
    Rprintf("\n--------------------------------------------------------------\n");
    Rprintf("(nY, nZ, nX, N, adjZ) = (%d, %d, %d, %d, %d)\t", nY, nZ, nX, N, adjZR);
    Rprintf("P_cut=%e\n", P_cut);
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
  fprintf(fo, "GeneRowID\tMarkerRowID\tFamily\tChisq\tPvalue\n");
  
  /***
   * identifify eQTL gene by gene
   */
  
  /* pY is the pointer to gene expression data */
  pY = Y;
  
  for(i=0; i<nY; i++,pY+=N){
    
    if(*trace == 1){
      if(i%100 == 0){ Rprintf("\ni=%d\n", i); }
    }else if(*trace > 1){
      Rprintf("\ni=%d\n", i);
    }
    
    /* **********************************************************
     * fit a baseline model using only the confouding covariates 
     * family is assigned to a value of either Poisson or NB
     * **********************************************************/
    
    dimsNew[1] = nX;
    /** 
     * no initial values, Note this value will be changed 
     * in glmNB after initail iteration 
     */
    dimsNew[3] = 0; 

    convBase = glmNB(dimsNew, &nIter, pY, z1, &linkR, offset, X, convGLM,
                     &rank0, Xb, fitted0, resid, weights, &phi0, &scale, 
                     &df0, &family, &twologlik0, scoreTestP, trace, &beta);
        
    if (!convBase) {
      if(*trace){
        Rprintf("\n  glmEQTL: i=%d, initial model (H0) fail to converge\n", i);
      }

      yFailBaselineModel[i] = 1; 
      continue;
    }
    
    if(*trace > 1){
      Rprintf("\n  glmEQTL: i=%d, finish fitting initial model: ");
      Rprintf("phi0 = %e, ", phi0);
      if (family==2) {
        Rprintf("family=Poisson\n");
      }else if (family==5) {
        Rprintf("family=NB\n");
      }else {
        Rprintf("family=%d\n", family);
      }
    }
    
    /* **********************************************************
     * Now start to fit models with both confounding covariates
     * and each of the SNPs 
     * **********************************************************/
    
    pZ = Z;
    
    for(j=0; j<nZ; j++,pZ+=N){

      if(*cis_only){
        if(eChr[i] != mChr[j]) continue;
        
        pos_diff = fabs(ePos[i] - mPos[j]);
        
        if(pos_diff > *cis_distance) continue;
      }
      
      if(*trace > 1){
        Rprintf("\ni=%d, j=%d\n", i, j);
      }
            
      /* *
       * fill the last column of X by Z[j], genotype of one SNP
       * start with the fitted values of the confounder-only model
       */
      
      for (k=0; k<N; k++) {
        pXlast[k]  = pZ[k];
        fitted1[k] = fitted0[k];
      }
      
      phi1  = phi0;
      scale = 1.0;
      
      /* family has been decided in the previous fitting of baseline model */
      
      if (adjZR) {
        dimsNew[1] = nX;
        dimsNew[3] = 1; /* use initial values */
        
        convSNPj = glmNBlog(dimsNew, &nIter, pY, z1, &linkR, offset, X, conv, convGLM, 
                            &rank1, Xb, fitted1, resid, weights, &phi1, &scale, 
                            &df1, &family, &twologlik1, scoreTestP, trace, &beta,
                            fitted2, offsetN);        
      }else {
        dimsNew[1] = nX + 1;
        dimsNew[3] = 1; /* use initial values */
        
        convSNPj = glmNB(dimsNew, &nIter, pY, z1, &linkR, offset, X, convGLM, 
                         &rank1, Xb, fitted1, resid, weights, &phi1, &scale, 
                         &df1, &family, &twologlik1, scoreTestP, trace, &beta);        
      }
      
      if(convSNPj == 0){
        if(*trace){
          Rprintf("\n  Fail GLM: i=%d, j=%d, adjZR=%d, ", i, j, adjZR);
          if (family==2) {
            Rprintf("family=Poisson\n");
          }else if (family==5) {
            Rprintf("family=NB\n");
          }else {
            Rprintf("family=%d\n", family);
          }          
        }
        continue;
      }
      
      /**
       * it is possible that df0 - df1 != 1
       * in glmFit, df = Nu - x_rank. It is possible that Nu is smaller than sample size
       * due to invalid fitted values for certain glm, e.g., negative values for Poisson
       */
      
      if(df0 - df1 != 1){
        Rprintf("i=%d, j=%d, df0=%d, df1=%d, rank0=%d, rank1=%d\n", i, j, df0, df1, rank0, rank1);
        
        if(df0 - df1 < 0.5) continue;
      }
      
      chisq = twologlik1 - twologlik0;

      if (chisq < -1e-5) {
        error("wrong twologlik! i=%d, j=%d, twologlik=(%.2e, %.2e)\n", 
              i, j, twologlik0, twologlik1);
      }

      if (chisq < 1e-5) { 
        pval = 1.0; 
      }else{
        /* pchisq(double x, double df, int lower_tail, int give_log) */
        pval  = pchisq(chisq, (double)(df0 - df1), 0, 0);
      }
      
      k = (int) (pval / 0.01);
      freqs[k] += 1;
      
      if(pval < P_cut){
        freqs[100] += 1;
        
        /* gene ID and SNP ID */
        fprintf(fo, "%d\t%d\t", i+1, j+1);
        fprintf(fo, "%d\t%.3f\t%.2e\n", family, chisq, pval);
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
  

  /* end time */
  sec_e  = time(NULL);
  if(*trace){
    Rprintf("total time spent is %ld secs\n", sec_e-sec_s);
    Rprintf("\n--------------------------------------------------------------\n");
  }
    
  Free(Xb);
  Free(fitted0);
  Free(fitted1);
  Free(fitted2);
  Free(resid);
  Free(weights);
  Free(offsetN);

  *succeed = 1;
}


/**********************************************************************
 *
 * glmEQTL_max1
 *
 * map eQTL by generalized linear model, Could carry out 
 * cis- only eQTL mapping. 
 *
 * This function is a simplified version of glmEQTL to be used by 
 * function lmEQTL_permute
 *
 * the output is the best (associated) marker and the corresponding 
 * p-value for each gene
 *
 ** 
 ** THIS FUNCTION HAS NOT BEEN FULLY TESTED YET. BETTER USE TRECASE
 ** FOR PERMUTATIONS
 **
 
 Input:
 
 Y   Response variable
 X   Confounding factors
 Z   Covariates of interest
 **********************************************************************/

void glmEQTL_max1 (int* dims, double* Y, double* X, double* Z, double* z1, int* link,  
                  double* offset, int* adjZ, int* cis_only, int* cis_distance, 
                  int* eChr, int* ePos, int* mChr, int* mPos, double* conv, 
                  double* convGLM, int* yFailBaselineModel, double* scoreTestP, 
                  int* best_m, double* pval, int* trace, int* succeed, 
                  double* Xb, double* fitted0, double* fitted1, double* resid, 
                  double* weights, double* fitted2, double* offsetN)
{
  int i, j, k, c, nIter, df0, df1, convBase, convSNPj;
  double chisq, phi0, phi1, scale, twologlik0, twologlik1, beta;
  int family = NB;
  int linkR  = *link;
  int adjZR  = *adjZ;
  
  /* pointers to Y, the last column of X, and Z */
  double *pY, *pXlast, *pZ;
  
  /* rank of model matrix */
  int rank0, rank1;
  
  /* position difference between gene and marker */
  int pos_diff;
  
  /* p-value for each gene */
  double pv1, pv_min=1.0;
  
  /* index of the marker with best association */
  int bestM;
  
  int nY = dims[0];
  int nX = dims[1];
  int nZ = dims[2];
  int N  = dims[3];
  int maxit = dims[4];
  int useOffset = dims[5];
    
  /* dimsNew is passed to function glmNB, and then function glmFit */
  int dimsNew[5];
  dimsNew[0] = N;
  dimsNew[1] = nX;
  dimsNew[2] = maxit;
  dimsNew[3] = 0; /* whetehr to use initial values */
  dimsNew[4] = useOffset;
  
    
  /* point to the last column of X */
  pXlast  = X + N*nX;
  
  if(*trace){
    Rprintf("\n--------------------------------------------------\n");
    Rprintf("(nY, nZ, nX, N) = (%d, %d, %d, %d)\t", nY, nZ, nX, N);
    Rprintf("\n--------------------------------------------------\n");
  }
  
  /* time records */
  time_t sec_s;
  time_t sec_e;
  
  /* starting time */
  sec_s = time(NULL);
  
  /***
   * identifify eQTL gene by gene
   */
  
  /* pY is the pointer to gene expression data */
  pY = Y;
  
  for(i=0; i<nY; i++,pY+=N){
    pv_min = 1.0;
    bestM  = -9;
    
    pval[i]   = pv_min;
    best_m[i] = bestM;

    if(*trace==1){
      if(i%100 == 0){ Rprintf("\ni=%d\n", i); }
    }else if(*trace > 1){
      Rprintf("\ni=%d\n", i);
    }
    
    /* **********************************************************
     * fit a baseline model using only the confouding covariates 
     * family is assigned to a value of either Poisson or NB
     * **********************************************************/
    
    dimsNew[1] = nX;
    /** 
     * no initial values, Note this value  will be changed 
     * in glmNB after initail iteration 
     */
    dimsNew[3] = 0; 
    
    convBase = glmNB(dimsNew, &nIter, pY, z1, &linkR, offset, X, convGLM, 
                     &rank0, Xb, fitted0, resid, weights, &phi0, &scale, 
                     &df0, &family, &twologlik0, scoreTestP, trace, &beta);
    
    if (!convBase) {
      if(*trace){
        Rprintf("\n  The initial model with only confounding covaraites fail to converge\n");
        Rprintf("  This transcript is skipped\n");
      }
      
      yFailBaselineModel[i] = 1; 
      continue;
    }
    
    if(*trace > 1){
      Rprintf("\n  finish fitting initial model with only confounding covaraites,");
      Rprintf("family=%d (2: Poisson, 5: NB)\n", family);
    }
    
    /* **********************************************************
     * Now start to fit models with both confounding covariates
     * and each of the SNPs 
     * **********************************************************/
        
    pZ = Z;
    
    for(j=0; j<nZ; j++,pZ+=N){
      
      if(*cis_only){
        if(eChr[i] != mChr[j]) continue;
        
        pos_diff = fabs(ePos[i] - mPos[j]);
        
        if(pos_diff > *cis_distance) continue;
      }
      
      if(*trace > 1){
        Rprintf("\ni=%d, j=%d\n", i, j);
      }
      
      /* *
       * fill the last column of X by Z[j], genotype of one SNP
       * start with the fitted values of the confounder-only model
       */
      
      for (k=0; k<N; k++) {
        pXlast[k]  = pZ[k];
        fitted1[k] = fitted0[k];
      }
      
      phi1  = phi0;
      scale = 1.0;
      
      /*
       * family has been decided in the previous fitting of baseline model
       */
      
      if (adjZR) {
        dimsNew[1] = nX;
        dimsNew[3] = 1; /* use initial values */
        convSNPj = glmNBlog(dimsNew, &nIter, pY, z1, &linkR, offset, X, conv, convGLM, 
                            &rank1, Xb, fitted1, resid, weights, &phi1, &scale, 
                            &df1, &family, &twologlik1, scoreTestP, trace, &beta, 
                            fitted2, offsetN);
      }else {
        dimsNew[1] = nX + 1;
        dimsNew[3] = 1; /* use initial values */
        convSNPj = glmNB(dimsNew, &nIter, pY, z1, &linkR, offset, X, convGLM, 
                         &rank1, Xb, fitted1, resid, weights, &phi1, &scale, 
                         &df1, &family, &twologlik1, scoreTestP, trace, &beta);
      }
            
      if(convSNPj == 0){
        if(*trace){
          Rprintf("\nFail to fit GLM model for i=%d, j=%d, family=%d (Poission:2, NB:5)\n", i, j, family);
        }
        continue;
      }
      
      /**
       * it is possible that df0 - df1 != 1
       * in glmFit, df = Nu - x_rank. It is possible that Nu is smaller than sample size
       * due to invalid fitted values for certain glm, e.g., negative values for Poisson
       */
      
      if(df0 - df1 != 1){
        Rprintf("i=%d, j=%d, df0=%d, df1=%d, rank0=%d, rank1=%d\n", i, j, df0, df1, rank0, rank1);
        
        if(df0 - df1 < 0.5) continue;
      }
      
      chisq = twologlik1 - twologlik0;
      
      if (chisq < -1e-5) {
        error("wrong twologlik! i=%d, j=%d, twologlik=(%.2e, %.2e)\n", 
              i, j, twologlik0, twologlik1);
      }
      
      if (chisq < 1e-5) { continue; }
      
      /* pchisq(double x, double df, int lower_tail, int give_log)
       */
      pv1  = pchisq(chisq, (double)(df0 - df1), 0, 0);
      
      if(pv1 < pv_min){
        pv_min = pv1;
        bestM  = j; /* here j is markerID - 1*/
      }
    }
    
    pval[i]   = pv_min;
    best_m[i] = bestM;
  }
  
  /* end time */
  sec_e  = time(NULL);
  
  if(*trace){
    Rprintf("total time spent is %ld secs\n", sec_e-sec_s);
    Rprintf("------------------------------------------------------\n");
  }
    
  *succeed = 1;
}


/**********************************************************************
 *
 * glmEQTL_permute
 *
 * map eQTL by generalized linear model, Could carry out 
 * cis- only eQTL mapping. 
 *
 * Carry out fixed number of permutations, and calculate permutation
 * p-value of the best association of each gene expression trait
 *
 
 Input:
 
 Y   Response variable
 X   Confounding factors
 Z   Covariates of interest
 **********************************************************************/


void glmEQTL_permute (int* dims, double* Y, double* X, double* Z, double *z1, 
                      int* npermute, int* link, double* offset, int* adjZ, 
                      int* cis_only, int* cis_distance, int* eChr, int* ePos, 
                      int* mChr, int* mPos, double* conv, double* convGLM, 
                      int* yFailBaselineModel, double* scoreTestP, 
                      int* best_m, double* pval, double* perPval, 
                      int* trace, int* succeed)
{
  
  int p, q, gap;
  time_t timer;

  int nY = dims[0];
  int nX = dims[1];
  int nZ = dims[2];
  int N  = dims[3];
  
  int nPer = *npermute;
  
  int* best_m0, *perm1, *nPer0;
  double* pval0, *perY;
  
  /***
   * variables to be used by glmEQTL_max1
   */  
  
  double *Xb, *fitted0, *fitted1, *fitted2, *resid, *weights, *offsetN;
  
  Xb      = (double *) Calloc(N*(nX+1), double);
  fitted0 = (double *) Calloc(N, double);
  fitted1 = (double *) Calloc(N, double);
  fitted2 = (double *) Calloc(N, double);
  resid   = (double *) Calloc(N, double);
  weights = (double *) Calloc(N, double);
  offsetN = (double *) Calloc(N, double);

  /***
   * best (associated) marker and p-value for each gene 
   * in each permutation
   */  
  best_m0 = (int *)    Calloc(nY, int);
  pval0   = (double *) Calloc(nY, double);
  
  /* permuted expression values */
  perY    = (double *) Calloc(nY*N, double);

  /* number of permutations for each gene */
  nPer0   = (int *) Calloc(nY, double);

  /* pernutation index */
  perm1   = (int *) Calloc(N, int);
  

  /* initial permutation p-value */

  for(q=0; q<nY; q++){
    perPval[q] = 0.0;
    nPer0[q]   = 0;
  }
  
  *succeed = 1;

  /***
   * find the smallest p-value for each gene with unpermuted data 
   */
  timer=time(NULL);

  Rprintf("\nun-permuted data %s\n", asctime(localtime(&timer)));
  
  glmEQTL_max1 (dims, Y, X, Z, z1, link, offset, adjZ, cis_only, cis_distance, 
                eChr, ePos, mChr, mPos, conv, convGLM, yFailBaselineModel, 
                scoreTestP, best_m, pval, trace, succeed, Xb, fitted0, 
                fitted1, resid, weights, fitted2, offsetN);
  
  if(!succeed){
    error("fail glmEQTL_max1\n");
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
    
    if(p % gap == 0){
      timer=time(NULL);
      Rprintf("\n%d-th permutation %s\n", p, asctime(localtime(&timer)));
    }
    
    getPermute(perm1, N);
    // Rprint_vi(perm1, 0, N-1);

    permute(Y, perY, perm1, nY, N);
    
    glmEQTL_max1 (dims, perY, X, Z, z1, link, offset, adjZ, cis_only, cis_distance, 
                  eChr, ePos, mChr, mPos, conv, convGLM, yFailBaselineModel, 
                  scoreTestP, best_m0, pval0, trace, succeed, Xb, fitted0, 
                  fitted1, resid, weights, fitted2, offsetN);
    
    if(!succeed){
      error("fail glmEQTL_max1\n");
    }
    
    for(q=0; q<nY; q++){
      if(best_m0[q] > -0.01){
        nPer0[q] += 1;
        if(pval0[q] <= pval[q]) perPval[q] += 1.0; 
      }
    }
  }
  
  PutRNGstate();
  
  for(q=0; q<nY; q++){
    perPval[q] /= nPer0[q];
  }
  
  /***
   * variables to be used by glmEQTL_max1
   */  
  
  Free(Xb);
  Free(fitted0);
  Free(fitted1);
  Free(fitted2);
  Free(resid);
  Free(weights);
  Free(offsetN);

  Free(best_m0);
  Free(pval0);
  Free(perY);
  Free(nPer0);
  Free(perm1);
  
}
