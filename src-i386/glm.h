/*
 *  glm.h
 *  
 *
 *  Created by Wei Sun on 5/11/10.
 *  Copyright 2010 UNC. All rights reserved.
 *
 */
#ifndef _GLM_H_
#define _GLM_H_

#define BINOMIAL  1
#define POISSON   2
#define GAUSSIAN  3
#define GAMMA     4
#define NB        5
#define ZINB      6

/* Link */

#define LOGIT     1
#define LOG       2
#define IDENTITY  3
#define INVERSE   4

/* GLM definition functions */

int     wcenter(double*, int, double*, int, double*);
int     wresid(double*,  int, double*, double*, double*, double*);
double  wssq(double *y, int n, double *weights);
double  varfun(int, double, double);
int     muvalid(int, double);
double  linkfun(int, double);
double  invlink(int, double);
double  dlink(int, double);
void     initialize(int, double*, double*, int, double*);

/* Fit a base model */

int glmFit(int* familyR, int* linkR, int* dims, int* nIter,
           double *y, double *offset, double *z,
           double *X, double *nTotal_binom, double *convR, 
           int *rank, double *Xb, double *fitted, double *resid, 
           double *weights, double *phi, int* trace, 
           double *scale, int *df_resid, double* beta);

/*  log likelihood of Poisson */
double loglik_Poisson(int N, double* mu, double* y);

/* log likelihood of negative binomial */

double loglik_NB(int N, double phi, double* mu, double* y);

/* Score test for additional terms */

void glm_score_test(int* dims, double *Z, double *resid, 
                    double *weights, double *Xb, double* scaleR,
                    double *chi2, int *df);

/* score and infor for solving MLE of phi */

void score_info(int N, double theta, double* mu, double *y, 
                double* score, double* info);

/* MLE of phi */

int phi_ml(double* y, double* mu, int N,  
           int limit, double eps, double* phi, int initPhi, int trace);

/* glmNB */

int glmNB(int *dims, int *nIter, double *y, double *z, 
          int *linkR, double *offset, double *X, double *convR, 
          int *rank, double *Xb, double *fitted, double *resid, 
          double *weights, double *phi, double *scale, 
          int *df_resid, int* family, double *twologlik, 
          double *scoreTestP, int *trace, double *beta);

int glmNBlog(int *dimsNew, int *nIter, double *pY, double *z, 
             int *linkR, double *offset, double *pX, double *conv, 
             double *convGLM, int *rank, double *Xb, double *fitted, 
             double *resid, double *weights, double *phi, double *scale, 
             int *df_resid, int* family, double *twoLL_trec, 
             double *scoreTestP, int *trace, double *beta,
             double *fitted2, double *offsetN);

#endif
