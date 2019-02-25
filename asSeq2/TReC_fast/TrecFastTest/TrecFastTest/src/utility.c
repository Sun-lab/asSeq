/**********************************************************************
 * 
 * utility.c
 *
 * copyright (c) 2006, Wei Sun, UCLA
 *
 * last modified Sep 27, 2006
 * first written Aug 29, 2006
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * utilities C functions
 *  
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include "utility.h"

/**********************************************************************
 * 
 * max(x1, x2), min(x1, x2)
 *
 **********************************************************************/
int max_int (const int x1, const int x2){
    if(x1 > x2)
        return x1;
    else
        return x2;
}

int min_int (const int x1, const int x2){
    if(x1 < x2)
        return x1;
    else
        return x2;
}

double max_double (const double x1, const double x2){
  if(x1 > x2)
    return x1;
  else
    return x2;
}

double min_double (const double x1, const double x2){
  if(x1 < x2)
    return x1;
  else
    return x2;
}

/**********************************************************************
 * 
 * max(v, nl, nh, val, idx)
 * in vector x[nl:nh], find the maximum value and the corresponding index
 *
 **********************************************************************/

void max(double *v, int nl, int nh, double *val, int *idx)
{
    int i;
    *idx = nl;
    *val = v[nl];
    
    for(i=nl+1; i<=nh; i++){
        if(v[i] > *val){
            *val = v[i];
            *idx = i;
        }
    }
}

/**********************************************************************
 * 
 * is_infinite(x)
 *
 * return +1 if x is positive infinity, -1 if x is negative infinity
 * and 0 otherwise
 *
 **********************************************************************/
int is_infinite (const double x){
    double y = x - x;
    int s = (y!=y);
    
    if(s && x >0)
        return +1;
    else if(s && x < 0)
        return -1;
    else
        return 0;
}

/**********************************************************************
 * 
 * reorg
 *
 * Reorganize a vector to a matrix of given size. 
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void reorg(double *v, double ***m, int nrow, int ncol)
{
    int i;
    
    *m = (double **)R_alloc(nrow, sizeof(double*));
    
    (*m)[0] = v;
    if(nrow>1){
        for(i=1; i<nrow; i++){
            (*m)[i] = (*m)[i-1] + ncol;
        }
    }
}

void reorg_int(int *v, int ***m, int nrow, int ncol)
{
    int i;
    
    *m = (int **)R_alloc(nrow, sizeof(int*));
    
    (*m)[0] = v;
    if(nrow>1){
        for(i=1; i<nrow; i++){
            (*m)[i] = (*m)[i-1] + ncol;
        }
    }
}


/**********************************************************************
 * 
 * readfile
 *
 * read data into a matrix with given rows and columns
 *
 **********************************************************************/

void readfile(double** mat, char *str, int nrl, int nrh, int ncl, int nch) {
    FILE *file;
    int i,j;
    char temp[255];
    file = fopen(str,"r+t");
    for (i = nrl; i <= nrh; i++) {
        for (j = ncl; j <= nch; j++) {
            fscanf(file,"%s",temp);
            mat[i][j] = (double) atof(temp);
        }
    }
    fclose(file);
}

/**********************************************************************
 * 
 * print_v
 *
 * print out a vector of type double
 *
 **********************************************************************/

void print_v(double* v, int nrl, int nrh)
{
	int i;
	for (i = nrl; i < nrh; i++){
		printf ("%f\t", v[i]);
	}
	printf ("%f\n", v[i]);
}

void Rprint_v(double* v, int nrl, int nrh)
{
	int i;
	for (i = nrl; i < nrh; i++){
		Rprintf ("%f\t", v[i]);
	}
	Rprintf ("%f\n", v[i]);
}

void Rprint_ve(double* v, int nrl, int nrh)
{
	int i;
	for (i = nrl; i < nrh; i++){
		Rprintf ("%.2e\t", v[i]);
	}
	Rprintf ("%.2e\n", v[i]);
}

void Rprint_vi(int* v, int nrl, int nrh)
{
	int i;
	for (i = nrl; i < nrh; i++){
		Rprintf ("%d\t", v[i]);
	}
	Rprintf ("%d\n", v[i]);
}

/**********************************************************************
 * 
 * print_me
 *
 * print out a matrix of format %e
 *
 **********************************************************************/

void print_me(double** m, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i <= nrh; i++){
        for(j = ncl; j <= nch; j++){
            printf ("%.2e\t", m[i][j]);
        }
        printf("\n");
	}
}

void Rprint_me(double** m, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i <= nrh; i++){
        for(j = ncl; j <= nch; j++){
            Rprintf ("%.2e\t", m[i][j]);
        }
        Rprintf("\n");
	}
}

void Rprint_mi(int** m, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i <= nrh; i++){
        for(j = ncl; j <= nch; j++){
            Rprintf ("%i\t", m[i][j]);
        }
        Rprintf("\n");
	}
}

/**********************************************************************
 * 
 * print_mf
 *
 * print out a matrix of format %f
 *
 **********************************************************************/

void print_mf(double** m, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i <= nrh; i++){
        for(j = ncl; j <= nch; j++){
            printf ("%f\t", m[i][j]);
        }
        printf("\n");
	}
}

void Rprint_mf(double** m, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i <= nrh; i++){
        for(j = ncl; j <= nch; j++){
            Rprintf ("%f\t", m[i][j]);
        }
        Rprintf("\n");
	}
}

/**********************************************************************
 * 
 * comp_double, com_double_rev
 *
 * compare two numbers of mode double 
 *
 **********************************************************************/

int comp_double(const double *num1, const double *num2)
{
  if (*num1 <  *num2){
    return -1;
  }else if (*num1 == *num2){
    return  0;
  }else{
    return  1;
  }
}

int comp_double_rev(const double *num1, const double *num2)
{
  if (*num1 >  *num2){
    return -1;
  }else if (*num1 == *num2){
    return  0;
  }else{
    return  1;
  }
}

/**********************************************************************
 * 
 * dnorm
 *
 * normal density, can be log(density)
 *
 **********************************************************************/

double dnorm(double x, double mean, double sd, int logP)
{
  double dx = 0.0;
  double pi = 3.1415926536;
  double log2PI  = log(2*pi);
  double sqrt2PI = sqrt(2*pi);
  if(logP){
    dx = -0.5*(x - mean)*(x - mean)/(sd*sd) - log(sd) - 0.5*log2PI;
  }else{
    dx = exp(-0.5*(x - mean)*(x - mean)/(sd*sd))/(sqrt2PI*sd);
  }
  return(dx);
}

/**********************************************************************
 * 
 * logL
 *
 * log Likelihood, assume normal density
 *
 **********************************************************************/

double logL(double *r, int n, double mean_r, double sd_r){
  int i;
  double l = 0.0;
  
  for(i=0; i<n; i++){
    l += dnorm(r[i], mean_r, sd_r, 1);
  }
  
  return(l);
}

/**********************************************************************
 * 
 * mean
 *
 * take mean value of x[nl:nh]
 *
 **********************************************************************/

double mean(double* x, int nl, int nh)
{
    double sum = 0.0;
    int i;
    for(i=nl; i<=nh; i++){
        sum += x[i];
    }
    sum /= (nh-nl+1);
    return(sum);
}

/**********************************************************************
 * 
 * sd
 *
 * calculate sd of  x[nl:nh], given mean_x
 *
 * use MLE of sd, i.e. 
 * sd(x) = sum(x - mean_x)^2/n instead of /(n-1)
 *
 **********************************************************************/

double sd(double* x, int nl, int nh, double mean_x)
{
  double sum2=0.0, t, sd;
  int i;
  int n = nh - nl + 1;
  for(i=nl; i<=nh; i++){
      t = x[i] - mean_x;
      sum2 += t*t;
  }
  sum2 /= n;
  sd = sqrt(sum2);
  
  if (sd < 1e-10) { sd = 1e-10; }

  return(sd);
}

/**********************************************************************
 * 
 * sd
 *
 * calculate sd of  x[nl:nh], given mean_x
 *
 * use MLE of sd, i.e. 
 * sd(x) = sum(x - mean_x)^2/n instead of /(n-1)
 *
 **********************************************************************/

double var(double* x, int nl, int nh, double mean_x)
{
  double sum2=0.0, t;
  int i;
  int n = nh - nl + 1;
  for(i=nl; i<=nh; i++){
    t = x[i] - mean_x;
    sum2 += t*t;
  }
  sum2 /= n;
    
  return(sum2);
}

/**********************************************************************
 * 
 * mean_sd
 *
 * take mean value and standard deviation of x[nl:nh]
 *
 * use MLE of sd, i.e. 
 * sd(x) = sum(x - mean(x))^2/n instead of /(n-1)
 *
 **********************************************************************/

int mean_sd(double* x, int nl, int nh, double* m_sd)
{
    double sum = 0.0, sum2=0.0;
    int i;
    int n = nh - nl + 1;
    for(i=nl; i<=nh; i++){
        sum += x[i];
        sum2 += x[i]*x[i];
    }
    sum2 /= n;
    sum  /= n;
    m_sd[0] = sum;
    m_sd[1] = sqrt(sum2 - sum*sum);
    return(1);
}

/**********************************************************************
 * 
 * logsumexp
 *
 * log(sum(exp(v)))
 *
 **********************************************************************/

double logsumexp(double* v, int nl, int nh)
{
  int i, idx=0;
  double res, val=0.0;
  
  if(nl > nh){
    error("nl > nh in logsumexp\n");
  }else if(nl == nh){
    val = v[0];
  }else{
    max(v, nl, nh, &val, &idx);
    if(val==1.0/0.0){
      error("positive infinite value in v\n");
    }
    res = 0;
    for(i=nl; i<=nh; i++){
      if(i==idx || v[i]==-1.0/0.0){
        continue;
      }
      res = res + exp(v[i] - val);
    }
    
    val = val + log(1+res);
  }
  return(val);
}

/**********************************************************************
 * 
 * slope
 *
 * slope of linear model lm(y[nl:nh] ~ x[nl:nh])
 *
 **********************************************************************/

double slope(double* nx, double* np, int n)
{
    int i;
    double sum_xp = 0.0;
    double sum_pp = 0.0;
    
    for(i=0; i<n; i++){
        sum_xp += nx[i]*np[i];
        sum_pp += np[i]*np[i];
    }
            
    return(sum_xp/sum_pp);
}


/**********************************************************************
 * 
 * getPermute
 *
 * generate a permutaion from 0 to n-1 using Knuth shuffle
 * http://en.wikipedia.org/wiki/Knuth_shuffle 
 *
 **********************************************************************/

int getPermute(int* per, int n) {
  int i, j, v;
  
  for (i=0; i<n; i++) {
    per[i] = i;
  }
  
  for (i=n-1; i>0; i--) {
    j = floor(unif_rand()*(i+1)); /* a integer from 0 to i */
    
    if(j > i){
      j=i;
      Rprintf("lucky you, your random number hit the boundary exactly!\n");
    }
    
    v = per[i];
    per[i] = per[j];
    per[j] = v;
  }
  
  return 1;
}


/**********************************************************************
 * 
 * permute Y,  which is a matrix of size nY*N, but stored as 
 * a vecter of length nY*N. We need to permute its columns.
 *
 **********************************************************************/

int permute(double* Y, double* perY, int* perm1, int nY, int N) {
  int i, j;
  double *pY, *ppY;
  
  /* pY is the pointer to gene expression data */
  pY  = Y;
  ppY = perY;
  
  for(i=0; i<nY; i++,pY+=N,ppY+=N){
    for(j=0; j<N; j++){
      ppY[j] = pY[perm1[j]];
    }
  }
  
  return 1;
}
