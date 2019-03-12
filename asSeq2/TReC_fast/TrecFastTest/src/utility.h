/**********************************************************************
 * 
 * utility.h
 *
 * copyright (c) 2006, Wei Sun, UCLA
 *
 * last modified Sep 27, 2006
 * first written Aug 29, 2006
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * utility C functions for the R/ss.hmm package
 *
 * Segmental Semi-Markov Hidden Markov Model
 *
 **********************************************************************/

int max_int (const int x1, const int x2);

int min_int (const int x1, const int x2);

double max_double (const double x1, const double x2);

double min_double (const double x1, const double x2);

void max(double *v, int nl, int nh, double *val, int *idx);

int is_infinite (const double x);

void reorg(double *v, double ***m, int nrow, int ncol);

void reorg_int(int *v, int ***m, int nrow, int ncol);

void readfile(double** mat, char *str, int nrl, int nrh, int ncl, int nch);

void print_v(double* v, int nrl, int nrh);
void Rprint_v(double* v, int nrl, int nrh);
void Rprint_ve(double* v, int nrl, int nrh);
void Rprint_vi(int* v, int nrl, int nrh);

void print_me(double** m, long nrl, long nrh, long ncl, long nch);
void Rprint_me(double** m, long nrl, long nrh, long ncl, long nch);
void Rprint_mi(int** m, long nrl, long nrh, long ncl, long nch);

void print_mf(double** m, long nrl, long nrh, long ncl, long nch);
void Rprint_mf(double** m, long nrl, long nrh, long ncl, long nch);

int comp_double(const double *num1, const double *num2);

int comp_double_rev(const double *num1, const double *num2);

double dnorm(double x, double mean, double sd, int logP);

double logL(double *r, int n, double mean_r, double sd_r);

double mean(double* x, int nl, int nh);

double sd(double* x, int nl, int nh, double mean_x);

double var(double* x, int nl, int nh, double mean_x);

int mean_sd(double* x, int nl, int nh, double* m_sd);

double logsumexp(double* vv, int nl, int nh);

double slope(double* nx, double* np, int n);

int getPermute(int* per, int n);

int permute(double* Y, double* perY, int* perm1, int nY, int N);
