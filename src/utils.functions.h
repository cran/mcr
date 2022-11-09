#include <R.h>
#include <Rdefines.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>


double quickselect(double *arr, int n, int k);
void get_Wmean(const double* x, const double* w, const int* nX, double* mean, double* sumW);
void get_mean(const double *x, const int *nX, double *mean);
double diff(double x, double y);
void calc_AngleMat(const double* X, const double* Y, const int* N,
					const int* pCor, double* ans, int* nAllItems,
					int* nNeg, int* nNeg2, int* nPos, int* nPos2,
					double* mu, double* sigma);
void calc_AngleMat_opt(const double* X, const double* Y, const int* N,
					const int* pCor, double* ans, int* nAllItems,
					int* nNeg, int* nNeg2, int* nPos, int* nPos2,
					double* mu, double* sigma);
void binapproxR(int* N, double *x, double *ans);
void calc_AngleMat_mod(const double* X, const double* Y, const int* N,
						const int* pCor, double* ans, int* nAllItems,
						int* nNeg, int* nNeg2, int* nPos, int* nPos2,
						double* mu, double* sigma);
