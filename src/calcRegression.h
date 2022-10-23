//
//C implementation of Deming regression
//


void calc_Deming(const double *X, const double *Y, int *nX,
				double *error_ratio, double *intercept, double *slope,
				double *se_intercept, double *se_slope, const int *Wmode,
				int *itermax, double *threshold, double *W, double *xw);
void calc_Linreg(const double *X, const double *Y, int *nX,
				double *intercept, double *slope,
				double *se_intercept, double *se_slope,
				double *W, double *xw);
void calc_PaBa(const double *X, const double *Y, const int *nX,
				double *intercept, double *slope,
				double *se_intercept, double *se_slope,
				double *pQuantile, int *pCor, int *tangent, int *Ncpu);
