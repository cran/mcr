
#include "utils.functions.h" /* for declarations for registration */

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;


double quickselect(double *arr, int n, int k)
{
	unsigned long i, ir, j, l, mid;
	double a, temp;

	l = 0;
	ir = n-1;
	for(;;) {
		if (ir <= l+1) {
			if (ir == l+1 && arr[ir] < arr[l]) {
				SWAP(arr[l],arr[ir]);
			}
			return arr[k];
		} else {
			mid=(l+ir) >> 1;
			SWAP(arr[mid],arr[l+1]);
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir]);
			}
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir]);
			}
			if (arr[l] > arr[l+1]) {
				SWAP(arr[l],arr[l+1]);
			}
			i = l+1;
			j = ir;
			a = arr[l+1];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j]);
			}
			arr[l+1] = arr[j];
			arr[j] = a;
			if (j >= k) ir = j-1;
			if (j <= k) l = i;
		}
	}
}




void get_Wmean(const double* x, const double* w, const int* nX, double* mean, double* sumW)
{
	// calc mean of a vector
	*mean = 0, *sumW = 0;
	int i;
	for(i = 0; i < *nX; i++){
		*mean += w[i] * x[i];
		*sumW += w[i];
	}
	*mean =  *mean / *sumW;
}


void get_mean(const double* x, const int* nX, double* mean)
{
	// calc mean of a vector
	*mean = 0;
	int i;
	for(i = 0; i < *nX; i++){
		*mean += x[i];
	}
	*mean =  *mean / *nX;
}



double diff(double x, double y)
{
	double eps = 1.e-12;
	double dRes = x - y;
	if(dRes == 0.) return 0.;
	if(fabs(dRes) < eps*(fabs(x) + fabs(y))/2.) return 0.;
	return dRes;
}




void calc_AngleMat(const double* X, const double* Y, const int* N,
					const int* pCor, double* ans, int* nAllItems,
					int* nNeg, int* nNeg2, int* nPos, int* nPos2,
					double* mu, double* sigma)
{
	// Calculate results
	int n = N[0];
	int item = 0, cnNeg = 0, cnNeg2 = 0, cnPos = 0, cnPos2 = 0;
	double dx, dy;
	int j, k;
	for(j = 0; j < n; j++) {
		for(k = j; k < n; k++){
		//for(k = 0; k < n; k++){
			//if(k > j) {
				dx = diff(X[k], X[j]);
				dy = diff(Y[k], Y[j]);
				if(dx != 0.){
					ans[j+n*k] = atan(dy/dx);
					*mu += ans[j+n*k];
					item++;
					if(ans[j+n*k] <= -M_PI_4){
						if(ans[j+n*k] < -M_PI_4){
							cnNeg2++;
						}else{
							cnNeg++;
						}
					}
					if(ans[j+n*k] >= M_PI_4){
						if(ans[j+n*k] > M_PI_4){
							cnPos2++;
						}else{
							cnPos++;
						}
					}
				}
				else if(dy != 0.) {
					// x==0, y!=0
					// only positive infinity for pos correlated
					if(pCor[0]){
						ans[j+n*k] = M_PI_2;
						cnPos2++;
					}
					else{
						ans[j+n*k] = -M_PI_2;
						cnNeg2++;
					}
					*mu += ans[j+n*k];
					item++;
				}
				else ans[j+n*k] = NA_REAL; // dx==0, dy==0, set to NA
			}
			//else ans[j+n*k] = NA_REAL; // not upper triangle, set to NA
		//}
	}
	nAllItems[0] = item;
	nNeg[0] = cnNeg;
	nNeg2[0] = cnNeg2;
	nPos[0] = cnPos;
	nPos2[0] = cnPos2;
	mu[0] = mu[0]/item;

}






void calc_AngleMat_opt(const double* X, const double* Y, const int* N,
					const int* pCor, double* ans, int* nAllItems,
					int* nNeg, int* nNeg2, int* nPos, int* nPos2,
					double* mu, double* sigma)
{
	// Calculate results
	int n = N[0];
	int cnNeg = 0, cnNeg2 = 0, cnPos = 0, cnPos2 = 0;
	double dx,dy;
	int j, k, ansN = 0;

	for(j = 0; j < n; j++) {
		for(k = j; k < n; k++){
			dx = diff(X[k], X[j]);
			dy = diff(Y[k], Y[j]);
			if(dx != 0.){
				ans[ansN] = atan(dy/dx);
				*mu += ans[ansN];
				if(ans[ansN] <= -M_PI_4){
					cnNeg++;
					if(ans[ansN] < -M_PI_4){
						cnNeg2++;
					}
				}
				if(ans[ansN] >= M_PI_4){
					cnPos++;
					if(ans[ansN] > M_PI_4){
						cnPos2++;
					}
				}
				ansN++;
			}else if(dy != 0.) {
				// x==0, y!=0
				// only positive infinity for pos correlated
				if(pCor[0]){
					ans[ansN] = M_PI_2;
					cnPos++;
					cnPos2++;
				}else{
					ans[ansN] = -M_PI_2;
					cnNeg++;
					cnNeg2++;
				}
				*mu += ans[ansN];
				ansN++;
			}
		}
	}
	nAllItems[0] = ansN;
	nNeg[0] = cnNeg;
	nNeg2[0] = cnNeg2;
	nPos[0] = cnPos;
	nPos2[0] = cnPos2;
	mu[0] = mu[0]/ansN;
}






void binapproxR(int* N, double *x, double *ans)
{
	int n = N[0];
	// Compute the mean and standard deviation
	double sum = 0;
	int i;
	for (i = 0; i < n; i++) {
		sum += x[i];
	}
	double mu = sum/n;

	sum = 0;
	for (i = 0; i < n; i++) {
		sum += (x[i]-mu)*(x[i]-mu);
	}
	double sigma = sqrt(sum/n);

	// Bin x across the interval [mu-sigma, mu+sigma]
	int bottomcount = 0;
	int bincounts[1001];
	for (i = 0; i < 1001; i++) {
		bincounts[i] = 0;
	}
	double scalefactor = 1000/(2*sigma);
	double leftend =  mu-sigma;
	double rightend = mu+sigma;
	int bin;

	for (i = 0; i < n; i++) {
		if (x[i] < leftend) {
			bottomcount++;
		}
		else if (x[i] < rightend) {
			bin = (int)((x[i]-leftend) * scalefactor);
			bincounts[bin]++;
		}
	}

	// If n is odd
	if (n & 1) {
		// Find the bin that contains the median
		int k = (n+1)/2;
		int count = bottomcount;

		for (i = 0; i < 1001; i++) {
			count += bincounts[i];

			if (count >= k) {
				//return (i+0.5)/scalefactor + leftend;
				ans[0] = (i+0.5)/scalefactor + leftend;
			}
		}
	}
	// If n is even
	else {
		// Find the bins that contains the medians
		int k = n/2;
		int count = bottomcount;

		for (i = 0; i < 1001; i++) {
			count += bincounts[i];

			if (count >= k) {
				int j = i;
				while (count == k) {
					j++;
					count += bincounts[j];
				}
				//return (i+j+1)/(2*scalefactor) + leftend;
				ans[0] = (i+j+1)/(2*scalefactor) + leftend;
			}
		}
	}
}






void calc_AngleMat_mod(const double* X, const double* Y, const int* N,
						const int* pCor, double* ans, int* nAllItems,
						int* nNeg, int* nNeg2, int* nPos, int* nPos2,
						double* mu, double* sigma)
{
	// Calculate results
	int n = N[0];
	int cnNeg = 0, cnNeg2 = 0, cnPos = 0, cnPos2 = 0;
	double dx, dy, tmpX;
	int j, k, ansN = 0;

	// Compute the mean
	for(j = 0; j < n; j++) {
		for(k = j; k < n; k++){
			dx = diff(X[k], X[j]);
			dy = diff(Y[k], Y[j]);
			if(dx != .0){
				*mu += atan(dy/dx);
				ansN++;
			}else if( dy != 0. ) {
				if(pCor[0]){
					*mu += M_PI_2;
				}else{
					*mu += -M_PI_2;
				}
				ansN++;
			}
		}
	}
	mu[0] = mu[0]/ansN;

	// Compute the standard deviation
	tmpX = 0.0;
	ansN = 0;
	for(j = 0; j < n; j++) {
			for(k = j; k < n; k++){
				dx = diff(X[k], X[j]);
				dy = diff(Y[k], Y[j]);
				if(dx != .0){
					tmpX = atan(dy/dx);
					*sigma += (tmpX-mu[0])*(tmpX-mu[0]);
					ansN++;
				}else if( dy != 0. ) {
					if(pCor[0]){
						*sigma += (M_PI_2-mu[0])*(M_PI_2-mu[0]);
					}else{
						*sigma += (-M_PI_2-mu[0])*(-M_PI_2-mu[0]);
					}
					ansN++;
				}
			}
		}
	sigma[0] = sqrt(sigma[0]/ansN);
	nAllItems[0] = ansN;


	// Calculate results
	n = N[0];
	ansN = 0;

	for(j = 0; j < n; j++) {
		for(k = j; k < n; k++){
			dx = diff(X[k], X[j]);
			dy = diff(Y[k], Y[j]);
			if(dx != .0){
				ans[ansN] = atan(dy/dx);
				if(ans[ansN] <= -M_PI_4){
					cnNeg++;
					if(ans[ansN] < -M_PI_4){
						cnNeg2++;
					}
				}
				if(ans[ansN] >= M_PI_4){
					cnPos++;
					if(ans[ansN] > M_PI_4){
						cnPos2++;
					}
				}
				ansN++;
			}else if( dy != 0. ) {
				// x==0, y!=0
				// only positive infinity for pos correlated
				if(pCor[0]){
					ans[ansN] = M_PI_2;
					cnPos++;
					cnPos2++;
				}else{
					ans[ansN] = -M_PI_2;
					cnNeg++;
					cnNeg2++;
				}
				ansN++;
			}
		}
	}
	nAllItems[0] = ansN;
	nNeg[0] = cnNeg;
	nNeg2[0] = cnNeg2;
	nPos[0] = cnPos;
	nPos2[0] = cnPos2;
}





