#include "mex.h"
#include "mxSWE3d.h"
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	double *height = mxGetPr(prhs[0]);
	double *hT = mxGetPr(prhs[1]);
	double *hS = mxGetPr(prhs[2]);
	double *z = mxGetPr(prhs[3]);
	double hcrit = mxGetScalar(prhs[4]);
	int Np = (int)mxGetScalar(prhs[5]);
	int K = (int)mxGetScalar(prhs[6]);

	plhs[0] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *rho = mxGetPr(plhs[0]);
	double *T = malloc(Np*K * sizeof(double));
	double *S = malloc(Np*K * sizeof(double));
	/*Depth to surface, in kilometers*/
	double *DTS = malloc(Np*K * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K; i++) {
		DotCriticalDivide(T + i*Np, hT + i*Np, &hcrit, height + i*Np, Np);
		DotCriticalDivide(S + i*Np, hS + i*Np, &hcrit, height + i*Np, Np);
		DotProduct(DTS + i*Np, height + i*Np, z + i*Np, Np);
		MultiplyByConstant(DTS + i*Np, DTS + i*Np, -0.001, Np);
		for (int p = 0; p < Np; p++) {
			rho[i*Np + p] = 999.83 + 5.053*DTS[i*Np + p] - 0.048*DTS[i*Np + p] * DTS[i*Np + p] + \
				(0.808 - 0.0085*DTS[i*Np + p])*S[i*Np + p] - \
				0.0708*(1.0 + 0.351*DTS[i*Np + p] + 0.068*(1 - 0.0683*DTS[i*Np + p])*T[i*Np + p])*T[i*Np + p] - \
				0.003*(1 - 0.059*DTS[i*Np + p] - 0.012*(1 - 0.064*DTS[i*Np + p])*T[i*Np + p])*(35 - S[i*Np + p])*T[i*Np + p];
		}
	}
	free(T);
	free(S);
	free(DTS);
}