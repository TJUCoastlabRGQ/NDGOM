#include "mex.h"
#include <math.h>
#include "../../../../../NdgMath/NdgMemory.h"
#include "../../../../../NdgMath/NdgMath.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern double *BaroclinicT, *BaroclinicS, *BaroclinicDTS;

extern char *BaroDensityInitialized;

void MyExit()
{
	if (!strcmp("True", BaroDensityInitialized)) {
		BaroDensityMemoryDeAllocation();
		BaroDensityInitialized = "False";
	}
	return;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	mexAtExit(&MyExit);
	double *height = mxGetPr(prhs[0]);
	double *hT = mxGetPr(prhs[1]);
	double *hS = mxGetPr(prhs[2]);
	double *z = mxGetPr(prhs[3]);
	double hcrit = mxGetScalar(prhs[4]);
	int Np = (int)mxGetScalar(prhs[5]);
	int K = (int)mxGetScalar(prhs[6]);

	plhs[0] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *rho = mxGetPr(plhs[0]);

	if (!strcmp("False", BaroDensityInitialized)) {
		BaroDensityMemoryAllocation(Np, K);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K; i++) {
		DotCriticalDivide(BaroclinicT + i*Np, hT + i*Np, &hcrit, height + i*Np, Np);
		DotCriticalDivide(BaroclinicS + i*Np, hS + i*Np, &hcrit, height + i*Np, Np);
		DotProduct(BaroclinicDTS + i*Np, height + i*Np, z + i*Np, Np);
		MultiplyByConstant(BaroclinicDTS + i*Np, BaroclinicDTS + i*Np, -0.001, Np);
		for (int p = 0; p < Np; p++) {
			rho[i*Np + p] = 999.83 + 5.053*BaroclinicDTS[i*Np + p] - 0.048*BaroclinicDTS[i*Np + p] * BaroclinicDTS[i*Np + p] + \
				(0.808 - 0.0085*BaroclinicDTS[i*Np + p])*BaroclinicS[i*Np + p] - \
				0.0708*(1.0 + 0.351*BaroclinicDTS[i*Np + p] + 0.068*(1 - 0.0683*BaroclinicDTS[i*Np + p])*BaroclinicT[i*Np + p])*BaroclinicT[i*Np + p] - \
				0.003*(1 - 0.059*BaroclinicDTS[i*Np + p] - 0.012*(1 - 0.064*BaroclinicDTS[i*Np + p])*BaroclinicT[i*Np + p])*(35 - BaroclinicS[i*Np + p])*BaroclinicT[i*Np + p];
		}
	}
}