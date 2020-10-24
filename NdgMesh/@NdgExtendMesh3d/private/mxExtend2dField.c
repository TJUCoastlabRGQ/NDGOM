#include "..\..\..\NdgMath\NdgMath.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int K2d = (int)mxGetScalar(prhs[0]);
	double *source = mxGetPr(prhs[1]);
	mwSize Np2d = mxGetM(prhs[1]);
	int NLayer = (int)mxGetScalar(prhs[2]);
	int Np3d = (int)mxGetScalar(prhs[3]);
	int K3d = (int)mxGetScalar(prhs[4]);
	int Nz = (int)mxGetScalar(prhs[5]);
	plhs[0] = mxCreateDoubleMatrix(Np3d, K3d, mxREAL);
	double *dest = mxGetPr(plhs[0]);
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d; i++){
		NdgExtend2dField(dest, source, Np2d, i, Np3d, NLayer, Nz);
	}
}