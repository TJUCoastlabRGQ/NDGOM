#include <mex.h>
#include "../../../NdgMath/NdgSWE3D.h"
#include "../../../NdgMath/NdgMath.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int Np2d = (int)mxGetScalar(prhs[0]);
	int K2d = (int)mxGetScalar(prhs[1]);
	double *V2d = mxGetPr(prhs[2]);
	double *V3d = mxGetPr(prhs[3]);
	double *Jz = mxGetPr(prhs[4]);
	double *field3d = mxGetPr(prhs[5]);
	int NLayer = (int)mxGetScalar(prhs[6]);
	int Np3d = (int)mxGetScalar(prhs[7]);
	int K3d = (int)mxGetScalar(prhs[8]);

	plhs[0] = mxCreateDoubleMatrix(Np2d, K2d, mxREAL);
	double *field2d = mxGetPr(plhs[0]);
	double *Tempfield2d = malloc(Np2d*K2d*sizeof(double));
	memset(Tempfield2d, 0, Np2d*K2d*sizeof(double));
	double *Tempfield3d = malloc(Np3d*K3d*sizeof(double)); 
	double *fmod = malloc(Np3d*K3d*sizeof(double));
	double *InvV3d = malloc(Np3d*Np3d*sizeof(double));
	memcpy(InvV3d, V3d, Np3d*Np3d*sizeof(double));
	MatrixInverse(InvV3d, (ptrdiff_t)Np3d);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d; i++){
		VerticalColumnIntegralField3d(field2d + i*Np2d, Np2d, V2d, Tempfield2d + i*Np2d, \
			Tempfield3d + i*Np3d*NLayer, field3d + i*Np3d*NLayer, Jz + i*Np3d*NLayer, \
			fmod + i*Np3d*NLayer, InvV3d, Np3d, NLayer);
	}
	free(Tempfield3d);
	free(Tempfield2d);
	free(fmod);
	free(InvV3d);
}