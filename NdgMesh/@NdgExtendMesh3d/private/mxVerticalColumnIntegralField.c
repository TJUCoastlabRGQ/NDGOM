#include <mex.h>
#include "../../../NdgMath/NdgMath.h"

void GetModCoefficient(double *dest, double *source, double *InvV, int Np, int NLayer){
	ptrdiff_t RowOPA = (ptrdiff_t)Np;
	ptrdiff_t ColOPB = (ptrdiff_t)NLayer;
	ptrdiff_t ColOPA = (ptrdiff_t)Np;
	double alpha = 1.0;
	ptrdiff_t LDA = (ptrdiff_t)Np;
	ptrdiff_t LDB = (ptrdiff_t)Np;
	double Beta = 0.0;
	ptrdiff_t LDC = (ptrdiff_t)Np;
	double *B = source;
	dgemm("N", "N", &RowOPA, &ColOPB, &ColOPA, &alpha, InvV, &LDA, B, &LDB, &Beta, dest, &LDC);
}

void GetIntegralValue(double *dest, int Np2d, double *V, double *source){
	ptrdiff_t RowOPA = (ptrdiff_t)Np2d;
	ptrdiff_t ColOPA = (ptrdiff_t)Np2d;
	ptrdiff_t ColOPB = 1;
	double alpha = 1.0;
	ptrdiff_t LDA = (ptrdiff_t)Np2d;
	double *B = source;
	ptrdiff_t LDB = (ptrdiff_t)Np2d;
	double Beta = 0.0;
	ptrdiff_t LDC = (ptrdiff_t)Np2d;
	dgemm("N", "N", &RowOPA, &ColOPB, &ColOPA, &alpha, V, &LDA, B, &LDB, &Beta, dest, &LDC);
	MultiplyByConstant(dest, dest, sqrt(2), Np2d);
}

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
		DotProduct(Tempfield3d + i*Np3d*NLayer, field3d + i*Np3d*NLayer, Jz + i*Np3d*NLayer, Np3d*NLayer);
		GetModCoefficient(fmod + i*Np3d*NLayer, Tempfield3d + i*Np3d*NLayer, InvV3d, Np3d, NLayer);
		for (int L = 0; L < NLayer; L++){
		//	void Add(double *dest, double *sourcea, double *sourceb, int size){
			Add(Tempfield2d + i*Np2d, Tempfield2d + i*Np2d, fmod + i*Np3d*NLayer + L*Np3d, Np2d);
		}
		GetIntegralValue(field2d + i*Np2d, Np2d, V2d, Tempfield2d + i*Np2d);
	}
	free(Tempfield3d);
	free(Tempfield2d);
	free(fmod);
	free(InvV3d);
}