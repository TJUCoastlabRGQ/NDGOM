#include "mex.h"
#include "../../../NdgMath/NdgSWE3D.h"
#include "../../../NdgMath/NdgMath.h"

/*This function is mainly used for testing purpose*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	const mxArray *mesh = prhs[0];
	mxArray *TempJz = mxGetField(mesh, 0, "Jz");
	double *Jz = mxGetPr(TempJz);
	mxArray *TempNLayer = mxGetField(mesh, 0, "Nz");
	int NLayer = (int)mxGetScalar(TempNLayer);
	mxArray *TempK = mxGetField(mesh, 0, "K");
	int K = (int)mxGetScalar(TempK);

	const mxArray *cell = prhs[1];
	mxArray *TempNpz = mxGetField(cell, 0, "Npz");
	int Npz = (int)mxGetScalar(TempNpz);
	mxArray *TempNp2d = mxGetField(cell, 0, "Nph");
	int Np2d = (int)mxGetScalar(TempNp2d);
	mxArray *TempNp = mxGetField(cell, 0, "Np");
	int Np = (int)mxGetScalar(TempNp);
	mxArray *TempVint = mxGetField(cell, 0, "Vint");
	double *Vint = mxGetPr(TempVint);
	mxArray *TempV = mxGetField(cell, 0, "V");
	double *V = mxGetPr(TempV);

	double *InvV3d = malloc(Np*Np*sizeof(double));
	memcpy(InvV3d, V, Np*Np*sizeof(double));
	MatrixInverse(InvV3d, (ptrdiff_t)Np);

	const mxArray *mesh2d = prhs[2];
	mxArray *TempK2d = mxGetField(mesh2d, 0, "K");
	int K2d = (int)mxGetScalar(TempK2d);

	double *field3d = mxGetPr(prhs[3]);


	plhs[0] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *DepIntegral = mxGetPr(plhs[0]);

	double *fmod = malloc(K2d*Np*sizeof(double));

	//void VerticalFaceColumnIntegral(double *dest, double *source, double *fmod, double *InvV2d, int Nfp2d, double *Jz, int Nlayer, double *V1d, int LNfp2d, int FToF)
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		VerticalIntegralFromBottom(DepIntegral + k*NLayer*Np, field3d + k*NLayer*Np, Jz + k*NLayer*Np, fmod + k*Np, NLayer, (ptrdiff_t)Np, InvV3d, Np2d, Npz, Vint);
	}

	free(InvV3d);
	free(fmod);
}