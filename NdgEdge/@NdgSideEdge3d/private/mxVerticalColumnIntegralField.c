#include "mex.h"
#include "../../../NdgMath/NdgSWE3D.h"
#include "../../../NdgMath/NdgMath.h"

/*This function is mainly used for testing purpose*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	const mxArray *edge = prhs[0]; 
	mxArray *TempV1d = mxGetField(edge, 0, "V1d");
	double *V1d = mxGetPr(TempV1d);
	mxArray *TempV2d = mxGetField(edge, 0, "V2d");
	double *V2d = mxGetPr(TempV2d);
	int LNfp2d = (int)mxGetM(TempV1d);
	int Nfp2d = (int)mxGetM(TempV2d);
	mxArray *TempNe = mxGetField(edge, 0, "Ne");
	int Ne = (int)mxGetScalar(TempNe);
	double *field3d = mxGetPr(prhs[1]);
	mxArray *TempNLayer = mxGetField(edge, 0, "Nz");
	int NLayer = (int)mxGetScalar(TempNLayer);
	mxArray *TempFToF = mxGetField(edge, 0, "FToF");
	double *FToF = mxGetPr(TempFToF);
	mxArray *TempJz = mxGetField(edge, 0, "Jz");
	double *Jz = mxGetPr(TempJz);

	int Ne2d = Ne / NLayer;
	
	plhs[0] = mxCreateDoubleMatrix(LNfp2d, Ne2d, mxREAL);
	double *OutAverage = mxGetPr(plhs[0]);

	double *InvV2d = malloc(Nfp2d*Nfp2d*sizeof(double));
	memcpy(InvV2d, V2d, Nfp2d*Nfp2d*sizeof(double));
	MatrixInverse(InvV2d, (ptrdiff_t)Nfp2d);

	double *fmod = malloc(Ne2d*Nfp2d*sizeof(double));
	memset(fmod, 0, Ne2d*Nfp2d*sizeof(double));

//void VerticalFaceColumnIntegral(double *dest, double *source, double *fmod, double *InvV2d, int Nfp2d, double *Jz, int Nlayer, double *V1d, int LNfp2d, int FToF)
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < Ne2d; e++){
		VerticalFaceColumnIntegral(OutAverage + e*LNfp2d, field3d + e*NLayer*Nfp2d, fmod + e*Nfp2d, InvV2d, (ptrdiff_t)Nfp2d, Jz + e*NLayer*Nfp2d, NLayer, V1d, (ptrdiff_t)LNfp2d, (int)(*(FToF + e*NLayer * 2)));
	}

	free(InvV2d);
	free(fmod);
}