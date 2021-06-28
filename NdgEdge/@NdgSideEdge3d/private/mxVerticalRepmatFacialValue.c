#include "mex.h"
#include "../../../NdgMath/NdgSWE3D.h"
#include "../../../NdgMath/NdgMath.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	const mxArray *edge = prhs[0];
	mxArray *TempNfp = mxGetField(edge, 0, "Nfp");
	int Nfp = (int)mxGetScalar(TempNfp);
	mxArray *TempNe = mxGetField(edge, 0, "Ne");
	int Ne = (int)mxGetScalar(TempNe);
	mxArray *TempNLayer = mxGetField(edge, 0, "Nz");
	int NLayer = (int)mxGetScalar(TempNLayer);
	mxArray *TempFToF = mxGetField(edge, 0, "FToF");
	double *FToF = mxGetPr(TempFToF);

	const mxArray *cell = prhs[1];
	mxArray *TempNpz = mxGetField(cell, 0, "Npz");
	int Npz = (int)mxGetScalar(TempNpz);
	int LNfp = Nfp / Npz;

	double *field2d = mxGetPr(prhs[2]);
	int Ne2d = (int)mxGetN(prhs[2]);

	plhs[0] = mxCreateDoubleMatrix(Nfp, Ne, mxREAL);
	double *DepthExtendValue = mxGetPr(plhs[0]);
	//void VerticalRepmatFacialValue(double *dest, double *source, int NLayer, int Nfp, int LNfp, int Npz, int FToF)
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < Ne2d; e++){
		VerticalRepmatFacialValue(DepthExtendValue + e*NLayer*Nfp, field2d + e*LNfp, NLayer, Nfp, LNfp, Npz, (int)(*(FToF + 2 * e*NLayer)));
	}
}