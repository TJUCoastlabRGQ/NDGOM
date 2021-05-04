#include <mex.h>
#include "../../../../../../NdgMath/NdgSWE.h"
#include "../../../../../../NdgMath/NdgSWE3D.h"
#include "../../../../../../NdgMath/NdgMath.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	mxArray *TempNe = mxGetField(prhs[0], 0, "Ne");
	int BENe = (int)mxGetScalar(TempNe);

	mxArray *TempNfp = mxGetField(prhs[0], 0, "Nfp");
	int BENfp = (int)mxGetScalar(TempNfp);

	mxArray *Tempnx = mxGetField(prhs[0], 0, "nx");
	double *BEnx = mxGetPr(Tempnx);

	mxArray *Tempny = mxGetField(prhs[0], 0, "ny");
	double *BEny = mxGetPr(Tempny);

	double *VSBEfm = mxGetPr(prhs[1]);
	double *VSBEfp = mxGetPr(prhs[2]);

	double gra = mxGetScalar(prhs[3]);
	double Hcrit = mxGetScalar(prhs[4]);

	double *VSBEzM2d = mxGetPr(prhs[5]);

	double *VSBEzP2d = mxGetPr(prhs[6]);

	signed char *ftype2d = (signed char *)mxGetData(prhs[7]);

	double *fext2d = mxGetPr(prhs[8]);

	
	plhs[0] = mxCreateDoubleMatrix(BENfp, BENe, mxREAL);
	double *FluxS = mxGetPr(plhs[0]);

	int Nfield = 3;

	/*
	double *BEhuM2d = VSBEfm2d, *BEhvM2d = VSBEfm2d + BENe2d * BENfp2d, \
		*BEhM2d = VSBEfm2d + 2 * BENe2d * BENfp2d;
	*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe; e++){
		NdgEdgeType type = (NdgEdgeType)ftype2d[e];  // boundary condition
		ImposeBoundaryCondition(&gra, type, BEnx + e*BENfp, BEny + e*BENfp, VSBEfm + e*BENfp, VSBEfp + e*BENfp, \
			VSBEzM2d + e*BENfp, VSBEzP2d + e*BENfp, fext2d + e*BENfp, BENfp, Nfield, BENe);
		EvaluateHydroStaticReconstructValue(Hcrit, VSBEfm + e*BENfp, VSBEfp + e*BENfp, VSBEzM2d + e*BENfp, VSBEzP2d + e*BENfp, BENfp, Nfield, BENe);
		GetPCENumericalFluxTerm_HLLC_LAI(FluxS + e*BENfp, VSBEfm + e*BENfp, VSBEfp + e*BENfp, BEnx + e*BENfp, BEny + e*BENfp, &gra, Hcrit, BENfp, BENe);
	}
}