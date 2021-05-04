#include <mex.h>
#include "../../../../../../NdgMath/NdgSWE.h"
#include "../../../../../../NdgMath/NdgSWE3D.h"
#include "../../../../../../NdgMath/NdgMath.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	mxArray *TempNe = mxGetField(prhs[0], 0, "Ne");
	int IENe = (int)mxGetScalar(TempNe);

	mxArray *TempNfp = mxGetField(prhs[0], 0, "Nfp");
	int IENfp = (int)mxGetScalar(TempNfp);

	mxArray *Tempnx = mxGetField(prhs[0], 0, "nx");
	double *IEnx = mxGetPr(Tempnx);

	mxArray *Tempny = mxGetField(prhs[0], 0, "ny");
	double *IEny = mxGetPr(Tempny);

	double *VSIEfm = mxGetPr(prhs[1]);
	double *VSIEfp = mxGetPr(prhs[2]);

	double gra = mxGetScalar(prhs[3]);
	double Hcrit = mxGetScalar(prhs[4]);

	
	plhs[0] = mxCreateDoubleMatrix(IENfp, IENe, mxREAL);
	double *FluxS = mxGetPr(plhs[0]);

	/*
	double *IEhuM2d = VSIEfm2d, *IEhvM2d = VSIEfm2d + IENfp2d * IENe2d, \
		*IEhM2d = VSIEfm2d + 2 * IENfp2d*IENe2d;
	double *IEhuP2d = VSIEfp2d, *IEhvP2d = VSIEfp2d + IENfp2d * IENe2d, \
		*IEhP2d = VSIEfp2d + 2 * IENfp2d*IENe2d;
		*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe; e++){
		GetPCENumericalFluxTerm_HLLC_LAI(FluxS + e*IENfp, VSIEfm + e*IENfp, VSIEfp + e*IENfp, IEnx + e*IENfp, IEny + e*IENfp, &gra, Hcrit, IENfp, IENe);
	}
}