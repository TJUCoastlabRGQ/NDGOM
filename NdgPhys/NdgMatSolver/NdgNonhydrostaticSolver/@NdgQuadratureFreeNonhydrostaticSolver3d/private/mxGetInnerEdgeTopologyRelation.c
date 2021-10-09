#include "SWENonhydrostatic3d.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	double *FToE = mxGetPr(prhs[0]);

	double *FToF = mxGetPr(prhs[1]);

	double *nx = mxGetPr(prhs[2]);

	double *ny = mxGetPr(prhs[3]);

	int Ne = (int)mxGetScalar(prhs[4]);

	int Nface = (int)mxGetScalar(prhs[5]);

	int Nface2d = Nface - 2;

	int Nfp = (int)mxGetScalar(prhs[6]);

	double *EToE = mxGetPr(prhs[7]);

	int K = (int)mxGetScalar(prhs[8]);

	plhs[0] = mxCreateDoubleMatrix(1,K,mxREAL);
	double *OutNum = mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(Nface + 1, K, mxREAL);
	double *UniElement = mxGetPr(plhs[1]);
	memset(UniElement, 0, (Nface + 1)*K * sizeof(double));

	plhs[2] = mxCreateDoubleMatrix(Nface2d*Nfp, K, mxREAL);
	double *SortedIEnx = mxGetPr(plhs[2]);
	memset(SortedIEnx, 0, Nface2d*Nfp*K*sizeof(double));

	plhs[3] = mxCreateDoubleMatrix(Nface2d*Nfp, K, mxREAL);
	double *SortedIEny = mxGetPr(plhs[3]);
	memset(SortedIEny, 0, Nface2d*Nfp*K * sizeof(double));

	plhs[4] = mxCreateDoubleMatrix(Nface2d, K, mxREAL);
	double *SortedIEGlobalFace = mxGetPr(plhs[4]);
	memset(SortedIEGlobalFace, 0, Nface2d*K * sizeof(double));

	plhs[5] = mxCreateDoubleMatrix(Nface2d, K, mxREAL);
	double *SortedIEAdjEle = mxGetPr(plhs[5]);
	memset(SortedIEAdjEle, 0, Nface2d*K * sizeof(double));

	plhs[6] = mxCreateDoubleMatrix(Nface2d, K, mxREAL);
	double *SortedIEReverseFlag = mxGetPr(plhs[6]);
	memset(SortedIEReverseFlag, 0, Nface2d*K * sizeof(double));

	plhs[7] = mxCreateDoubleMatrix(1, K, mxREAL);
	double *SortedIEInternalFace = mxGetPr(plhs[7]);
	memset(SortedIEInternalFace, 0, K * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < K; e++) {

		FindUniqueElementAndSortOrder(UniElement + e*(Nface + 1), EToE + e*Nface, OutNum + e, Nface, e + 1);

		FindFaceAndDirectionVector(SortedIEnx + e*Nface2d*Nfp, SortedIEGlobalFace + e*Nface2d, SortedIEAdjEle + e*Nface2d, \
			SortedIEInternalFace + e, SortedIEReverseFlag + e*Nface2d, Nfp, e + 1, FToE, FToF, nx, Ne, Nface2d);

		FindFaceAndDirectionVector(SortedIEny + e*Nface2d*Nfp, SortedIEGlobalFace + e*Nface2d, SortedIEAdjEle + e*Nface2d, \
			SortedIEInternalFace + e, SortedIEReverseFlag + e*Nface2d, Nfp, e + 1, FToE, FToF, ny, Ne, Nface2d);
	}

}