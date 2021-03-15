#include "../NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgQuadratureFreeNonhydrostaticSolver3d/private/SWENonhydrostatic3d.h"


//void FetchDataInSparseMatrix(double *dest, double *src, int NonzeroNum, int Np)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	//The first dimension of the sparse matrix
	double *source = mxGetPr(prhs[0]);

	int NonzeroNum = (int)mxGetScalar(prhs[1]);
	int Np = (int)mxGetScalar(prhs[2]);

	plhs[0] = mxCreateDoubleMatrix(Np, Np, mxREAL);
	double *dest = mxGetPr(plhs[0]);

	FetchDataInSparseMatrix(dest, source, NonzeroNum, Np);

}