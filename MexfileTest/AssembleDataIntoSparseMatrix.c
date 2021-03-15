#include "../NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgQuadratureFreeNonhydrostaticSolver3d/private/SWENonhydrostatic3d.h"

//void AssembleContributionIntoSparseMatrix(double *dest, double *src, int NonzeroNum, int Np);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	//The first dimension of the sparse matrix
	int Row = (int)mxGetScalar(prhs[0]);
	int InsertRow = (int)mxGetScalar(prhs[1]);
	//The data to be insert into the sparse matrix
	double *data = mxGetPr(prhs[2]);
	int TotalNonzero = mxGetNumberOfElements(prhs[2]);
	//The sparse matrix
	mwIndex *jcPNPS = mxGetJc(prhs[3]);
	mwIndex *irPNPS = mxGetIr(prhs[3]);

	plhs[0] = mxCreateSparse(Row, InsertRow, TotalNonzero, mxREAL);
	double *sr = mxGetPr(plhs[0]);
	mwIndex *irs = mxGetIr(plhs[0]);
	memcpy(irs, irPNPS, TotalNonzero*sizeof(mwIndex));
	mwIndex *jcs = mxGetJc(plhs[0]);
	memcpy(jcs, jcPNPS, (InsertRow + 1)*sizeof(mwIndex));
	AssembleContributionIntoSparseMatrix(sr, data, InsertRow, InsertRow);
}