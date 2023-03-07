#include <mex.h>
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *VertNode = mxGetPr(prhs[0]);
	mwSize Nv = mxGetN(prhs[0]);
	mwSize maxPoint = mxGetM(prhs[0]);
	double *NodeWeight = mxGetPr(prhs[1]);
	double *Data = mxGetPr(prhs[2]);
	mwSize Np = mxGetM(prhs[2]);
	mwSize K = mxGetN(prhs[2]);
	plhs[0] = mxCreateDoubleMatrix(Np, K, mxREAL);
	double *dest = mxGetPr(plhs[0]);
	memcpy(dest, Data, Np*K * sizeof(double));
	int i,j;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS) private(j)
#endif
	for (i = 0; i < Nv; i++) {
		double TempNodeData = 0.0;
		//Get the weighted data
		for (j = 0; j < (int)VertNode[i*maxPoint]; j++) {
			TempNodeData += Data[(int)VertNode[i*maxPoint + j+1]-1] * NodeWeight[i*(maxPoint-1) + j];
		}
		//Replace the Data
		for (j = 0; j < (int)VertNode[i*maxPoint]; j++) {
			dest[(int)VertNode[i*maxPoint + j + 1]-1] = TempNodeData;
		}
	}
}