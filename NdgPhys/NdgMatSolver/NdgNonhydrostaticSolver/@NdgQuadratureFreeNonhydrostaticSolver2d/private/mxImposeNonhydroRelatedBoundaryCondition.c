#include "mex.h"
#include <math.h>

typedef enum{
	ZeroGrad = 0,
	Zero     = 1
} NdgNonhydroBoundaryType;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// first, Confirm the location of the corresponding wet to dry face 
	double *EToV = mxGetPr(prhs[0]);
	double *FToV = mxGetPr(prhs[1]);
	double *ZeroFluxBoundary = mxGetPr(prhs[2]);
	double *fm = mxGetPr(prhs[3]);
	double *fp = mxGetPr(prhs[4]);
	signed char *btype = (signed char *)mxGetData(prhs[5]);
	
	
	
	
	
	
}