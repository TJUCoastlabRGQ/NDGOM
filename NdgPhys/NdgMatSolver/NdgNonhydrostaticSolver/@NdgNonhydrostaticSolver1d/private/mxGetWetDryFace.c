#include "../../../../../Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.h"

# include "mex.h"
# include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 2
#define NLHS 1

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	    /* check input & output */
    if (nrhs != NRHS) {
        mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
        mexPrintf("%d inputs required.\n", NRHS);
    }
    
    if (nlhs != NLHS) {
        mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
        mexPrintf("%d inputs required.\n", NLHS);
    }
    
     signed char* regType = (signed char*)mxGetData(prhs[0]);
     double *FToE = mxGetPr(prhs[1]);
     mwSize Ne = mxGetN(prhs[1]);
	
    plhs[0] = mxCreateDoubleMatrix(1, Ne, mxREAL);
    double *DryFaceFlag = mxGetPr(plhs[0]);	
    
    evaluateWetDryInterface( regType, prhs[1], DryFaceFlag);
}
