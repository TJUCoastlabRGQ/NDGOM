#include "Date.h"

#define NLHS 0
#define NRHS 2

int K2d, Np2d;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* check input & output */
	if (nrhs != NRHS) {
		mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
		mexPrintf("%d inputs required.\n", NRHS);
	}

	if (nlhs != NLHS) {
		mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
		mexPrintf("%d inputs required.\n", NLHS);
	}
		K2d = (int)mxGetScalar(prhs[0]);
		Np2d = (int)mxGetScalar(prhs[1]);
}